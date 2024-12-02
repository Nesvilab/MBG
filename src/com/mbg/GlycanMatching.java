package com.mbg;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import umich.ms.glyco.Glycan;

import static com.mbg.Peak.findFirstPriorityFile;
import static java.lang.System.out;

public class GlycanMatching {
    public PSM_File psmFile;
    public ArrayList<PSM> psms;
    private final Params params;
    private HashMap<String, GlycoSite> observedGlycoSites;
    private String[] qualifiedGlycans;  // glycans in psm.tsv with a specific minimum count

    public GlycanMatching(Params params, PSM_File psmFile, ArrayList<PSM> psms) {
        this.psmFile = psmFile;
        this.psms = psms;
        this.params = params;
    }

    public void groupGlycosites() {
        observedGlycoSites = new HashMap<>();
        HashMap<String, Integer> allObservedGlycans = new HashMap<>();

        for (PSM psm: psms) {
            if (psm.glycan == null) {
                continue;
            }
            if (observedGlycoSites.containsKey(psm.peptide)) {
                observedGlycoSites.get(psm.peptide).addPSM(psm);
            } else {
                // new site
                observedGlycoSites.put(psm.peptide, new GlycoSite(psm, psm.peptide));
            }
            // counts observed glycans
            if (allObservedGlycans.containsKey(psm.glycan.toString())) {
                allObservedGlycans.put(psm.glycan.toString(), allObservedGlycans.get(psm.glycan.toString()) + 1);
            } else {
                allObservedGlycans.put(psm.glycan.toString(), 1);
            }
        }
        // filter to keep entries where the glycan psm counts >5
        qualifiedGlycans = allObservedGlycans.entrySet()
                .stream()
                .filter(entry -> entry.getValue() > 5)
                .map(Map.Entry::getKey).toArray(String[]::new);
    }

    /**
     * Given a glycan database and rules from Params, fill in the unobserved glycans for each glycosite
     * and generate dummy PSMs for them.
     */
    public ArrayList<PSM> fillInGlycoSites(String rawFileName) throws Exception {
        ArrayList<PSM> newPSMs = new ArrayList<>();
        String msFilePath = matchMSfilename(rawFileName, params.msFilePathsMap);
//        msFilePath = "Z:\\dpolasky\\data\\glyco\\pGlyco2\\yeast_PXD005565\\cwq_mix2-1_726.raw";
        Path msFileDir = Paths.get(msFilePath).getParent();

        // read msFile
        out.printf("Reading spectrum file %s\n", rawFileName);
        MSMS ms2 = new MSMS(rawFileName, msFilePath);
        ms2.msFileReader(params.numThreads);
        // if PASEF data, get im for each MS2
        if (!params.noPASEF) {
            Path mzBinPath = msFileDir.resolve(rawFileName+".mzBIN");
            File mzBinFile = mzBinPath.toFile();
            if (!mzBinFile.exists()) {
                System.out.printf("No mzBIN file found for %s\n", rawFileName);
                System.exit(-1);
            }
            MSMS.getMS2IMs(mzBinFile, params.numThreads);
        }
        // find file for IonQuant
        String[] extensions = {".mzML", ".raw", ".d"};
        String foundFile = findFirstPriorityFile(msFileDir.toString(), rawFileName, extensions);
        if (foundFile != null) {
            // build index for IonQuant
            Peak.indexBuilder(foundFile, params);
        } else {
            System.out.println("No file found with extension .mzML, .raw or .d.");
            System.exit(-1);
        }
        // get apex rt for each psm and generate internal median RT shift database; get base PSMs
        RT_adjuster cRT_adjuster = new RT_adjuster(this.observedGlycoSites);
        // init decoys and targets
        Scoring.init();

        for (Map.Entry<String, GlycoSite> siteEntry: observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            if (Objects.equals(site.id, "TKPREEQYNSTYR")) {
                System.out.println(1);
            }
            if (site.psms.size() >= params.minPSMs && site.observedGlycans.size() >= params.minGlycans) {
                for (PSM base_PSM : site.basePSMs) {
                    getInferredPSMs(base_PSM, params, site, cRT_adjuster, newPSMs, params.expandDB);
                }
            }
        }
        ArrayList<PSM> qualifiedPSMs = new ArrayList<>();
        // if decoy list is empty or decoy or target list is too small to scoring, return all newPSMs
        if (Scoring.decoyData.isEmpty() || Scoring.decoyData.size() < 2 || Scoring.targetData.size() < 2) {
            for (PSM p: newPSMs) {
                if (p.intensity <0.01) {
                    continue;  // remove unmatched PSMs
                }
                qualifiedPSMs.add(p);
            }
            return qualifiedPSMs;
        }
        // if TMT labeled sample, only keep PSMs with MS2
        if (params.runTMT) {
            for (PSM p: newPSMs) {
                if (p.yIntensities > 0) {  // if detected Y1 ion
                    qualifiedPSMs.add(p);
                }
            }
            return qualifiedPSMs;
        }
        // get threshold score and coefficient of score function
        Scoring.getCutOffScore(params.fdr);
        // get qualified PSMs
        for (PSM p: newPSMs) {
            if (p.intensity <0.01) {
                continue;  // remove unmatched PSMs
            }
            // get updated score
            double cScore = Scoring.getMBGScore(p, psmFile.PSMcolumns, params.noPASEF);
            if (Scoring.overThreshold(cScore)) {
                qualifiedPSMs.add(p);
            }
        }
        return qualifiedPSMs;
    }

    private void getInferredPSMs(PSM base_PSM, Params params, GlycoSite site, RT_adjuster cRT_adjuster, ArrayList<PSM> newPSMs, Integer steps) {
        if (steps < 0) {
            return;
        }
        int rankCounter = 1;
        HashMap<Glycan, String> dummy_glycans = PSM.getDummyGlycansV2(base_PSM, params.toAddGlycans);
        for (Map.Entry<Glycan, String> glycanEntry: dummy_glycans.entrySet()) {
            Glycan new_glycan = glycanEntry.getKey();
            String addedResidual = glycanEntry.getValue();
            if (site.observedGlycans.containsKey(new_glycan.toString()) || site.basePSM == null) {
                // this glycan was observed for this site, no action required
                continue;
            }
            boolean expanded_glycan = false;
            if (!Arrays.asList(qualifiedGlycans).contains(new_glycan.toString())) {
                // this glycan not detected in psm.tsv ever, and do not expanding glycan database
                if (params.expandDB > 0) {
                    expanded_glycan = true;
                } else {
                    continue;
                }
            }
            // check if mass of dummy glycan is too close to glycans on current site
            if (isWithinDistance(site.glycanMasses, new_glycan.mass, 1.5)) {
                continue;
            }
            // glycan not observed. Make a PSM for it
            List<PSM> addedPSMs = PSM.makeDummyPSM(base_PSM, new_glycan, psmFile.PSMcolumns, params, rankCounter, cRT_adjuster, false, base_PSM.glycan.toString(), addedResidual);
            // decoy, not add into newPSMs
            PSM.makeDummyPSM(base_PSM, new_glycan, psmFile.PSMcolumns, params, rankCounter, cRT_adjuster, true, base_PSM.glycan.toString(), addedResidual);

            newPSMs.addAll(addedPSMs);
            // update observed glycans
            if (addedPSMs.get(0).intensity > 0.01) {
                site.observedGlycans.put(new_glycan.toString(), new ArrayList<>());
                site.glycanMasses.add(new_glycan.mass);
            }
            rankCounter++;
            // set current inferred PSMs as base PSM, if within path and has matched PSM
            if (steps>0 & addedPSMs.get(0).intensity > 0.01) {
                getInferredPSMs(addedPSMs.get(0), params, site, cRT_adjuster, newPSMs, steps-1);
            }
        }
    }

    private boolean isWithinDistance(List<Double> list, double inputValue, double distance) {
        for (double element : list) {
            if (Math.abs(element - inputValue) <= distance) {
                return true;
            }
        }
        return false;
    }

    // Match the name in the psm.tsv (base filename without the extension) to the file path from the manifest
    private String matchMSfilename(String rawfileName, HashMap<String, String> msFilePathsMap) {
        if (msFilePathsMap.containsKey(rawfileName)) {
            return msFilePathsMap.get(rawfileName);
        } else {
            out.printf("Error: Could not find MS spectrum file for rawfile %s from the psm.tsv. Please check the manifest file and psm.tsv\n", rawfileName);
            System.exit(1);
            return null;
        }
    }
}
