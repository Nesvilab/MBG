package com.mbg;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanResidue;

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
            String siteID = psm.getGlycoSiteID();
            if (observedGlycoSites.containsKey(siteID)) {
                observedGlycoSites.get(siteID).addPSM(psm);
            } else {
                // new site
                observedGlycoSites.put(siteID, new GlycoSite(psm, siteID));
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
    public ArrayList<PSM> fillInGlycoSites(String rawFileName, String msFilePath) throws Exception {
        ArrayList<PSM> newPSMs = new ArrayList<>();
        GlycanResidue mod_nh3 = params.glycanResidues.get("NH4");
        GlycanResidue neuAc = params.glycanResidues.get("NeuGc");  // for generating decoys

        // read msFile
        MSMS ms2 = new MSMS(rawFileName, msFilePath);
        ms2.msFileReader();
        // if PASEF data, get im for each MS2
        if (!params.noPASEF) {
            Path mzBinPath = Paths.get(msFilePath);
            mzBinPath = mzBinPath.resolve(rawFileName+".mzBIN");
            File mzBinFile = mzBinPath.toFile();
            if (!mzBinFile.exists()) {
                System.out.printf("No mzBIN file found for %s\n", rawFileName);
                System.exit(-1);
            }
            MSMS.getMS2IMs(mzBinFile, 7);
        }
        // find file for IonQuant
        String foundFile = findFirstPriorityFile(msFilePath, rawFileName);
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
            if (site.psms.size() >= params.minPSMs && site.observedGlycans.size() >= params.minGlycans) {
                for (PSM base_PSM : site.basePSMs) {
                    int rankCounter = 1;

                    ArrayList<Glycan> dummy_glycans = PSM.getDummyGlycans(base_PSM, params.toAddResidues);
                    for (Glycan new_glycan: dummy_glycans) {
                        if (site.observedGlycans.containsKey(new_glycan.toString()) || site.basePSM == null) {
                            // this glycan was observed for this site, no action required
                            continue;
                        }
                        if (!Arrays.asList(qualifiedGlycans).contains(new_glycan.toString())) {
                            // this glycan not detected in psm.tsv
                            continue;
                        }

                        // glycan not observed. Make a PSM for it
                        List<PSM> addedPSMs = PSM.makeDummyPSM(base_PSM, new_glycan, psmFile.PSMcolumns, params, rankCounter, cRT_adjuster, false);
                        // decoy, not add into newPSMs
                        PSM.makeDummyPSM(base_PSM, new_glycan, psmFile.PSMcolumns, params, rankCounter, cRT_adjuster, true);

                        newPSMs.addAll(addedPSMs);
                        // update observed glycans
                        site.observedGlycans.put(new_glycan.toString(), new ArrayList<>());
                        rankCounter++;
                    }

                }
                // make NH3 glycan
                int rankCounter = 1;
//                for (PSM base_psm : site.basePSMs) {
//                    Glycan new_glycan = PSM.addResidual(base_psm.glycan, mod_nh3);
//                    List<PSM> addedPSMs = PSM.makeDummyPSM(base_psm, new_glycan, psmFile.PSMcolumns, params, rankCounter, cRT_adjuster);
//                    newPSMs.addAll(addedPSMs);
//                }
                // generate decoy for every psm
//                for (PSM base_psm : site.basePSMs) {
//                    Glycan new_glycan = PSM.addResidual(base_psm.glycan, neuAc);
//                    List<PSM> addedPSMs = PSM.makeDummyPSM(base_psm, new_glycan, psmFile.PSMcolumns, params, rankCounter, cRT_adjuster, true);
//                    newPSMs.addAll(addedPSMs);
//                }
            }
        }
        // get threshold score and coefficient of score function
        Scoring.getCutOffScore();
        // get qualified PSMs
        ArrayList<PSM> qualifiedPSMs = new ArrayList<>();
        for (PSM p: newPSMs) {
            if (p.intensity <0.01) {
                continue;  // remove unmatched PSMs
            }
            if (p.hasY1Ion) {  // if detected Y1 ion
                qualifiedPSMs.add(p);
            } else {
                if (Scoring.overThreshold(p, params.noPASEF)) {
                    qualifiedPSMs.add(p);
                }
            }
        }
        return qualifiedPSMs;
    }
}
