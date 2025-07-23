package com.mbg;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import umich.ms.glyco.Glycan;

import static com.mbg.Peak.findFirstPriorityFile;
import static java.lang.System.out;

/**
 * GlycanMatching class handles the core functionality of matching and inferring glycans
 * for peptide-spectrum matches (PSMs). It processes observed glycans and generates
 * potential new glycopeptide matches based on known glycan modifications.
 */
public class GlycanMatching {
    // Class fields
    public PSM_File psmFile;                    // Reference to the PSM file being processed
    public ArrayList<PSM> psms;                 // List of PSMs to process
    private final Params params;                // Configuration parameters
    private HashMap<String, GlycoSite> observedGlycoSites;  // Map of observed glycosylation sites
    private String[] qualifiedGlycans;          // Glycans that meet minimum observation criteria

    /**
     * Constructor for GlycanMatching
     * @param params Configuration parameters
     * @param psmFile PSM file reference
     * @param psms List of PSMs to process
     */
    public GlycanMatching(Params params, PSM_File psmFile, ArrayList<PSM> psms) {
        this.psmFile = psmFile;
        this.psms = psms;
        this.params = params;
    }

    /**
     * Groups PSMs by their glycosylation sites and counts observed glycans.
     * This method populates the observedGlycoSites map and identifies qualified glycans
     * that meet the minimum observation criteria.
     */
    public void groupGlycosites() {
        observedGlycoSites = new HashMap<>();
        HashMap<String, Integer> glycanObservationCounts = new HashMap<>();

        // Process each PSM to group by glycosite and count glycan observations
        for (PSM psm: psms) {
            if (psm.glycan == null) {
                continue;
            }
            
            // Group PSMs by peptide sequence
            if (observedGlycoSites.containsKey(psm.peptide)) {
                observedGlycoSites.get(psm.peptide).addPSM(psm);
            } else {
                observedGlycoSites.put(psm.peptide, new GlycoSite(psm, psm.peptide));
            }
            
            // Count occurrences of each glycan
            String glycanString = psm.glycan.toString();
            glycanObservationCounts.merge(glycanString, 1, Integer::sum);
        }

        // Filter to keep only glycans that meet the minimum observation criteria
        qualifiedGlycans = glycanObservationCounts.entrySet()
                .stream()
                .filter(entry -> entry.getValue() > 0)
                .map(Map.Entry::getKey)
                .toArray(String[]::new);
    }

    /**
     * Processes glycosites and generates new PSMs for unobserved glycans.
     * This method:
     * 1. Reads and processes MS/MS data
     * 2. Generates dummy PSMs for potential glycans
     * 3. Scores and filters the results
     * 
     * @param rawFileName Name of the raw file being processed
     * @return List of qualified PSMs that passed filtering
     * @throws Exception if there are issues processing the MS data
     */
    public ArrayList<PSM> fillInGlycoSites(String rawFileName) throws Exception {
        ArrayList<PSM> inferredPSMs = new ArrayList<>();
        String msFilePath = matchMSfilename(rawFileName, params.msFilePathsMap);
        Path msFileDir = Paths.get(msFilePath).getParent();

        // Process MS/MS data
        out.printf("Reading spectrum file %s\n", rawFileName);
        MSMS ms2 = new MSMS(rawFileName, msFilePath);
        ms2.msFileReader(params.numThreads);

        // Handle PASEF data if available
        if (!params.noPASEF) {
            Path mzBinPath = msFileDir.resolve(rawFileName + ".mzBIN");
            File mzBinFile = mzBinPath.toFile();
            if (!mzBinFile.exists()) {
                System.out.printf("No mzBIN file found for %s\n", rawFileName);
                System.exit(-1);
            }
            MSMS.getMS2IMs(mzBinFile, params.numThreads);
        }

        // Prepare for IonQuant analysis
        String[] supportedExtensions = {".mzML", ".raw", ".d"};
        String foundFile = findFirstPriorityFile(msFileDir.toString(), rawFileName, supportedExtensions);
        if (foundFile != null) {
            Peak.indexBuilder(foundFile, params);
        } else {
            System.out.println("No file found with extension .mzML, .raw or .d.");
            System.exit(-1);
        }

        // Initialize RT adjustment and scoring
        RT_adjuster rtAdjuster = new RT_adjuster(this.observedGlycoSites);
        Scoring.init();

        // Process each glycosite
        for (Map.Entry<String, GlycoSite> siteEntry: observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            if (site.psms.size() >= params.minPSMs && site.observedGlycans.size() >= params.minGlycans) {
                for (PSM basePSM : site.basePSMs) {
                    getInferredPSMs(basePSM, params, site, rtAdjuster, inferredPSMs, params.expandDB, params.maxSkips);
                }
            }
        }

        // Filter and score PSMs
        ArrayList<PSM> qualifiedPSMs = new ArrayList<>();
        List<Double> glycanMassList = new ArrayList<>();

        if (Scoring.decoyData.isEmpty() || Scoring.decoyData.size() < 2 || Scoring.targetData.size() < 2) {
            // Simple filtering when scoring is not possible
            for (PSM psm: inferredPSMs) {
                if (psm.intensity < 0.01) continue;
                if (!glycanMassList.contains(psm.glycan.mass)) {
                    glycanMassList.add(psm.glycan.mass);
                }
                qualifiedPSMs.add(psm);
            }
        } else {
            // Full scoring and filtering
            if (params.noPASEF) {
                System.out.println("Coefficients Items: [Intensity, mass error, kl, Shift Frequency, Relative Intensity of Y0+Y1, RT Difference]");
            } else {
                System.out.println("Coefficients Items: [Intensity, mass error, kl, Shift Frequency, Relative Intensity of Y0+Y1, RT Difference, IM Difference]");
            }

            Scoring.getCutOffScore(params.fdr);
            for (PSM psm: inferredPSMs) {
                if (psm.intensity < 0.01) continue;
                double score = Scoring.getMBGScore(psm, psmFile.PSMcolumns, params.noPASEF);
                if (Scoring.overThreshold(score)) {
                    if (!glycanMassList.contains(psm.glycan.mass)) {
                        glycanMassList.add(psm.glycan.mass);
                    }
                    qualifiedPSMs.add(psm);
                }
            }
        }

        // Write glycan masses to file
        if (!glycanMassList.isEmpty()) {
            PSM_File.massWriter(glycanMassList);
        }
        
        // Additional filtering for TMT-labeled samples
        if (params.runTMT) {
            ArrayList<PSM> qualifiedPSMsTMT = new ArrayList<>();
            for (PSM p: qualifiedPSMs) {
                if (p.yIntensities > 0) {  // if detected Y1 ion
                    qualifiedPSMsTMT.add(p);
                }
            }
            System.out.printf("Filtered %d PSMs without MS2 hit, remained %d PSMs %n", qualifiedPSMs.size() - qualifiedPSMsTMT.size(), qualifiedPSMsTMT.size());
            return qualifiedPSMsTMT;
        }

        return qualifiedPSMs;
    }

    /**
     * Recursively generates inferred PSMs for a given base PSM and glycosite.
     * 
     * @param basePSM The base PSM to start from
     * @param params Configuration parameters
     * @param site The glycosite being processed
     * @param rtAdjuster RT adjustment calculator
     * @param newPSMs List to store new PSMs
     * @param steps Number of remaining expansion steps
     */
    private void getInferredPSMs(PSM basePSM, Params params, GlycoSite site, RT_adjuster rtAdjuster, 
                               ArrayList<PSM> newPSMs, Integer steps, Integer remainMissingSteps) {
        if (steps < 0) return;

        int rankCounter = basePSM.maxRank; // add rank to spectrum file to avoid bug in IonQuant and DDA+ (if >1 PSM per spectrum)
        HashMap<Glycan, String> potentialGlycans = PSM.getDummyGlycansV2(basePSM, params.toAddGlycans);

        for (Map.Entry<Glycan, String> glycanEntry : potentialGlycans.entrySet()) {
            Glycan newGlycan = glycanEntry.getKey();
            String addedResidue = glycanEntry.getValue();

            // Skip if glycan already observed or no base PSM
            if (site.observedGlycans.containsKey(newGlycan.toString()) || site.basePSM == null) {
                continue;
            }

            boolean isNewGlycan = !Arrays.asList(qualifiedGlycans).contains(newGlycan.toString());
            if (isNewGlycan && params.expandDB <= 0) {
                continue;
            }

            // Skip if mass is too close to existing glycans
            if (isWithinDistance(site.glycanMasses, newGlycan.mass, 1.5)) {
                continue;
            }

            // Generate new PSMs
            List<PSM> addedPSMs = PSM.makeDummyPSM(basePSM, newGlycan, psmFile.PSMcolumns, params,
                    rankCounter, rtAdjuster, false, basePSM.glycan.toString(),
                    addedResidue, isNewGlycan);

            // Generate decoy PSMs
            PSM.makeDummyPSM(basePSM, newGlycan, psmFile.PSMcolumns, params, rankCounter, rtAdjuster,
                    true, basePSM.glycan.toString(), addedResidue, isNewGlycan);

            newPSMs.addAll(addedPSMs);

            // Update observed glycans if match found
            if (addedPSMs.get(0).intensity > 0.01) {
                site.observedGlycans.put(newGlycan.toString(), new ArrayList<>());
                site.glycanMasses.add(newGlycan.mass);
                rankCounter++; // increase rankCounter only when a peak is found
                basePSM.maxRank++; // increment max rank for the next PSM
            }

            // Recursively process if within expansion limit and match found
            if (steps > 0) {
                if (addedPSMs.get(0).intensity < 0.001) {
                    if (remainMissingSteps > 0) {
                        // continue inferring even if missing iff remainMissingSteps > 0
                        getInferredPSMs(addedPSMs.get(0), params, site, rtAdjuster, newPSMs, steps - 1, remainMissingSteps - 1);
                    }
                } else {
                    // continue inferring if a peak is found
                    getInferredPSMs(addedPSMs.get(0), params, site, rtAdjuster, newPSMs, steps - 1, remainMissingSteps);
                }
            }
        }
    }

    /**
     * Checks if a value is within a specified distance of any element in a list.
     * 
     * @param list List of values to check against
     * @param inputValue Value to check
     * @param distance Maximum allowed distance
     * @return true if inputValue is within distance of any element in list
     */
    private boolean isWithinDistance(List<Double> list, double inputValue, double distance) {
        return list.stream().anyMatch(element -> Math.abs(element - inputValue) <= distance);
    }

    /**
     * Matches a raw file name to its corresponding MS file path from the manifest.
     * 
     * @param rawFileName Name of the raw file
     * @param msFilePathsMap Map of raw file names to their paths
     * @return The matched file path
     */
    private String matchMSfilename(String rawFileName, HashMap<String, String> msFilePathsMap) {
        if (msFilePathsMap.containsKey(rawFileName)) {
            return msFilePathsMap.get(rawFileName);
        } else {
            out.printf("Error: Could not find MS spectrum file for rawfile %s from the psm.tsv. Please check the manifest file and psm.ts.\n", rawFileName);
            System.exit(1);
            return null;
        }
    }
}
