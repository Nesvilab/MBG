package com.mbg;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.IntStream;

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
    public HashSet<String> newGlycans = new HashSet<>();    // List of new inferred glycans (not in the input database)
    public HashSet<Integer> excludeScanNumbers;
    private String currentRawFileName;
    private int type1Count;
    private int type2Count;
    /** Decoy PSMs that found an MS1 hit, collected for diagnostic output after LDA fitting. */
    private ArrayList<PSM> decoyPSMsForDiagnostics = new ArrayList<>();

    /**
     * Constructor for GlycanMatching
     * @param params Configuration parameters
     * @param psmFile PSM file reference
     * @param psms List of PSMs to process
     */
    public GlycanMatching(Params params, PSM_File psmFile, ArrayList<PSM> psms, String fileName) {
        this.psmFile = psmFile;
        this.psms = psms;
        this.params = params;
        // if not allow chimeric spectra
        if (!this.params.allowChimeric) {
            this.excludeScanNumbers = params.excludeScanNumbersMAp.get(fileName);
        } else {
            this.excludeScanNumbers = new HashSet<>();
        }
    }

    /**
     * Groups PSMs by their glycosylation sites and counts observed glycans.
     * This method populates the observedGlycoSites map and identifies qualified glycans
     * that meet the minimum observation criteria.
     */
    public void groupGlycosites() {
        observedGlycoSites = new HashMap<>();
        HashMap<String, Integer> glycanObservationCounts = new HashMap<>();

        // Process each PSM to group by glycosite and count glycan observations.
        // Inference is seeded only by confidently identified glycopeptides; above-maxq
        // rows remain in allPSMs so that the four-class FDR lookups can see them.
        for (PSM psm: psms) {
            if (psm.glycanState != PSM.GlycanState.CONFIDENT) {
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
     * @param rawFileName        Name of the raw file being processed
     * @param combinedShiftedRTs
     * @param combinedShiftedIMs
     * @return List of qualified PSMs that passed filtering
     * @throws Exception if there are issues processing the MS data
     */
    public ArrayList<PSM> fillInGlycoSites(String rawFileName, HashMap<String, List<Float>> combinedShiftedRTs, HashMap<String, List<Float>> combinedShiftedIMs) throws Exception {
        ArrayList<PSM> inferredPSMs = new ArrayList<>();
        this.currentRawFileName = rawFileName;
        this.type1Count = 0;
        this.type2Count = 0;
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

        // Initialize RT adjustment and scoring; reset per-raw-file decoy collection.
        RT_adjuster rtAdjuster = new RT_adjuster(this.observedGlycoSites, params.minPSMs, combinedShiftedRTs, combinedShiftedIMs);
        Scoring.init();
        decoyPSMsForDiagnostics = new ArrayList<>();

        // Process each glycosite
        for (Map.Entry<String, GlycoSite> siteEntry: observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            if (site.observedGlycans.size() >= params.minGlycans) {
                for (PSM basePSM : site.basePSMs) {
                    getInferredPSMs(basePSM, params, site, rtAdjuster, inferredPSMs, params.expandDB, params.maxSkips);
                }
            }
        }

        System.out.printf("Candidate types for %s: Type 1 (novel) = %d, Type 2 (already identified) = %d%n",
                rawFileName, type1Count, type2Count);

        // Filter and score PSMs
        ArrayList<PSM> qualifiedPSMs = new ArrayList<>();
        List<Double> glycanMassList = new ArrayList<>();

        int totalTargets = Scoring.type2Targets.size() + Scoring.type1Targets.size();
        int totalDecoys  = Scoring.type2Decoys.size()  + Scoring.type1Decoys.size();
        if (totalTargets < 2 || totalDecoys < 2) {
            // Not enough training data — fall back to intensity-only pass-through.
            for (PSM psm : inferredPSMs) {
                if (psm.intensity < 0.01) continue;
                if (!glycanMassList.contains(psm.glycan.mass)) glycanMassList.add(psm.glycan.mass);
                qualifiedPSMs.add(psm);
            }
        } else {
            if (params.noPASEF) {
                System.out.println("Coefficients Items: [Intensity, mass error, kl, Shift Frequency, Relative Intensity of Y0+Y1, RT Difference]");
            } else {
                System.out.println("Coefficients Items: [Intensity, mass error, kl, Shift Frequency, Relative Intensity of Y0+Y1, RT Difference, IM Difference]");
            }

            Scoring.getCutOffScore(params.qvalue);

            // Record decoy PSMs that found an MS1 hit (collected during candidate generation).
            for (PSM decoy : decoyPSMsForDiagnostics) {
                double decoyScore = Scoring.getMBGScore(decoy, psmFile.PSMcolumns, params.noPASEF);
                Diagnostics.addDecoyRecord(decoy, decoyScore, params.noPASEF);
            }

            // Score every PSM; route TYPE_2 to always-accept, accumulate TYPE_1 for EM.
            List<PSM>    type1PSMs    = new ArrayList<>();
            List<Double> type1Scores  = new ArrayList<>();

            for (PSM psm : inferredPSMs) {
                if (psm.intensity < 0.01) continue;
                double score = Scoring.getMBGScore(psm, psmFile.PSMcolumns, params.noPASEF);

                if (psm.candidateType == PSM.CandidateType.TYPE_2) {
                    Diagnostics.addTargetRecord(psm, score, Double.NaN, Double.NaN, params.noPASEF);
                    if (!glycanMassList.contains(psm.glycan.mass)) glycanMassList.add(psm.glycan.mass);
                    qualifiedPSMs.add(psm);
                } else {
                    type1PSMs.add(psm);
                    type1Scores.add(score);
                }
            }

            // TYPE_1 accept: mixture-model EM q-value, or legacy LDA threshold.
            List<Double>[] emScores = Scoring.scoreType1Candidates();
            int nType1 = emScores[0].size() + emScores[1].size();

            if (!params.legacyFdr && nType1 >= MixtureModelEM.N_MIN) {
                MixtureModelEM em = new MixtureModelEM(emScores[0], emScores[1]);
                List<Double> fitData = new ArrayList<>(emScores[0]);
                fitData.addAll(emScores[1]);
                em.fit(fitData);

                // Sort indices by LDA score descending for q-value sweep
                Integer[] order = IntStream.range(0, type1PSMs.size()).boxed()
                        .sorted((a, b) -> Double.compare(type1Scores.get(b), type1Scores.get(a)))
                        .toArray(Integer[]::new);

                double[] fdrs = new double[order.length];
                double[] post = new double[order.length];
                double fpSum = 0;
                for (int rank = 0; rank < order.length; rank++) {
                    double p = em.posterior(type1Scores.get(order[rank]));
                    post[rank] = p;
                    fpSum += (1.0 - p);
                    fdrs[rank] = fpSum / (rank + 1);
                }

                // Monotonize q-values: q_i = min_{k >= i} FDR_k
                double[] qv = new double[order.length];
                double minFDR = fdrs[order.length - 1];
                for (int rank = order.length - 1; rank >= 0; rank--) {
                    minFDR = Math.min(minFDR, fdrs[rank]);
                    qv[rank] = minFDR;
                }

                // Map back to PSM index space
                double[] qByPSM   = new double[type1PSMs.size()];
                double[] postByPSM = new double[type1PSMs.size()];
                for (int rank = 0; rank < order.length; rank++) {
                    qByPSM[order[rank]]    = qv[rank];
                    postByPSM[order[rank]] = post[rank];
                }

                int accepted = 0;
                for (int i = 0; i < type1PSMs.size(); i++) {
                    PSM psm = type1PSMs.get(i);
                    Diagnostics.addTargetRecord(psm, type1Scores.get(i), postByPSM[i], qByPSM[i], params.noPASEF);
                    if (qByPSM[i] <= params.qvalue) {
                        if (!glycanMassList.contains(psm.glycan.mass)) glycanMassList.add(psm.glycan.mass);
                        qualifiedPSMs.add(psm);
                        accepted++;
                    }
                }
                System.out.printf("TYPE_1 EM: %d/%d accepted at q <= %.3f%n",
                        accepted, type1PSMs.size(), params.qvalue);

            } else {
                // Legacy or insufficient data: gate on LDA threshold
                if (params.legacyFdr) {
                    System.out.println("Legacy LDA threshold mode (--legacy-fdr).");
                } else {
                    System.out.printf("TYPE_1 sample size %d < %d; using LDA threshold.%n",
                            nType1, MixtureModelEM.N_MIN);
                }
                for (int i = 0; i < type1PSMs.size(); i++) {
                    PSM psm = type1PSMs.get(i);
                    Diagnostics.addTargetRecord(psm, type1Scores.get(i), Double.NaN, Double.NaN, params.noPASEF);
                    if (Scoring.overThreshold(type1Scores.get(i))) {
                        if (!glycanMassList.contains(psm.glycan.mass)) glycanMassList.add(psm.glycan.mass);
                        qualifiedPSMs.add(psm);
                    }
                }
            }
        }

        // Write glycan masses to file
        if (!glycanMassList.isEmpty()) {
            Utils.massWriter(glycanMassList);
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
        if (steps < 1) return;

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

            // Skip if mass is too close to existing glycans (preserves pre-Commit-2 behavior;
            // the 1.5 Da deduplication will move to output-time in Commit 6).
            if (isWithinDistance(site.glycanMasses, newGlycan.mass, 1.5)) {
                continue;
            }

            // Classify as Type 2 if (peptide, glycan) was confidently identified in ANY
            // raw file of this experiment — a replicate hit is a known true positive even
            // when absent from the current file's observedGlycans (e.g., missing value).
            // Uses the cross-file global index. Plumbing only — Scoring still uses the
            // merged targetData / decoyData until the four-list refactor in Commit 5.
            PSM.CandidateType candidateType = psmFile.isConfidentlyIdentified(
                    basePSM.peptide, newGlycan.toString())
                    ? PSM.CandidateType.TYPE_2
                    : PSM.CandidateType.TYPE_1;
            if (candidateType == PSM.CandidateType.TYPE_2) {
                type2Count++;
            } else {
                type1Count++;
            }

            // Generate target PSMs
            List<PSM> addedPSMs = PSM.makeDummyPSM(basePSM, newGlycan, psmFile.PSMcolumns, params,
                    rankCounter, rtAdjuster, false, basePSM.glycan.toString(),
                    addedResidue, isNewGlycan, excludeScanNumbers, candidateType, psmFile, currentRawFileName);

            // Generate decoy PSMs; collect those with an MS1 hit for diagnostic output.
            List<PSM> decoyList = PSM.makeDummyPSM(basePSM, newGlycan, psmFile.PSMcolumns, params, rankCounter, rtAdjuster,
                    true, basePSM.glycan.toString(), addedResidue, isNewGlycan, excludeScanNumbers, candidateType, psmFile, currentRawFileName);
            if (!decoyList.isEmpty() && decoyList.get(0).intensity > 0.01) {
                decoyPSMsForDiagnostics.add(decoyList.get(0));
            }

            newPSMs.addAll(addedPSMs);

            // Update observed glycans if match found
            if (addedPSMs.get(0).intensity > 0.01) {
                site.observedGlycans.put(newGlycan.toString(), new ArrayList<>());
                site.glycanMasses.add(newGlycan.mass);
                rankCounter++; // increase rankCounter only when a peak is found
                basePSM.maxRank++; // increment max rank for the next PSM
                if (isNewGlycan) {
                    // found new glycan: add to list for printing to downstream tools
                    newGlycans.add(newGlycan.toString());
                }
            }

            // Recursively process if within expansion limit and match found
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
