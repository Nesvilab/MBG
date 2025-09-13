package com.mbg;

import com.mbg.msreader.Spectrum;
import ionquant.api.Entry;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanParser;
import umich.ms.glyco.GlycanResidue;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PSM implements Comparable<PSM>{

    double rt;
    double im;
    double apexRT;
    double apexIM;
    float averageMass;
    int scanNum;
    String peptide;
    Glycan glycan;
    String protein;
    int proteinGlycositeIndex;
    int peptideGlycositeIndex;
    ArrayList<Integer> peptideGlycositeIndexes;
    String[] splits;
    char glycoSiteAA;
    double calcMass;
    double calcMz;
    int charge;
    double glycoQval;
    double deltaMass;
    boolean hasMS2Spectrum;
    public double yIntensities;
    float aveMassDiff;
    float imDiff;
    float rtDiff;
    int rank;
    double intensity;
    double avePPM;
    double kl; // isotope distribution score
    int shiftFrequency;
    String confidenceLevel;
    double featureIntensity;
    public String rawfile;
    public int maxRank;     // if there are multiple PSMs for a given scan (DDA+, DIA), the max rank found in the PSM table

    private static final Pattern modMassPattern = Pattern.compile("[\\d]+[A-Z]\\(([\\d.]+)\\)");
    private static final Pattern modSitePattern = Pattern.compile("([\\d]+)[A-Z]\\([\\d.]+\\)");
    private static final Pattern rankPattern = Pattern.compile("_rank([\\d]+)\\.pep\\.xml");


    /**
     * Initialize a PSM object from a line of the psm.tsv table
     */
    public PSM(String psmLine, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params, boolean isInferred, boolean hasMS2, double hasY1, String parentGlycan, String addedResidual, int addCols) {
        splits = psmLine.split("\t");

        // expand to hold new columns if needed
        if (addCols > 0) {
            int add = addCols;
            ArrayList<String> newSplits = new ArrayList<>(Arrays.asList(splits));
            while (add > 0) {
                newSplits.add(PSM_File.insertIndex, "");
                add--;
            }
            splits = newSplits.toArray(new String[0]);
            hasMS2Spectrum = hasMS2;
            yIntensities = hasY1;
            confidenceLevel = "";
        } else {
            if (!isInferred) {
                // if this is a re-analysis (i.e., not adding new columns) and this is not an inferred PSM from this analysis, read its inferred status from the previous result
                hasMS2Spectrum = Boolean.parseBoolean(splits[psmCols.get(PSM_File.PSMattr.has_ms2_spec)]);
                yIntensities = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.has_Y1)]);
                confidenceLevel = "";
            } else {
                hasMS2Spectrum = hasMS2;
                yIntensities = hasY1;
                confidenceLevel = "";
            }
        }

        // update new fields and information that is needed for all PSMs (including non-glyco)
        splits[psmCols.get(PSM_File.PSMattr.has_ms2_spec)] = String.valueOf(hasMS2Spectrum);
        splits[psmCols.get(PSM_File.PSMattr.has_Y1)] = String.valueOf(yIntensities);
        splits[psmCols.get(PSM_File.PSMattr.MBG_score)] = String.valueOf(0);
        splits[psmCols.get(PSM_File.PSMattr.parent_glycan)] = parentGlycan;
        splits[psmCols.get(PSM_File.PSMattr.added_residual)] = addedResidual;
        splits[psmCols.get(PSM_File.PSMattr.confidence_level)] = confidenceLevel;
        Matcher rankMatch = rankPattern.matcher(splits[psmCols.get(PSM_File.PSMattr.spectrum_file)]);
        if (rankMatch.find()) {
            rank = Integer.parseInt(rankMatch.group(1));
        } else {
            // no ranks provided, so need to make the base PSM rank 1 so added PSMs can be other ranks
            updateRank(1, psmCols);
        }
        scanNum = getScanNum(psmCols);
        rawfile = getRawfile(psmCols);

        // parse glycan and scores
        String glycan = splits[psmCols.get(PSM_File.PSMattr.glycan)];
        if (glycan.isEmpty()) {
            // non-glyco PSM
            this.glycan = null;
            return;
        }
        if (!psmCols.containsKey(PSM_File.PSMattr.glyco_q_val)) {
            glycoQval = 0;
        }
        else {
            try {
                glycoQval = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.glyco_q_val)]);
            } catch (NumberFormatException ex) {
                glycoQval = 1;
            }
        }
        if (glycoQval > params.maxGlycanQval && addCols > 0) {
            this.glycan = null;
            return;
        }
        this.glycan = GlycanParser.parseGlycanString(glycan, params.glycanResidues);
        // parse other needed PSM info
        peptide = splits[psmCols.get(PSM_File.PSMattr.peptide)];
        protein = splits[psmCols.get(PSM_File.PSMattr.protein)];
        if (!params.noPASEF) {
            im = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.IM)]);
        }
        rt = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.RT)]);
        int proteinStart = Integer.parseInt(splits[psmCols.get(PSM_File.PSMattr.protein_start)]);
        String assignedMods = splits[psmCols.get(PSM_File.PSMattr.assigned_mods)];
        peptideGlycositeIndexes = getGlycoSiteIndices(assignedMods);
        peptideGlycositeIndex = peptideGlycositeIndexes.get(0);
        glycoSiteAA = peptide.charAt(peptideGlycositeIndex - 1);
        proteinGlycositeIndex = proteinStart + peptideGlycositeIndex;
        calcMass = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.Calculated_MASS)]);
        calcMz = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.Calculated_MZ)]);
        charge = Integer.parseInt(splits[psmCols.get(PSM_File.PSMattr.charge)]);
        deltaMass = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.delta_mass)]);
        intensity = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.intensity)]);
        maxRank = rank;
    }

    public PSM(String psmLine, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params, int addCols) {
        splits = psmLine.split("\t");

        // expand to hold new columns if needed
        if (addCols > 0) {
            int add = addCols;
            ArrayList<String> newSplits = new ArrayList<>(Arrays.asList(splits));
            while (add > 0) {
                newSplits.add(PSM_File.insertIndex, "");
                add--;
            }
            splits = newSplits.toArray(new String[0]);
        }

        Matcher rankMatch = rankPattern.matcher(splits[psmCols.get(PSM_File.PSMattr.spectrum_file)]);
        if (rankMatch.find()) {
            rank = Integer.parseInt(rankMatch.group(1));
        } else {
            // no ranks provided, so need to make the base PSM rank 1 so added PSMs can be other ranks
            updateRank(1, psmCols);
        }
        scanNum = getScanNum(psmCols);
        rawfile = getRawfile(psmCols);

        // parse glycan and scores
        String glycan = splits[psmCols.get(PSM_File.PSMattr.glycan)];
        if (glycan.isEmpty()) {
            // non-glyco PSM
            this.glycan = null;
            return;
        }
        if (!psmCols.containsKey(PSM_File.PSMattr.glyco_q_val)) {
            glycoQval = 0;
        }
        else {
            try {
                glycoQval = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.glyco_q_val)]);
            } catch (NumberFormatException ex) {
                glycoQval = 1;
            }
        }
        if (glycoQval > params.maxGlycanQval) {
            this.glycan = null;
            return;
        }
        this.glycan = GlycanParser.parseGlycanString(glycan, params.glycanResidues);
        // parse other needed PSM info
        peptide = splits[psmCols.get(PSM_File.PSMattr.peptide)];
        protein = splits[psmCols.get(PSM_File.PSMattr.protein)];
        if (!params.noPASEF) {
            im = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.IM)]);
        }
        rt = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.RT)]);
        int proteinStart = Integer.parseInt(splits[psmCols.get(PSM_File.PSMattr.protein_start)]);
        String assignedMods = splits[psmCols.get(PSM_File.PSMattr.assigned_mods)];
        peptideGlycositeIndexes = getGlycoSiteIndices(assignedMods);
        peptideGlycositeIndex = peptideGlycositeIndexes.get(0);
        glycoSiteAA = peptide.charAt(peptideGlycositeIndex - 1);
        proteinGlycositeIndex = proteinStart + peptideGlycositeIndex;
        calcMass = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.Calculated_MASS)]);
        calcMz = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.Calculated_MZ)]);
        charge = Integer.parseInt(splits[psmCols.get(PSM_File.PSMattr.charge)]);
        deltaMass = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.delta_mass)]);
        intensity = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.intensity)]);
        maxRank = rank;
    }

    /**
     * Update the PSM's (internal) rank and spectrum file column (psm table)
     */
    private void updateRank(int newRank, HashMap<PSM_File.PSMattr, Integer> psmCols) {
        rank = newRank;
        String oldSpecFile = splits[psmCols.get(PSM_File.PSMattr.spectrum_file)].replaceAll("_rank[\\d]+\\.pep\\.xml", ".pep.xml").split("\\.pep\\.xml")[0];
        splits[psmCols.get(PSM_File.PSMattr.spectrum_file)] = String.format("%s_rank%d.pep.xml", oldSpecFile, rank);
    }

    /**
     * Update the PSM's (internal) rank and spectrum file column (psm table)
     */
    private void updateSplits(HashMap<PSM_File.PSMattr, Integer> psmCols, boolean noPASEF, boolean if_new_glycan) {
        splits[psmCols.get(PSM_File.PSMattr.mass_error)] = String.valueOf(avePPM);
        splits[psmCols.get(PSM_File.PSMattr.KL)] = String.valueOf(kl);
        splits[psmCols.get(PSM_File.PSMattr.shift_frequency)] = String.valueOf(shiftFrequency);
        splits[psmCols.get(PSM_File.PSMattr.rt_diff)] = String.valueOf(rtDiff);
        if (!noPASEF) {
            splits[psmCols.get(PSM_File.PSMattr.im_diff)] = String.valueOf(imDiff);
        }
        splits[psmCols.get(PSM_File.PSMattr.if_new_glycan)] = String.valueOf(if_new_glycan);
        if (featureIntensity > 0) {
            confidenceLevel = "3";
        } else if (hasMS2Spectrum) {
            confidenceLevel = "2";
        } else {
            confidenceLevel = "1";
        }
        splits[psmCols.get(PSM_File.PSMattr.confidence_level)] = confidenceLevel;
//        splits[psmCols.get(PSM_File.PSMattr.glyco_q_val)] = "";  // leave MBG q-value blank
    }


    /**
     * Generate a 'dummy' PSM for a newly matched glycan from an existing base PSM. Edit the base PSM to
     * convert to the new glycan
     *
     * @param existingPSM
     * @param newGlycan
     * @param psmCols
     * @param params
     * @param excludeScanNumbers
     */
    public static List<PSM> makeDummyPSM(PSM existingPSM, Glycan newGlycan, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params, int rankCounter, RT_adjuster rtAdjuster, boolean ifdecoy, String parentGlycan, String addedResidual, boolean if_new_glycan, HashSet<Integer> excludeScanNumbers) {
        String[] newSplits = new String[existingPSM.splits.length];
        System.arraycopy(existingPSM.splits, 0, newSplits, 0, newSplits.length);

        // replace values with new glycan
        newSplits[psmCols.get(PSM_File.PSMattr.glycan)] = newGlycan.toPSMString();
        // assigned mods
        String oldAssignedMods = existingPSM.splits[psmCols.get(PSM_File.PSMattr.assigned_mods)];
        ArrayList<String> newAssignedMods = new ArrayList<>();
        for (String mod : oldAssignedMods.split(",")) {
            Matcher siteMatch = modSitePattern.matcher(mod);
            if (siteMatch.find()) {
                int site = Integer.parseInt(siteMatch.group(1));
                if (site != 1) {
                    System.out.println(1);
                }
                if (existingPSM.peptideGlycositeIndexes.contains(site)) {
                    // edit this mod to have the mass of the new glycan
                    if (site == existingPSM.peptideGlycositeIndex) {
                        newAssignedMods.add(String.format("%d%s(%.4f)", site, existingPSM.peptide.charAt(site - 1), newGlycan.mass));
                    }
                } else {
                    newAssignedMods.add(mod);
                }
            } else {
                // include TMT to assigned mods
                newAssignedMods.add(mod);
            }
        }
        // todo: fix the assigned_mod
        newSplits[psmCols.get(PSM_File.PSMattr.assigned_mods)] = String.join(",", newAssignedMods);

        // mod peptide
        String oldGlycEntry = String.format("[%d]", Math.round(existingPSM.glycan.mass + AA_masses.aa_masses.get(existingPSM.peptide.charAt(existingPSM.peptideGlycositeIndex - 1))));
        String newGlycEntry = String.format("[%d]", Math.round(newGlycan.mass + AA_masses.aa_masses.get(existingPSM.peptide.charAt(existingPSM.peptideGlycositeIndex - 1))));
        String newModPep = existingPSM.splits[psmCols.get(PSM_File.PSMattr.mod_peptide)].replace(oldGlycEntry, newGlycEntry);
        newSplits[psmCols.get(PSM_File.PSMattr.mod_peptide)] = newModPep;

        // update calculated masses and delta mass for new glycan
        double diff_mass = newGlycan.mass - existingPSM.glycan.mass;
        double diff_mz = (newGlycan.mass - existingPSM.glycan.mass) / (double) existingPSM.charge;
        newSplits[psmCols.get(PSM_File.PSMattr.Calculated_MASS)] = String.format("%.4f", existingPSM.calcMass + diff_mass);
        newSplits[psmCols.get(PSM_File.PSMattr.Calculated_MZ)] = String.format("%.4f", existingPSM.calcMz + diff_mz);
        // update calibrated masses/mz, observed mass/mz
        newSplits[psmCols.get(PSM_File.PSMattr.CALIBRATED_OBSERVED_MASS)] = String.format("%.4f", Double.parseDouble(existingPSM.splits[psmCols.get(PSM_File.PSMattr.CALIBRATED_OBSERVED_MASS)]) + diff_mass);
        newSplits[psmCols.get(PSM_File.PSMattr.CALIBRATED_OBSERVED_MZ)] = String.format("%.4f", Double.parseDouble(existingPSM.splits[psmCols.get(PSM_File.PSMattr.CALIBRATED_OBSERVED_MZ)]) + diff_mz);
        newSplits[psmCols.get(PSM_File.PSMattr.Observed_MASS)] = String.format("%.4f", Double.parseDouble(existingPSM.splits[psmCols.get(PSM_File.PSMattr.Observed_MASS)]) + diff_mass);
        newSplits[psmCols.get(PSM_File.PSMattr.Observed_Mz)] = String.format("%.4f", Double.parseDouble(existingPSM.splits[psmCols.get(PSM_File.PSMattr.Observed_Mz)]) + diff_mz);

        if (ifdecoy) {
            newSplits[psmCols.get(PSM_File.PSMattr.Calculated_MZ)] = String.format("%.4f", existingPSM.calcMz + ((newGlycan.mass - existingPSM.glycan.mass + 1.0005*11) / (double) existingPSM.charge));
        }
        newSplits[psmCols.get(PSM_File.PSMattr.delta_mass)] = String.format("%.4f", existingPSM.deltaMass);

        // adjust RT and IM
        double newRT = rtAdjuster.getShiftedRT(existingPSM.glycan, newGlycan, existingPSM.apexRT);  // in minute
        double newIM = rtAdjuster.getShiftedIM(existingPSM.glycan, newGlycan, existingPSM.apexIM);
        double newMZ = Double.parseDouble(newSplits[psmCols.get(PSM_File.PSMattr.Calculated_MZ)]);
        newSplits[psmCols.get(PSM_File.PSMattr.RT)] = String.format("%.3f", newRT * 60);
        if (!params.noPASEF) {
            newSplits[psmCols.get(PSM_File.PSMattr.IM)] = String.format("%.3f", newIM);
        }

        List<PSM> addedPSMs = new ArrayList<>();
        float aveMassDiff = 0;
        float imDiff = 0;
        float rtDiff = 0;

        // match against MS1, get intensity, rt duration and average mass
        Peak peakMatcher = new Peak();
        peakMatcher.matchMS1((float) newMZ, (float) newRT, (float) newIM, existingPSM.charge, 0);
        if (peakMatcher.quantifiedEntry == null) {
            newSplits[psmCols.get(PSM_File.PSMattr.intensity)] = String.valueOf(0);
            PSM addedPSM = new PSM(String.join("\t", newSplits), psmCols, params, true, false, 0, parentGlycan, addedResidual, 0);
            // update apex RT IM for iteration
            addedPSM.apexRT = newRT;
            addedPSM.apexIM = newIM;
            addedPSMs.add(addedPSM);
            return addedPSMs;
        }
        // to find MS2
        List<Spectrum> ms2Candidates;
        if (params.noPASEF) {
            ms2Candidates = MSMS.spectraSearcher(peakMatcher.quantifiedEntry.rtStart, peakMatcher.quantifiedEntry.rtEnd, newMZ, excludeScanNumbers);
        } else {
            ms2Candidates = MSMS.spectraSearcher(peakMatcher.quantifiedEntry.rtStart, peakMatcher.quantifiedEntry.rtEnd, newMZ, peakMatcher.quantifiedEntry.imStart, peakMatcher.quantifiedEntry.imEnd, excludeScanNumbers);
        }
        // update intensity
        newSplits[psmCols.get(PSM_File.PSMattr.intensity)] = String.valueOf(peakMatcher.quantifiedEntry.intensity);
        // fill in averageMassDifference, RT and IM difference
        aveMassDiff = rtAdjuster.getAverageMassDiff(existingPSM.glycan, newGlycan, existingPSM.averageMass, peakMatcher.averageMass, existingPSM.charge);
        rtDiff = Math.abs((float) (peakMatcher.quantifiedEntry.apexRt - newRT));
        imDiff = Math.abs((float) (peakMatcher.quantifiedEntry.apexIm - newIM));
        // get average mass error in ppm for scoring
        double c_avePPM = getAvePPM(peakMatcher.quantifiedEntry);

        // todo update spectra information if ms2 candidates found
        int shiftFrequency = rtAdjuster.getShiftFrequency(existingPSM.glycan, newGlycan);
        double featureIntensity = 0; // intensity sum of y0 and y1

        // leaving specified columns blank
        setBlankColumns(newSplits, psmCols);

        if (!ms2Candidates.isEmpty()) {
            // match Y1, Y0
            List<Double> yIntensities;
            yIntensities = MSMS.y1Matcher(ms2Candidates, existingPSM, params);
            featureIntensity = yIntensities.stream()
                    .mapToDouble(Double::doubleValue) // Convert to double
                    .average()
                    .orElse(Double.NaN); // Handle empty list case

            for (int i = 0; i < ms2Candidates.size(); i++) {
                newSplits[psmCols.get(PSM_File.PSMattr.spectrum)] = String.join(".", new String[]{existingPSM.rawfile, String.format("%0" + MSMS.maxScanNumLength + "d", ms2Candidates.get(i).scanNum),
                        String.format("%0" + MSMS.maxScanNumLength + "d", ms2Candidates.get(i).scanNum), String.valueOf(ms2Candidates.get(i).charge)});
                // generate new PSM
                PSM addedPSM = new PSM(String.join("\t", newSplits), psmCols, params, true, true, yIntensities.get(i), parentGlycan, addedResidual, 0);
                addedPSM.updateRank(existingPSM.rank + rankCounter, psmCols);
                addedPSM.rtDiff = rtDiff;
                addedPSM.imDiff = imDiff;
                addedPSM.avePPM = c_avePPM;
                addedPSM.kl = peakMatcher.quantifiedEntry.kl;
                addedPSM.shiftFrequency = shiftFrequency;
                addedPSM.featureIntensity = yIntensities.get(i);  // use real intensity not average all MS2 hits
                // update MBG related columns (some columns already updated in generating new psm)
                addedPSM.updateSplits(psmCols, params.noPASEF, if_new_glycan);
                // update apex RT IM for iteration
                addedPSM.apexRT = newRT;
                addedPSM.apexIM = newIM;
                addedPSMs.add(addedPSM);
            }

        } else {
            // generate new PSM
            PSM addedPSM = new PSM(String.join("\t", newSplits), psmCols, params, true, false, 0, parentGlycan, addedResidual,  0);
            addedPSM.updateRank(existingPSM.rank + rankCounter, psmCols);
            addedPSM.rtDiff = rtDiff;
            addedPSM.imDiff = imDiff;
            addedPSM.avePPM = c_avePPM;
            addedPSM.kl = peakMatcher.quantifiedEntry.kl;
            addedPSM.shiftFrequency = shiftFrequency;
            addedPSM.featureIntensity = featureIntensity;
            // update MBG related columns (some columns already updated in generating new psm)
            addedPSM.updateSplits(psmCols, params.noPASEF, if_new_glycan);
            // update apex RT IM for iteration
            addedPSM.apexRT = newRT;
            addedPSM.apexIM = newIM;
            addedPSMs.add(addedPSM);
        }
        // add intensity, rtDiff, imDiff, glycan shift frequency to decoys or targets for scoring

        if (!params.noPASEF) {
            if (ifdecoy) {
                Scoring.decoyData.add(new double[] {(double) peakMatcher.quantifiedEntry.intensity, c_avePPM, peakMatcher.quantifiedEntry.kl, shiftFrequency, featureIntensity, (double) rtDiff, (double) imDiff});
            } else {
                Scoring.targetData.add(new double[] {(double) peakMatcher.quantifiedEntry.intensity, c_avePPM, peakMatcher.quantifiedEntry.kl, shiftFrequency, featureIntensity, (double) rtDiff, (double) imDiff});
            }
        } else {
            if (ifdecoy) {
                Scoring.decoyData.add(new double[] {(double) peakMatcher.quantifiedEntry.intensity, c_avePPM, peakMatcher.quantifiedEntry.kl, shiftFrequency, featureIntensity, (double) rtDiff});
            } else {
                Scoring.targetData.add(new double[] {(double) peakMatcher.quantifiedEntry.intensity, c_avePPM, peakMatcher.quantifiedEntry.kl, shiftFrequency, featureIntensity, (double) rtDiff});
            }
        }
        return addedPSMs;
    }

    private static double getAvePPM(Entry quantifiedEntry) {
        double avePPM = 0;
        int nMatchedIsotope = 0;
        if (!Double.isNaN(quantifiedEntry.ppm0)) {
            avePPM += Math.abs(quantifiedEntry.ppm0);
            nMatchedIsotope += 1;
        }
        if (!Double.isNaN(quantifiedEntry.ppm1)) {
            avePPM += Math.abs(quantifiedEntry.ppm1);
            nMatchedIsotope += 1;
        }
        if (!Double.isNaN(quantifiedEntry.ppm2)) {
            avePPM += Math.abs(quantifiedEntry.ppm2);
            nMatchedIsotope += 1;
        }
        if (nMatchedIsotope == 0) {
            return 0;
        } else {
            return avePPM / nMatchedIsotope;
        }
    }

    static Map<String, Map<String, Float>> glycanMassShift = new HashMap<>();
//    {{
//        put("HexNAc(2)Hex(5)", 109.9F);
//        put("HexNAc(2)Hex(6)", 68.1F);
//        put("HexNAc(2)Hex(7)", 42.4F);
//        put("HexNAc(2)Hex(8)", 19.9F);
//        put("HexNAc(2)Hex(9)", 11.3F);
//        put("HexNAc(2)Hex(10)", 0F);
//    }};

    private static void setBlankColumns(String[] splits, HashMap<PSM_File.PSMattr, Integer> psmCols){
        for (PSM_File.PSMattr attr : PSM_File.PSMattr_blank)
            if (psmCols.containsKey(attr)) {
                splits[psmCols.get(attr)] = "";
            }
//        splits[psmCols.get(PSM_File.PSMattr.hyperscore)] = String.valueOf(0);
//        splits[psmCols.get(PSM_File.PSMattr.nextscore)] = String.valueOf(0);
//        splits[psmCols.get(PSM_File.PSMattr.probability)] = String.valueOf(0);
        if (psmCols.containsKey(PSM_File.PSMattr.glyco_q_val)) {
            splits[psmCols.get(PSM_File.PSMattr.glyco_q_val)] = String.valueOf(0);
        }
//        splits[psmCols.get(PSM_File.PSMattr.expectation)] = String.valueOf(1);
    }

    private static double getShiftedRT(Glycan oldGlycan, Glycan newGlycan, double oldRT, RT_adjuster rtAdjuster) {
        //
        double newRT;
        if (RT_adjuster.monosaccharideRtShifts.containsKey(newGlycan.name)) {
            newRT = oldRT + glycanMassShift.get(newGlycan.name).getOrDefault(oldGlycan.name, 0f);
        } else
            newRT = oldRT;
//        double newRT = (oldRT + glycanMassShift.getOrDefault(newGlycan.name, 0f)) - glycanMassShift.getOrDefault(oldGlycan.name, 0f);
        return newRT;
    }

    public static void getBasePSM(GlycoSite site, Glycan newGlycan) {
        // get base glycan
        boolean monoMassDiff = false;
        PSM base_psm = null;
        for (PSM c_psm: site.psms) {
            Integer glycanDiff = glycanComparator(c_psm.glycan.composition, newGlycan.composition);
            if (glycanDiff == 1) {
                monoMassDiff = true;
                base_psm = c_psm;
            }
        site.basePSM = base_psm;
        return;

        }

        // determine RT shift
    }

    private static int glycanComparator(Map <GlycanResidue, Integer> glycan1, Map <GlycanResidue, Integer> glycan2) {
        int difference = 0;
        for (GlycanResidue key : glycan1.keySet()) {
            difference += Math.abs(glycan1.getOrDefault(key, 0) - glycan2.getOrDefault(key, 0));
        }
        return difference;
    }

    public void add_MBG_score(HashMap<PSM_File.PSMattr, Integer> psmCols, double score) {
        splits[psmCols.get(PSM_File.PSMattr.MBG_score)] = String.format("%.3f", score);
    }

    /**
     * Given a glycan, adding NH3 to generate a new glycan
     *
     * @return
     */
    public static Glycan addResidual(Glycan baseGlycan, GlycanResidue mod) {
        TreeMap<GlycanResidue, Integer> newComposition = new TreeMap(baseGlycan.composition);
        newComposition.put(mod, baseGlycan.composition.getOrDefault(mod, 0) + 1);
        return new Glycan(newComposition);
    }

    /**
    * Generate all glycans that have a given residual difference from the current psm
    * */
    public static ArrayList<Glycan> getDummyGlycans(PSM psm, HashMap<String, GlycanResidue> residues) {
        ArrayList<Glycan> dummy_glycans = new ArrayList<>();
        for (Map.Entry<String, GlycanResidue> residualEntry: residues.entrySet()) {
            GlycanResidue residual = residualEntry.getValue();
            if (psm.glycan.composition.containsKey(residual)) {
                int count = psm.glycan.composition.get(residual);
                TreeMap<GlycanResidue, Integer> newCompositionAdd = new TreeMap<>(psm.glycan.composition);
                newCompositionAdd.put(residual, count + 1);
                dummy_glycans.add(new Glycan(newCompositionAdd));
                // loss
                TreeMap<GlycanResidue, Integer> newCompositionLoss = new TreeMap<>(psm.glycan.composition);
                if (count > 1) {
                    newCompositionLoss.put(residual, count - 1);
                    dummy_glycans.add(new Glycan(newCompositionLoss));
                }
                else {
                    newCompositionLoss.remove(residual);
                    Glycan new_glycan = new Glycan(newCompositionLoss);
                    // in case of glycan with none glycan residue
                    if (new_glycan.toString().length() < 2) {
                        continue;
                    }
                    dummy_glycans.add(new_glycan);
                }
            }
            else {
                TreeMap<GlycanResidue, Integer> newCompositionAdd = new TreeMap<>(psm.glycan.composition);
                newCompositionAdd.put(residual, 1);
                dummy_glycans.add(new Glycan(newCompositionAdd));
            }
        }
        return dummy_glycans;
    }

    public static HashMap<Glycan, String> getDummyGlycansV2(PSM psm, List<Glycan> toAddglycans) {
        ArrayList<Glycan> dummy_glycans = new ArrayList<>();
        ArrayList<String> addedResiduals = new ArrayList<>();
        for (Glycan cGlycan: toAddglycans) {
            // generate dummy PSMs with gained or lost glycan
            // gained
            saccharidesAdder(psm.glycan, cGlycan, dummy_glycans, true);
            addedResiduals.add(cGlycan.toString());
            // loss, only add to addedResiduals when eligible
            if (saccharidesAdder(psm.glycan, cGlycan, dummy_glycans, false)) {
                // make the composition counts negative to get the loss residual
                TreeMap<GlycanResidue, Integer> lossComposition = new TreeMap<>(cGlycan.composition);
                for (Map.Entry<GlycanResidue, Integer> entry: lossComposition.entrySet()) {
                    entry.setValue(entry.getValue() * -1);
                }
                Glycan lossResidualGlycan = new Glycan(lossComposition);
                addedResiduals.add(lossResidualGlycan.toString());
            }
        }
        HashMap<Glycan, String> glycanMap = new HashMap<>();
        for (int i = 0; i < dummy_glycans.size(); i++) {
            glycanMap.put(dummy_glycans.get(i), addedResiduals.get(i));
        }
        return glycanMap;
    }

    private static boolean saccharidesAdder(Glycan oldGlycan, Glycan glycanToAdd, ArrayList<Glycan> dummy_glycans, boolean ifAdd) {
        TreeMap<GlycanResidue, Integer> newCompositionAdd = new TreeMap<>(oldGlycan.composition);
        if (ifAdd) {
            for (Map.Entry<GlycanResidue, Integer> residueIntegerEntry: glycanToAdd.composition.entrySet()) {
                GlycanResidue cResidual = residueIntegerEntry.getKey();
                if (oldGlycan.composition.containsKey(cResidual)) {
                    newCompositionAdd.put(cResidual, residueIntegerEntry.getValue() + oldGlycan.composition.get(cResidual));
                } else {
                    newCompositionAdd.put(cResidual, residueIntegerEntry.getValue());
                }
            }
        } else {
            for (Map.Entry<GlycanResidue, Integer> residueIntegerEntry : glycanToAdd.composition.entrySet()) {
                GlycanResidue cResidual = residueIntegerEntry.getKey();
                int lossN = residueIntegerEntry.getValue();
                int oldN = oldGlycan.composition.getOrDefault(cResidual, 0);
                if (oldN < lossN) {
                    return false;
                } else if (oldN == lossN) {
                    newCompositionAdd.remove(cResidual);
                } else {
                    newCompositionAdd.put(cResidual, oldN - lossN);
                }
            }
        }
        // if empty after subtraction, ignore this empty glycan and return False
        if (newCompositionAdd.isEmpty()) {
            return false;
        }
        dummy_glycans.add(new Glycan(newCompositionAdd));
        return true;
    }


    /**
     * Match the glycan composition to one of the assigned mods and use that to get the site within the peptide.
     * Assumes a single glycosite because MBG only makes sense for singly-glycosylated peptides
     * @param assignedMods Assigned Modifications string
     */
    private int getGlycoSiteIndex(String assignedMods) {
        String[] mods = assignedMods.split(",");
        for (String mod: mods) {
            Matcher massMatch = modMassPattern.matcher(mod);
            if (massMatch.find()) {
                double mass = Double.parseDouble(massMatch.group(1));
                if (Math.abs(mass - glycan.mass) < 0.001) {
                    // found the assigned mod entry for the glycan. Use its position to compute the glycosite
                    Matcher siteMatch = modSitePattern.matcher(mod);
                    if (siteMatch.find()) {
                        return Integer.parseInt(siteMatch.group(1));
                    }
                }
            }
        }
        // todo return correct glycosite Index when there is multiple glycan on peptide
        return 1;
    }

    private ArrayList<Integer> getGlycoSiteIndices(String assignedMods) {
        // Step 1: Parse the string to get a list of (mass, position) pairs.
        List<MassSitePair> allMods = new ArrayList<>();
        String[] mods = assignedMods.split(",");

        for (String mod : mods) {
            Matcher massMatch = modMassPattern.matcher(mod);
            Matcher siteMatch = modSitePattern.matcher(mod);

            if (massMatch.find() && siteMatch.find()) {
                double mass = Double.parseDouble(massMatch.group(1));
                int position = Integer.parseInt(siteMatch.group(1));
                allMods.add(new MassSitePair(mass, position));
            }
        }

        // Step 2: Use a recursive helper function to find the first valid combination.
        return findFirstValidCombination(
                glycan.mass,
                allMods,
                0,
                0.0,
                new ArrayList<>()
        );
    }

    // A helper class to group mass and its corresponding position.
    private class MassSitePair {
        double mass;
        int position;

        MassSitePair(double mass, int position) {
            this.mass = mass;
            this.position = position;
        }
    }

    // Recursive backtracking method to find the first combination that sums to the target mass.
    private ArrayList<Integer> findFirstValidCombination(
            double targetMass,
            List<MassSitePair> allMods,
            int startIndex,
            double currentSum,
            List<MassSitePair> currentCombination) {

        // Base case: Success
        // If the current sum is very close to the target, we've found a valid combination.
        if (Math.abs(currentSum - targetMass) < 0.001) {
            // Collect and return the positions from this first valid combination.
            ArrayList<Integer> result = new ArrayList<>();
            for (MassSitePair pair : currentCombination) {
                result.add(pair.position);
            }
            return result;
        }

        // Base case: Failure
        // If we've exceeded the target or run out of masses to check, stop this path.
        if (currentSum > targetMass || startIndex >= allMods.size()) {
            return null; // Return null to signal that this path didn't lead to a solution.
        }

        // Recursive step: explore all possibilities.
        for (int i = startIndex; i < allMods.size(); i++) {
            MassSitePair nextMod = allMods.get(i);

            // Include the current modification in the combination.
            currentCombination.add(nextMod);

            // Recurse with the new sum and the next index.
            ArrayList<Integer> result = findFirstValidCombination(
                    targetMass,
                    allMods,
                    i + 1,
                    currentSum + nextMod.mass,
                    currentCombination
            );

            // If a result was found in a deeper call, pass it up immediately.
            if (result != null) {
                return result;
            }

            // Backtrack: remove the modification to explore other combinations.
            currentCombination.remove(currentCombination.size() - 1);
        }

        // If no combination was found in any path, return araylist[1].
        ArrayList<Integer> result = new ArrayList<>();
        result.add(1);
        return result;
    }

    public String getGlycoSiteID() {
        return String.format("%s_%d", protein, proteinGlycositeIndex);
    }

    public int getScanNum(HashMap<PSM_File.PSMattr, Integer> psmColumns) {
        String specStr = splits[psmColumns.get(PSM_File.PSMattr.spectrum)];
        String[] specSplits = specStr.split("\\.");
        return Integer.parseInt(specSplits[1]);
    }

    public String getRawfile(HashMap<PSM_File.PSMattr, Integer> psmColumns) {
        String specStr = splits[psmColumns.get(PSM_File.PSMattr.spectrum)];
        String[] specSplits = specStr.split("\\.");
        return specSplits[0];
    }

    @Override
    public int compareTo(PSM o) {
        if (rawfile.equals(o.rawfile)) {
            return Integer.compare(scanNum, o.scanNum);
        } else {
            return rawfile.compareTo(o.rawfile);
        }
    }
}
