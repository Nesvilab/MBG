package com.mbg;

import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanResidue;

import java.util.*;

import static java.lang.System.out;

/**
 * RT_adjuster handles the calculation and management of retention time (RT) and 
 * ion mobility (IM) shifts between different glycan compositions.
 * 
 * This class computes median shifts in RT and IM values for glycan transitions, 
 * which are then used to predict RT and IM values for inferred glycopeptides 
 * during the glycoform inference process.
 */
public class RT_adjuster {
    /** Maps from source glycan to target glycan to list of observed RT shifts */
    static HashMap<String, HashMap<String, List<Float>>> allShiftedRTs;
    
    /** Maps from source glycan to target glycan to median RT shift */
    static HashMap<String, HashMap<String, Float>> shiftedRTs;
    
    /** Maps from source glycan to target glycan to list of observed IM shifts */
    static HashMap<String, HashMap<String, List<Float>>> allShiftedIMs;
    
    /** Maps from source glycan to target glycan to median IM shift */
    static HashMap<String, HashMap<String, Float>> shiftedIMs;
    
    /** Maps from monosaccharide difference to list of observed RT shifts */
    static HashMap<String, List<Float>> allShiftedMonoRTs;
    
    /** Maps from monosaccharide difference to list of observed IM shifts */
    static HashMap<String, List<Float>> allShiftedMonoIMs;
    
    /** Maps from monosaccharide difference to median RT shift */
    static HashMap<String, Float> monosaccharideRtShifts;
    
    /** Maps from monosaccharide difference to median IM shift */
    static HashMap<String, Float> monosaccharideImShifts;

    /**
     * Constructs an RT_adjuster based on observed glycosites.
     *
     * @param observedGlycoSites Map of peptide sequences to observed glycosites
     * @param combinedShiftedRTs
     * @param combinedShiftedIMs
     */
    public RT_adjuster(HashMap<String, GlycoSite> observedGlycoSites, Integer minPSMs, HashMap<String, List<Float>> combinedShiftedRTs, HashMap<String, List<Float>> combinedShiftedIMs) {
        // Initialize data structures
        allShiftedRTs = new HashMap<>();
        shiftedRTs = new HashMap<>();
        allShiftedIMs = new HashMap<>();
        shiftedIMs = new HashMap<>();
        allShiftedMonoRTs = new HashMap<>();
        allShiftedMonoIMs = new HashMap<>();
        monosaccharideRtShifts = new HashMap<>();
        monosaccharideImShifts = new HashMap<>();
        
        // Process observed glycosites to extract RT and IM shifts
        matchMS1forPSMs(observedGlycoSites, minPSMs);
        
        // Calculate median shifts for use in predictions
        calculateMedians();
        // combine_shifted RT/IMs
        combineAllShiftedRTIMs(allShiftedRTs, combinedShiftedRTs);
        combineAllShiftedRTIMs(allShiftedIMs, combinedShiftedIMs);
        combineAllShiftedMonoRTIMs(allShiftedMonoRTs, combinedShiftedRTs);
        combineAllShiftedMonoRTIMs(allShiftedMonoIMs, combinedShiftedIMs);
    }

    /**
     * Matches MS1 features for PSMs and collects RT and IM shift data.
     * 
     * @param observedGlycoSites Map of peptide sequences to observed glycosites
     */
    public void matchMS1forPSMs(HashMap<String, GlycoSite> observedGlycoSites, Integer minPSMs) {
        Peak peakMatcher = new Peak();
        
        // Process each glycosite
        for (Map.Entry<String, GlycoSite> siteEntry : observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            HashMap<Glycan, Float> glycanRTs = new HashMap<>();
            HashMap<Glycan, Float> glycanIMs = new HashMap<>();
            
            // For each glycan observed at this site
            for (Map.Entry<String, ArrayList<PSM>> glycanEntry: site.observedGlycans.entrySet()) {
                ArrayList<PSM> glycanPsms = glycanEntry.getValue();
                if (glycanPsms.size() < minPSMs) {
                    continue;
                }
                Glycan currentGlycan = glycanPsms.get(0).glycan;
                // Find the PSM with the highest intensity for this glycan
                float topApexRT = 0;
                float topIntensity = 0;
                PSM topPSM = null;
                
                for (PSM psm : glycanPsms) {
                    // Match MS1 feature
                    peakMatcher.matchMS1(
                        (float) psm.calcMz, 
                        (float) psm.rt / 60, 
                        (float) psm.im, 
                        psm.charge, 
                        0
                    );
                    
                    // Update PSM with matched feature information
                    if (peakMatcher.quantifiedEntry != null) {
                        psm.apexRT = peakMatcher.quantifiedEntry.apexRt;
                        psm.apexIM = peakMatcher.quantifiedEntry.apexIm;
                        psm.averageMass = peakMatcher.averageMass;
                        psm.intensity = peakMatcher.quantifiedEntry.intensity;
                        
                        // Track the highest intensity PSM
                        if (psm.intensity > topIntensity) {
                            topApexRT = (float) psm.apexRT;
                            topIntensity = (float) psm.intensity;
                            topPSM = psm;
                        }
                    }
                }
                
                // Store RT and IM for the highest intensity PSM
                if (topApexRT > 0) {
                    glycanRTs.put(currentGlycan, topApexRT);
                    glycanIMs.put(currentGlycan, (float) topPSM.apexIM);
                    site.basePSMs.add(topPSM);
                }
            }
            
            // Collect RT and IM shifts between different glycans
            collectShiftedRtIm(glycanRTs, "RT");
            collectShiftedRtIm(glycanIMs, "IM");
        }
    }

    /**
     * Collects observed shifts in RT or IM between glycan pairs.
     * 
     * @param glycanValues Map of glycans to their RT or IM values
     * @param mode "RT" or "IM" to specify which type of shift to collect
     */
    private static void collectShiftedRtIm(HashMap<Glycan, Float> glycanValues, String mode) {
        // Compare each glycan with every other glycan
        for (Map.Entry<Glycan, Float> entry1: glycanValues.entrySet()) {
            Glycan glycan1 = entry1.getKey();
            Float value1 = entry1.getValue();
            
            for (Map.Entry<Glycan, Float> entry2: glycanValues.entrySet()) {
                Glycan glycan2 = entry2.getKey();
                
                // Ensure glycan2 contains all monosaccharides in glycan1
                if (!glycan2.composition.keySet().containsAll(glycan1.composition.keySet())) {
                    continue;
                }
                
                Float value2 = entry2.getValue();
                
                // Check if glycan2 has more of each monosaccharide than glycan1
                boolean validComparison = true;
                int totalDifference = 0;
                
                for (Map.Entry<GlycanResidue, Integer> compEntry : glycan2.composition.entrySet()) {
                    GlycanResidue residue = compEntry.getKey();
                    Integer count2 = compEntry.getValue();
                    Integer count1 = glycan1.composition.getOrDefault(residue, 0);
                    
                    // Glycan2 must have more of each monosaccharide than glycan1
                    if (count2 < count1) {
                        validComparison = false;
                        break;
                    }
                    
                    totalDifference += count2 - count1;
                }
                
                // Store the shift if valid
                if (validComparison && totalDifference > 0) {
                    float shift = value2 - value1;
                    putData(glycan1.toString(), glycan2.toString(), shift, mode);
                    putMonoData(glycan1, glycan2, shift, mode);
                }
            }
        }
    }

    /**
     * Records observed shifts for specific monosaccharide changes.
     * 
     * @param g1 Source glycan
     * @param g2 Target glycan
     * @param value Observed shift value
     * @param mode "RT" or "IM"
     */
    private static void putMonoData(Glycan g1, Glycan g2, Float value, String mode) {
        // Get the monosaccharide difference between glycans
        String shiftedSaccharides = getShiftedSaccharides(g1, g2);
        
        // Store in appropriate collection
        if (Objects.equals(mode, "RT")) {
            allShiftedMonoRTs.computeIfAbsent(shiftedSaccharides, k -> new ArrayList<>());
            allShiftedMonoRTs.get(shiftedSaccharides).add(value);
        } else {
            allShiftedMonoIMs.computeIfAbsent(shiftedSaccharides, k -> new ArrayList<>());
            allShiftedMonoIMs.get(shiftedSaccharides).add(value);
        }
    }

    /**
     * Identifies the single monosaccharide difference between two glycans.
     * 
     * @param g1 Source glycan (smaller)
     * @param g2 Target glycan (larger)
     * @return String representing the monosaccharide that differs by exactly 1
     */
    private static String getShiftedMonosaccharide(Glycan g1, Glycan g2) {
        // g2 should always be larger than g1
        String shiftedMonosaccharide = "";
        
        for (Map.Entry<GlycanResidue, Integer> compEntry2: g2.composition.entrySet()) {
            GlycanResidue glycanResidue = compEntry2.getKey();
            int residueNum2 = compEntry2.getValue();
            
            if (!g1.composition.containsKey(glycanResidue) || 
                residueNum2 - g1.composition.get(glycanResidue) == 1) {
                shiftedMonosaccharide = glycanResidue.toString();
                break;
            }
        }
        
        if (shiftedMonosaccharide.isEmpty()) {
            out.printf("Error: shifted monosaccharide difference not equal to 1 between %s and %s\n", 
                      g1.toString(), g2);
            System.exit(1);
        }
        
        return shiftedMonosaccharide;
    }

    /**
     * Calculates the complete glycan composition difference between two glycans.
     * 
     * @param g1 Source glycan (smaller)
     * @param g2 Target glycan (larger)
     * @return String representing the difference in glycan composition
     */
    private static String getShiftedSaccharides(Glycan g1, Glycan g2) {
        // Create a map to store the composition difference
        TreeMap<GlycanResidue, Integer> compositionDiff = new TreeMap<>();
        
        // Calculate the difference for each monosaccharide
        for (Map.Entry<GlycanResidue, Integer> compEntry2: g2.composition.entrySet()) {
            GlycanResidue glycanResidue = compEntry2.getKey();
            int residueNum2 = compEntry2.getValue();
            int diff = residueNum2 - g1.composition.getOrDefault(glycanResidue, 0);
            
            if (diff > 0) {
                compositionDiff.put(glycanResidue, diff);
            }
        }
        
        // Convert the difference to a glycan string
        String shiftedSaccharides = "";
        if (!compositionDiff.isEmpty()) {
            Glycan glycanDiff = new Glycan(compositionDiff);
            shiftedSaccharides = glycanDiff.toString();
        }
        
        return shiftedSaccharides;
    }

    /**
     * Records an observed shift between two glycans.
     * 
     * @param sourceGlycan Source glycan
     * @param targetGlycan Target glycan
     * @param value Observed shift value
     * @param mode "RT" or "IM"
     */
    private static void putData(String sourceGlycan, String targetGlycan, Float value, String mode) {
        HashMap<String, List<Float>> innerMap;
        
        if (Objects.equals(mode, "RT")) {
            // Ensure the outer map exists
            allShiftedRTs.computeIfAbsent(sourceGlycan, k -> new HashMap<>());
            innerMap = allShiftedRTs.get(sourceGlycan);
        } else {
            allShiftedIMs.computeIfAbsent(sourceGlycan, k -> new HashMap<>());
            innerMap = allShiftedIMs.get(sourceGlycan);
        }

        // Ensure the inner list exists
        innerMap.computeIfAbsent(targetGlycan, k -> new ArrayList<>());
        
        // Add the value to the list
        innerMap.get(targetGlycan).add(value);
    }

    /**
     * Calculates median shifts for all collected RT and IM observations.
     */
    public static void calculateMedians() {
        // Calculate RT median shifts for glycan pairs
        for (String glycanName : allShiftedRTs.keySet()) {
            HashMap<String, List<Float>> innerMap = allShiftedRTs.get(glycanName);
            HashMap<String, Float> medianShifts = new HashMap<>();

            for (String targetGlycan : innerMap.keySet()) {
                List<Float> shifts = innerMap.get(targetGlycan);
                Float median = calculateMedian(shifts);
                medianShifts.put(targetGlycan, median);
            }
            
            shiftedRTs.put(glycanName, medianShifts);
        }
        
        // Calculate IM median shifts for glycan pairs
        for (String glycanName : allShiftedIMs.keySet()) {
            HashMap<String, List<Float>> innerMap = allShiftedIMs.get(glycanName);
            HashMap<String, Float> medianShifts = new HashMap<>();

            for (String targetGlycan : innerMap.keySet()) {
                List<Float> shifts = innerMap.get(targetGlycan);
                Float median = calculateMedian(shifts);
                medianShifts.put(targetGlycan, median);
            }
            
            shiftedIMs.put(glycanName, medianShifts);
        }
        
        // Calculate RT median shifts for monosaccharide changes
        for (String residue : allShiftedMonoRTs.keySet()) {
            monosaccharideRtShifts.put(residue, calculateMedian(allShiftedMonoRTs.get(residue)));
        }
        
        // Calculate IM median shifts for monosaccharide changes
        for (String residue : allShiftedMonoIMs.keySet()) {
            monosaccharideImShifts.put(residue, calculateMedian(allShiftedMonoIMs.get(residue)));
        }
    }

    /**
     * Calculates the median value of a list of floats.
     * 
     * @param values List of float values
     * @return Median value or null if the list is empty
     */
    public static Float calculateMedian(List<Float> values) {
        if (values.isEmpty()) {
            return null;
        }
        
        Collections.sort(values);
        int size = values.size();
        
        if (size % 2 == 0) {
            // Even number of elements - return average of middle two
            return (values.get(size / 2 - 1) + values.get(size / 2)) / 2;
        } else {
            // Odd number of elements - return middle element
            return values.get(size / 2);
        }
    }

    /**
     * Predicts the retention time for a new glycan based on an existing glycan.
     * 
     * @param oldGlycan Source glycan
     * @param newGlycan Target glycan
     * @param oldRT Retention time of the source glycan
     * @return Predicted retention time for the target glycan
     */
    Float getShiftedRT(Glycan oldGlycan, Glycan newGlycan, double oldRT) {
        Float rtShift;
        
        if (oldGlycan.mass < newGlycan.mass) {
            // Try to find a direct glycan-to-glycan shift
            if (shiftedRTs.containsKey(oldGlycan.toString()) && 
                shiftedRTs.get(oldGlycan.toString()).containsKey(newGlycan.toString())) {
                rtShift = shiftedRTs.get(oldGlycan.toString()).get(newGlycan.toString());
            } else {
                // Fall back to monosaccharide-based shift
                String shiftedSaccharides = getShiftedSaccharides(oldGlycan, newGlycan);
                rtShift = monosaccharideRtShifts.getOrDefault(shiftedSaccharides, 0.0f);
            }
        } else {
            // Reverse direction, so negate the shift
            if (shiftedRTs.containsKey(newGlycan.toString()) && 
                shiftedRTs.get(newGlycan.toString()).containsKey(oldGlycan.toString())) {
                rtShift = -shiftedRTs.get(newGlycan.toString()).get(oldGlycan.toString());
            } else {
                // Fall back to monosaccharide-based shift
                String shiftedSaccharides = getShiftedSaccharides(newGlycan, oldGlycan);
                rtShift = -monosaccharideRtShifts.getOrDefault(shiftedSaccharides, 0.0f);
            }
        }
        
        return (float) oldRT + rtShift;
    }

    Float getShiftedIM(Glycan oldGlycan, Glycan newGlycan, double oldIM) {
        // Mirrors getShiftedRT's bidirectional lookup. shiftedIMs is populated only with
        // smaller->larger composition keys, so the loss branch must use newGlycan (the
        // smaller side) as outer key and negate the stored shift. The fallback path uses
        // getShiftedSaccharides (handles multi-residue diffs) rather than
        // getShiftedMonosaccharide, which only supports single-residue diffs and exits
        // the process on anything else.
        Float shiftedIM;
        if (oldGlycan.mass < newGlycan.mass) {
            if (shiftedIMs.containsKey(oldGlycan.toString()) && shiftedIMs.get(oldGlycan.toString()).containsKey(newGlycan.toString())) {
                shiftedIM = shiftedIMs.get(oldGlycan.toString()).get(newGlycan.toString());
            } else {
                String shiftedSaccharides = getShiftedSaccharides(oldGlycan, newGlycan);
                shiftedIM = monosaccharideImShifts.getOrDefault(shiftedSaccharides, (float) 0);
            }
        } else {
            if (shiftedIMs.containsKey(newGlycan.toString()) && shiftedIMs.get(newGlycan.toString()).containsKey(oldGlycan.toString())) {
                shiftedIM = -shiftedIMs.get(newGlycan.toString()).get(oldGlycan.toString());
            } else {
                String shiftedSaccharides = getShiftedSaccharides(newGlycan, oldGlycan);
                shiftedIM = -monosaccharideImShifts.getOrDefault(shiftedSaccharides, (float) 0);
            }
        }
        float newIM = (float) (oldIM + shiftedIM);
        if (newIM < 0) {
            return (float) 0;
        } else {
            return newIM;
        }
    }

    Float getAverageMassDiff(Glycan oldGlycan, Glycan newGlycan, float oldAveMass, float newAveMass, int charge) {
        return (float) (((newGlycan.mass - oldGlycan.mass) - (newAveMass - oldAveMass) * charge)) ;
    }

    Integer getShiftFrequency(Glycan oldGlycan, Glycan newGlycan) {
        // The evidence table allShiftedRTs is populated only with smaller->larger composition
        // pairs, so a direct lookup misses every loss-type candidate (parent->smaller). The
        // count of supporting peptides is direction-agnostic, so also fall back to the reverse
        // direction. Mirrors getShiftedRT, which already handles both directions.
        String oldKey = oldGlycan.toString();
        String newKey = newGlycan.toString();
        if (allShiftedRTs.containsKey(oldKey) && allShiftedRTs.get(oldKey).containsKey(newKey)) {
            return allShiftedRTs.get(oldKey).get(newKey).size();
        }
        if (allShiftedRTs.containsKey(newKey) && allShiftedRTs.get(newKey).containsKey(oldKey)) {
            return allShiftedRTs.get(newKey).get(oldKey).size();
        }
        return 0;
    }

    public static void combineAllShiftedRTIMs(
            HashMap<String, HashMap<String, List<Float>>> c_allShiftedRTs,
            HashMap<String, List<Float>> combinedShifts) {


        for (Map.Entry<String, HashMap<String, List<Float>>> outerEntry : c_allShiftedRTs.entrySet()) {
            String outerKey = outerEntry.getKey();
            HashMap<String, List<Float>> innerMap = outerEntry.getValue();

//            combinedShifts.putIfAbsent(outerKey, new HashMap<>());
//            HashMap<String, List<Float>> combinedInner = combinedShifts.get(outerKey);

            for (Map.Entry<String, List<Float>> innerEntry : innerMap.entrySet()) {
                String innerKey = innerEntry.getKey();
                String c_shiftKey = outerKey + "->" + innerKey;
                List<Float> values = innerEntry.getValue();

                combinedShifts.putIfAbsent(c_shiftKey, new ArrayList<>());
                combinedShifts.get(c_shiftKey).addAll(values);
            }
        }

    }

    public static void combineAllShiftedMonoRTIMs(
            HashMap<String, List<Float>> c_shifts,
            HashMap<String, List<Float>> combinedShifts) {

        for (Map.Entry<String, List<Float>> innerEntry : c_shifts.entrySet()) {
            String innerKey = innerEntry.getKey();
            List<Float> values = innerEntry.getValue();

            combinedShifts.putIfAbsent(innerKey, new ArrayList<>());
            combinedShifts.get(innerKey).addAll(values);
        }
    }
}
