package com.mbg;

import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanResidue;

import java.util.*;

import static java.lang.System.out;

public class RT_adjuster {
    static HashMap<String, HashMap<String, List<Float>>> allShiftedRTs;
    static HashMap<String, HashMap<String, Float>> shiftedRTs;
    static HashMap<String, HashMap<String, List<Float>>> allShiftedIMs;
    static HashMap<String, HashMap<String, Float>> shiftedIMs;
    static HashMap<String, List<Float>> allShiftedMonoRTs;
    static HashMap<String, List<Float>> allShiftedMonoIMs;
    static HashMap<String, Float> monosaccharideRtShifts;
    static HashMap<String, Float> monosaccharideImShifts;

    public RT_adjuster(HashMap<String, GlycoSite> observedGlycoSites) {
        allShiftedRTs = new HashMap<>(); // storage RT shifts from one glycan comp to another
        shiftedRTs = new HashMap<>();
        allShiftedIMs = new HashMap<>(); // storage IM shifts from one glycan comp to another
        shiftedIMs = new HashMap<>();
        allShiftedMonoRTs = new HashMap<>(); // storage (a set of) monosaccharide difference RTs
        allShiftedMonoIMs = new HashMap<>(); // storage (a set of) monosaccharide difference IMs
        monosaccharideRtShifts = new HashMap<>();
        monosaccharideImShifts = new HashMap<>();
        matchMS1forPSMs(observedGlycoSites);
        calculateMedians();
    }

    public void matchMS1forPSMs(HashMap<String, GlycoSite> observedGlycoSites) {
        Peak peakMatcher = new Peak();
        for (Map.Entry<String, GlycoSite> siteEntry : observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            HashMap<Glycan, Float> glycanRTs = new HashMap<>();
            HashMap<Glycan, Float> glycanIMs = new HashMap<>();
            for (Map.Entry<String, ArrayList<PSM>> glycanEntry: site.observedGlycans.entrySet()) {
                ArrayList<PSM> glycan_psm = glycanEntry.getValue();
                Glycan cGlycan = glycan_psm.get(0).glycan;

                float topApexRT = 0;  // get top intensity apex RT (in minute) of this glycopeptide
                float topIntensity = 0;
                PSM topPSM = null;
                for (PSM psm : glycan_psm) {
                    peakMatcher.matchMS1((float) psm.calcMz, (float) psm.rt / 60, (float) psm.im, psm.charge, 0);
                    if (peakMatcher.quantifiedEntry != null) {
                        psm.apexRT = peakMatcher.quantifiedEntry.apexRt;
                        psm.apexIM = peakMatcher.quantifiedEntry.apexIm;
                        psm.averageMass = peakMatcher.averageMass;
                        psm.intensity = peakMatcher.quantifiedEntry.intensity;
                        if (psm.intensity > topIntensity) {
                            topApexRT = (float) psm.apexRT;
                            topIntensity = (float) psm.intensity;
                            topPSM = psm;
                        }
                    }
                }
                if (topApexRT > 0) {
                    glycanRTs.put(cGlycan, topApexRT);
                    glycanIMs.put(cGlycan, (float) topPSM.apexIM);
                    site.basePSMs.add(topPSM);
                }
            }
            // storage into RT adjuster for calculating RT median shift
            collectShiftedRtIm(glycanRTs, "RT");
            collectShiftedRtIm(glycanIMs, "IM");
        }
    }

    private static void collectShiftedRtIm(HashMap<Glycan, Float> glycanRTs, String mode) {
        // iterate
        for (Map.Entry<Glycan, Float> rtShiftEntry1: glycanRTs.entrySet()) {
            Glycan glycan1 = rtShiftEntry1.getKey();
            Float shifedRT1 = rtShiftEntry1.getValue();
            for (Map.Entry<Glycan, Float> rtShiftEntry2: glycanRTs.entrySet()) {
                Glycan glycan2 = rtShiftEntry2.getKey();
                // if glycan1 contain monosaccharide that not in glycan2
                if (!glycan2.composition.keySet().containsAll(glycan1.composition.keySet())) {
                    continue;
                }
                Float shifedRT2 = rtShiftEntry2.getValue();
                // monosaccharide difference
                Iterator<Map.Entry<GlycanResidue, Integer>> var2 = glycan2.composition.entrySet().iterator();
                boolean ifAdd = true;
                Integer totalDiff = 0;
                while (var2.hasNext()) {
                    Map.Entry<GlycanResidue, Integer> compEntry = var2.next();
                    GlycanResidue residueName = compEntry.getKey();
                    Integer count = glycan1.composition.getOrDefault(residueName, 0);
                    Integer count2 = compEntry.getValue();
                    // glycan2 must greater than glycan1 and contains all monosaccharides of glycan1
                    if (count2 < count) {
                        ifAdd = false;
                    }
                    totalDiff += count2 - count;
//                    if (totalDiff > 1) {
//                        ifAdd = false;
//                    }
                }
                if (ifAdd && totalDiff > 0) {
                    putData(glycan1.toString(), glycan2.toString(), shifedRT2-shifedRT1, mode);
                    putMonoData(glycan1, glycan2, shifedRT2-shifedRT1, mode);
                }
            }
        }
    }

    private static void putMonoData(Glycan g1, Glycan g2, Float value, String mode) {
        // get the shifted monosaccharide difference, g2 should always bigger than g1
        String shiftedSaccharides = getShiftedSaccharides(g1, g2);
        if (Objects.equals(mode, "RT")) {
            allShiftedMonoRTs.computeIfAbsent(shiftedSaccharides, k -> new ArrayList<>());
            allShiftedMonoRTs.get(shiftedSaccharides).add(value);
        } else {
            allShiftedMonoIMs.computeIfAbsent(shiftedSaccharides, k -> new ArrayList<>());
            allShiftedMonoIMs.get(shiftedSaccharides).add(value);
        }
    }

    private static String getShiftedMonosaccharide(Glycan g1, Glycan g2) {
        // g2 should always bigger than g1
        String shiftedMonosaccharide = "";
        for (Map.Entry<GlycanResidue, Integer> compEntry2: g2.composition.entrySet()) {
            GlycanResidue glycanResidue = compEntry2.getKey();
            int residueNum2 = compEntry2.getValue();
            if (!g1.composition.containsKey(glycanResidue) || residueNum2 - g1.composition.get(glycanResidue) == 1) {
                shiftedMonosaccharide = glycanResidue.toString();
                break;
            }
        }
        if (shiftedMonosaccharide.isEmpty()) {
            out.printf("Error: shifted monosaccharide difference not equal to 1 between %s and %s\n", g1.toString(), g2);
            System.exit(1);
        }
        return shiftedMonosaccharide;
    }

    /** get saccharide difference between two glycans
     *
     * @param g1 smaller glycan
     * @param g2 bigger glycan
     * @return
     */
    private static String getShiftedSaccharides(Glycan g1, Glycan g2) {
        String shiftedSaccharides = "";
        HashMap<GlycanResidue, Integer> compositionDiff = new HashMap<>();
        for (Map.Entry<GlycanResidue, Integer> compEntry2: g2.composition.entrySet()) {
            GlycanResidue glycanResidue = compEntry2.getKey();
            int residueNum2 = compEntry2.getValue();
            int diff = residueNum2 - g1.composition.getOrDefault(glycanResidue, 0);
            if (diff > 0) {
                compositionDiff.put(glycanResidue, diff);
            }
        }
        if (!compositionDiff.isEmpty()) {
            Glycan glycanDiff = new Glycan(compositionDiff);
            shiftedSaccharides = glycanDiff.toString();
        }
        return shiftedSaccharides;
    }

    private static void putData(String outerKey, String innerKey, Float value, String mode) {
        HashMap<String, List<Float>> innerMap;
        if (Objects.equals(mode, "RT")) {
            // Check if the outer key exists, if not create a new HashMap
            allShiftedRTs.computeIfAbsent(outerKey, k -> new HashMap<>());

            // Get the inner HashMap for the outer key
            innerMap = allShiftedRTs.get(outerKey);
        } else {
            allShiftedIMs.computeIfAbsent(outerKey, k -> new HashMap<>());
            innerMap = allShiftedIMs.get(outerKey);
        }

        // Check if the inner key exists, if not create a new List<Float>
        innerMap.computeIfAbsent(innerKey, k -> new ArrayList<>());

        // Get the list associated with the inner key and add the value
        List<Float> valueList = innerMap.get(innerKey);
        valueList.add(value);
    }

    public static void calculateMedians() {
        // calculate RT median shifts
        for (String glycanName : allShiftedRTs.keySet()) {
            HashMap<String, List<Float>> innerMap = allShiftedRTs.get(glycanName);
            HashMap<String, Float> playerMedians = new HashMap<>();

            for (String game : innerMap.keySet()) {
                List<Float> scores = innerMap.get(game);
                Float median = calculateMedian(scores);
                playerMedians.put(game, median);
            }
            shiftedRTs.put(glycanName, playerMedians);
        }
        // calculate IM median shifts
        for (String glycanName : allShiftedIMs.keySet()) {
            HashMap<String, List<Float>> innerMap = allShiftedIMs.get(glycanName);
            HashMap<String, Float> playerMedians = new HashMap<>();

            for (String game : innerMap.keySet()) {
                List<Float> scores = innerMap.get(game);
                Float median = calculateMedian(scores);
                playerMedians.put(game, median);
            }
            shiftedIMs.put(glycanName, playerMedians);
        }
        // calculate RT median shifts of monosaccharides
        for (String residual : allShiftedMonoRTs.keySet()) {
            monosaccharideRtShifts.put(residual,calculateMedian(allShiftedMonoRTs.get(residual)));
        }
        // calculate IM median shifts of monosaccharides
        for (String residual : allShiftedMonoIMs.keySet()) {
            monosaccharideImShifts.put(residual,calculateMedian(allShiftedMonoIMs.get(residual)));
        }
    }

    public static Float calculateMedian(List<Float> scores) {
        if (scores.isEmpty()) {
            return null; // or throw an exception
        }
        Collections.sort(scores);
        int size = scores.size();
        if (size % 2 == 0) {
            // Even number of elements
            return (scores.get(size / 2 - 1) + scores.get(size / 2)) / 2;
        } else {
            // Odd number of elements
            return scores.get(size / 2);
        }
    }


    Float getShiftedRT(Glycan oldGlycan, Glycan newGlycan, double oldRT) {
        Float shiftedRT;
        if (oldGlycan.mass < newGlycan.mass) {
            if (shiftedRTs.containsKey(oldGlycan.toString()) && shiftedRTs.get(oldGlycan.toString()).containsKey(newGlycan.toString())) {
                shiftedRT = shiftedRTs.get(oldGlycan.toString()).get(newGlycan.toString());
            } else {
                String shiftedMonosaccharide = getShiftedSaccharides(oldGlycan, newGlycan);
                shiftedRT = monosaccharideRtShifts.getOrDefault(shiftedMonosaccharide, (float) 0);
            }
        } else {
            if (shiftedRTs.containsKey(oldGlycan.toString()) && shiftedRTs.get(oldGlycan.toString()).containsKey(newGlycan.toString())) {
                shiftedRT = -shiftedRTs.get(oldGlycan.toString()).get(newGlycan.toString());
            } else {
                String shiftedMonosaccharide = getShiftedSaccharides(newGlycan, oldGlycan);
                shiftedRT = -monosaccharideRtShifts.getOrDefault(shiftedMonosaccharide, (float) 0);
            }
        }
        return (float) (oldRT + shiftedRT);
    }

    Float getShiftedIM(Glycan oldGlycan, Glycan newGlycan, double oldIM) {
        Float shiftedIM;
        if (oldGlycan.mass < newGlycan.mass) {
            if (shiftedIMs.containsKey(oldGlycan.toString()) && shiftedIMs.get(oldGlycan.toString()).containsKey(newGlycan.toString())) {
                shiftedIM = shiftedIMs.get(oldGlycan.toString()).get(newGlycan.toString());
            } else {
                String shiftedMonosaccharide = getShiftedSaccharides(oldGlycan, newGlycan);
                shiftedIM = monosaccharideImShifts.getOrDefault(shiftedMonosaccharide, (float) 0);
            }
        } else {
            if (shiftedIMs.containsKey(oldGlycan.toString()) && shiftedIMs.get(oldGlycan.toString()).containsKey(newGlycan.toString())) {
                shiftedIM = -shiftedIMs.get(oldGlycan.toString()).get(newGlycan.toString());
            } else {
                String shiftedMonosaccharide = getShiftedMonosaccharide(newGlycan, oldGlycan);
                shiftedIM = -monosaccharideImShifts.getOrDefault(shiftedMonosaccharide, (float) 0);
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
        if (allShiftedRTs.containsKey(oldGlycan.toString()) && allShiftedRTs.get(oldGlycan.toString()).containsKey(newGlycan.toString())) {
            return allShiftedRTs.get(oldGlycan.toString()).get(newGlycan.toString()).size();
        }
        else {
            return 0;
        }
    }
}
