package com.mbg;

import edu.umich.andykong.ptmshepherd.core.Spectrum;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanParser;
import umich.ms.glyco.GlycanResidue;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PSM implements Comparable<PSM>{

    double rt;
    double apexRT;
    int scanNum;
    String peptide;
    Glycan glycan;
    String protein;
    int proteinGlycositeIndex;
    int peptideGlycositeIndex;
    String[] splits;
    char glycoSiteAA;
    double calcMass;
    double calcMz;
    int charge;
    double glycoQval;
    double deltaMass;
    boolean inferredGlycoform;
    boolean hasMS2Spectrum;
    boolean hasY1Ion;
    int rank;
    double intensity;
    public String rawfile;

    private static final Pattern modMassPattern = Pattern.compile("[\\d]+[A-Z]\\(([\\d.]+)\\)");
    private static final Pattern modSitePattern = Pattern.compile("([\\d]+)[A-Z]\\([\\d.]+\\)");
    private static final Pattern rankPattern = Pattern.compile("_rank([\\d]+)\\.pep\\.xml");


    /**
     * Initialize a PSM object from a line of the psm.tsv table
     * @param psmLine
     */
    public PSM(String psmLine, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params, boolean isInferred, boolean hasMS2, boolean hasY1, int addCols) {
        splits = psmLine.split("\t");

        // expand to hold new columns if needed
        if (addCols > 0) {
            int add = addCols;
            ArrayList<String> newSplits = new ArrayList<>(Arrays.asList(splits));
            while (add > 0) {
                newSplits.add(psmCols.get(PSM_File.PSMattr.glyco_q_val) + 1, "");
                add--;
            }
            splits = newSplits.toArray(new String[0]);
            inferredGlycoform = isInferred;
            hasMS2Spectrum = hasMS2;
            hasY1Ion = hasY1;
        } else {
            if (!isInferred) {
                // if this is a re-analysis (i.e., not adding new columns) and this is not an inferred PSM from this analysis, read its inferred status from the previous result
                inferredGlycoform = Boolean.parseBoolean(splits[psmCols.get(PSM_File.PSMattr.inferred_glycoform)]);
                hasMS2Spectrum = Boolean.parseBoolean(splits[psmCols.get(PSM_File.PSMattr.has_ms2_spec)]);
                hasY1Ion = Boolean.parseBoolean(splits[psmCols.get(PSM_File.PSMattr.has_Y1)]);
            } else {
                inferredGlycoform = isInferred;
                hasMS2Spectrum = hasMS2;
                hasY1Ion = hasY1;
            }
        }

        // update new fields and information that is needed for all PSMs (including non-glyco)
        splits[psmCols.get(PSM_File.PSMattr.inferred_glycoform)] = String.valueOf(inferredGlycoform);
        splits[psmCols.get(PSM_File.PSMattr.has_ms2_spec)] = String.valueOf(hasMS2Spectrum);
        splits[psmCols.get(PSM_File.PSMattr.has_Y1)] = String.valueOf(hasY1Ion);
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
        if (glycan.length() == 0) {
            // non-glyco PSM
            this.glycan = null;
            return;
        }
        try {
            glycoQval = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.glyco_q_val)]);
        } catch (NumberFormatException ex) {
            glycoQval = 1;
        }
        if (glycoQval > params.maxGlycanQval) {
            this.glycan = null;
            return;
        }
        this.glycan = GlycanParser.parseGlycanString(glycan, params.glycanResidues);
        // parse other needed PSM info
        peptide = splits[psmCols.get(PSM_File.PSMattr.peptide)];
        protein = splits[psmCols.get(PSM_File.PSMattr.protein)];
        rt = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.RT)]);
        int proteinStart = Integer.parseInt(splits[psmCols.get(PSM_File.PSMattr.protein_start)]);
        String assignedMods = splits[psmCols.get(PSM_File.PSMattr.assigned_mods)];
        peptideGlycositeIndex = getGlycoSiteIndex(assignedMods);
        glycoSiteAA = peptide.charAt(peptideGlycositeIndex - 1);
        proteinGlycositeIndex = proteinStart + peptideGlycositeIndex;
        calcMass = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.calc_pep_mass)]);
        calcMz = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.calc_pep_mz)]);
        charge = Integer.parseInt(splits[psmCols.get(PSM_File.PSMattr.charge)]);
        deltaMass = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.delta_mass)]);
        intensity = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.intensity)]);
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
     * Generate a 'dummy' PSM for a newly matched glycan from an existing base PSM. Edit the base PSM to
     * convert to the new glycan
     *
     * @param existingPSM
     * @param newGlycan
     * @param psmCols
     * @param params
     */
    public static List<PSM> makeDummyPSM(PSM existingPSM, Glycan newGlycan, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params, int rankCounter) {
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
                if (site == existingPSM.peptideGlycositeIndex) {
                    // edit this mod to have the mass of the new glycan
                    newAssignedMods.add(String.format("%d%s(%.4f)", site, existingPSM.peptide.charAt(site - 1), newGlycan.mass));
                } else {
                    newAssignedMods.add(mod);
                }
            }
        }
        newSplits[psmCols.get(PSM_File.PSMattr.assigned_mods)] = String.join(",", newAssignedMods);

        // mod peptide
        String oldGlycEntry = String.format("[%d]", Math.round(existingPSM.glycan.mass + AA_masses.aa_masses.get(existingPSM.peptide.charAt(existingPSM.peptideGlycositeIndex - 1))));
        String newGlycEntry = String.format("[%d]", Math.round(newGlycan.mass + AA_masses.aa_masses.get(existingPSM.peptide.charAt(existingPSM.peptideGlycositeIndex - 1))));
        String newModPep = existingPSM.splits[psmCols.get(PSM_File.PSMattr.mod_peptide)].replace(oldGlycEntry, newGlycEntry);
        newSplits[psmCols.get(PSM_File.PSMattr.mod_peptide)] = newModPep;

        // update calculated masses and delta mass for new glycan
        newSplits[psmCols.get(PSM_File.PSMattr.calc_pep_mass)] = String.format("%.4f", existingPSM.calcMass + newGlycan.mass - existingPSM.glycan.mass);
        newSplits[psmCols.get(PSM_File.PSMattr.calc_pep_mz)] = String.format("%.4f", existingPSM.calcMz + ((newGlycan.mass - existingPSM.glycan.mass) / (double) existingPSM.charge));
        newSplits[psmCols.get(PSM_File.PSMattr.delta_mass)] = String.format("%.4f", existingPSM.deltaMass - existingPSM.glycan.mass + newGlycan.mass);

        // RT
        // todo: compute some RT adjustment...
        newSplits[psmCols.get(PSM_File.PSMattr.RT)] = String.format("%.3f", getShiftedRT(existingPSM.glycan, newGlycan, existingPSM.rt));
        double newRT = Double.parseDouble(newSplits[psmCols.get(PSM_File.PSMattr.RT)]);
        double newMZ = Double.parseDouble(newSplits[psmCols.get(PSM_File.PSMattr.calc_pep_mz)]);

        List<PSM> addedPSMs = new ArrayList<>();

        // match against MS1, get intensity, rt duration and average mass
        Peak peakMatcher = new Peak();
        peakMatcher.matchMS1((float) newMZ, (float) newRT / 60, 0, existingPSM.charge, 0);
        if (peakMatcher.quantifiedEntry123 == null) {
            addedPSMs.add(new PSM(String.join("\t", newSplits), psmCols, params, true, false, false, 0));
            return addedPSMs;
        }
        // to find MS2
        List<Spectrum> ms2Candidates = MSMS.spectraSearcher(peakMatcher.quantifiedEntry123.rtStart, peakMatcher.quantifiedEntry123.rtEnd, (float) (newMZ-1), (float) (newMZ+1));
        // update intensity
        newSplits[psmCols.get(PSM_File.PSMattr.intensity)] = String.valueOf(peakMatcher.quantifiedEntry123.intensity);

        // todo update spectra information if ms2 candidates found

        if (!ms2Candidates.isEmpty()) {
            // match Y1
            List<Boolean> y1booleans;
            y1booleans = MSMS.y1Matcher(ms2Candidates, existingPSM, params);
            for (int i = 0; i < ms2Candidates.size(); i++) {
                newSplits[psmCols.get(PSM_File.PSMattr.spectrum)] = String.join(".", new String[]{existingPSM.rawfile, String.valueOf(ms2Candidates.get(i).scanNum),
                        String.valueOf(ms2Candidates.get(i).scanNum), String.valueOf(ms2Candidates.get(i).charge)});
                PSM addedPSM = new PSM(String.join("\t", newSplits), psmCols, params, true, true, y1booleans.get(i), 0);
                addedPSM.updateRank(existingPSM.rank + rankCounter, psmCols);
                addedPSMs.add(addedPSM);
            }

        } else {
            PSM addedPSM = new PSM(String.join("\t", newSplits), psmCols, params, true, false, false, 0);
            addedPSM.updateRank(existingPSM.rank + rankCounter, psmCols);
            addedPSMs.add(addedPSM);
        }


        // Y1 matcher test

//        List<Spectrum> testCandidate = MSMS.spectraSearcher(existingPSM.rt-40, existingPSM.rt+40, (float) (existingPSM.calcMz-1), (float) (existingPSM.calcMz+1));
//        List<Spectrum> testY1 = new ArrayList<>();
//        if (!testCandidate.isEmpty()) {
//            testY1 = MSMS.y1Matcher(testCandidate, existingPSM, params);
//        }
        // construct the new PSM

//        PSM addedPSM = new PSM(String.join("\t", newSplits), psmCols, params, true, !ms2Candidates.isEmpty(), !y1booleans.isEmpty(), 0);
//        PSM addedPSM = new PSM(String.join("\t", existingPSM.splits), psmCols, params, true, !testCandidate.isEmpty(), !testY1.isEmpty(), 0);
        // update rank
//        addedPSM.updateRank(existingPSM.rank + rankCounter, psmCols);
        return addedPSMs;
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

    private static double getShiftedRT(Glycan oldGlycan, Glycan newGlycan, double oldRT) {
        double newRT;
        if (glycanMassShift.containsKey(newGlycan.name)) {
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

    private static int glycanComparator(HashMap <GlycanResidue, Integer> glycan1, HashMap <GlycanResidue, Integer> glycan2) {
        int difference = 0;
        for (GlycanResidue key : glycan1.keySet()) {
            difference += Math.abs(glycan1.getOrDefault(key, 0) - glycan2.getOrDefault(key, 0));
        }
        return difference;
    }

    /**
     * Given a glycan, adding NH3 to generate a new glycan
     *
     * @return
     */
    public static Glycan addResidual(Glycan baseGlycan, GlycanResidue mod) {
        HashMap<GlycanResidue, Integer> newComposition = new HashMap(baseGlycan.composition);
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
                // add
                HashMap<GlycanResidue, Integer> newCompositionAdd = new HashMap(psm.glycan.composition);
                newCompositionAdd.put(residual, count + 1);
                dummy_glycans.add(new Glycan(newCompositionAdd));
                // loss
                HashMap<GlycanResidue, Integer> newCompositionLoss = new HashMap(psm.glycan.composition);
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
                HashMap<GlycanResidue, Integer> newCompositionAdd = new HashMap(psm.glycan.composition);
                newCompositionAdd.put(residual, 1);
                dummy_glycans.add(new Glycan(newCompositionAdd));
            }
        }
        return dummy_glycans;
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
        return -1;
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
