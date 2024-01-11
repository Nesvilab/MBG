package com.mbg;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PSM implements Comparable<PSM>{

    double rt;
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
    int rank;

    private static final Pattern modMassPattern = Pattern.compile("[\\d]+[A-Z]\\(([\\d.]+)\\)");
    private static final Pattern modSitePattern = Pattern.compile("([\\d]+)[A-Z]\\([\\d.]+\\)");
    private static final Pattern rankPattern = Pattern.compile("_rank([\\d]+)\\.pep\\.xml");


    /**
     * Initialize a PSM object from a line of the psm.tsv table
     * @param psmLine
     */
    public PSM(String psmLine, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params, boolean isInferred, int addCols) {
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
        }
        // update new fields and information that is needed for all PSMs (including non-glyco)
        inferredGlycoform = isInferred;
        splits[psmCols.get(PSM_File.PSMattr.inferred_glycoform)] = String.valueOf(inferredGlycoform);
        Matcher rankMatch = rankPattern.matcher(splits[psmCols.get(PSM_File.PSMattr.spectrum_file)]);
        if (rankMatch.find()) {
            rank = Integer.parseInt(rankMatch.group(1));
        } else {
            // no ranks provided, so need to make the base PSM rank 1 so added PSMs can be other ranks
            updateRank(1, psmCols);
        }
        scanNum = getScanNum(psmCols);

        // parse glycan and scores
        String glycan = splits[psmCols.get(PSM_File.PSMattr.glycan)];
        if (glycan.length() == 0) {
            // non-glyco PSM
            this.glycan = null;
            return;
        }
        glycoQval = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.glyco_q_val)]);
        if (glycoQval > params.maxGlycanQval) {
            this.glycan = null;
            return;
        }
        this.glycan = new Glycan(glycan, params);
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
     * @param existingPSM
     * @param newGlycan
     * @param psmCols
     * @param params
     */
    public static PSM makeDummyPSM(PSM existingPSM, Glycan newGlycan, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params, int rankCounter) {
        String[] newSplits = new String[existingPSM.splits.length];
        System.arraycopy(existingPSM.splits, 0, newSplits, 0, newSplits.length);

        // replace values with new glycan
        newSplits[psmCols.get(PSM_File.PSMattr.glycan)] = newGlycan.toOutputString();
        // assigned mods
        String oldAssignedMods = existingPSM.splits[psmCols.get(PSM_File.PSMattr.assigned_mods)];
        ArrayList<String> newAssignedMods = new ArrayList<>();
        for (String mod: oldAssignedMods.split(",")) {
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

        // construct the new PSM
        PSM addedPSM = new PSM(String.join("\t", newSplits), psmCols, params, true, 0);
        // update rank
        addedPSM.updateRank(existingPSM.rank + rankCounter, psmCols);
        return addedPSM;
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

    @Override
    public int compareTo(PSM o) {
        return Integer.compare(scanNum, o.scanNum);
    }
}
