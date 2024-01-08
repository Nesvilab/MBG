import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PSM {

    double rt;
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

    private static final Pattern modMassPattern = Pattern.compile("[\\d]+[A-Z]\\(([\\d.]+)\\)");
    private static final Pattern modSitePattern = Pattern.compile("([\\d]+)[A-Z]\\([\\d.]+\\)");


    /**
     * Initialize a PSM object from a line of the psm.tsv table
     * @param psmLine
     */
    public PSM(String psmLine, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params) {
        splits = psmLine.split("\t");
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
     * Generate a 'dummy' PSM for a newly matched glycan from an existing base PSM. Edit the base PSM to
     * convert to the new glycan
     * @param existingPSM
     * @param newGlycan
     * @param psmCols
     * @param params
     */
    public static PSM makeDummyPSM(PSM existingPSM, Glycan newGlycan, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params) {
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
        return new PSM(String.join("\t", newSplits), psmCols, params);
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

}
