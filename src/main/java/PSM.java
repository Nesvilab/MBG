import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PSM {

    double rt;
    String peptide;
    Glycan glycanComp;
    String protein;
    int glycoSiteIndex;
    String[] splits;

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
            glycanComp = null;
            return;
        }
        glycanComp = new Glycan(glycan, params);
        peptide = splits[psmCols.get(PSM_File.PSMattr.peptide)];
        protein = splits[psmCols.get(PSM_File.PSMattr.protein)];
        rt = Double.parseDouble(splits[psmCols.get(PSM_File.PSMattr.RT)]);
        int proteinStart = Integer.parseInt(splits[psmCols.get(PSM_File.PSMattr.protein_start)]);
        String assignedMods = splits[psmCols.get(PSM_File.PSMattr.assigned_mods)];
        glycoSiteIndex = proteinStart + getGlycoSiteIndex(assignedMods);

    }

    /**
     * Generate a 'dummy' PSM for a newly matched glycan from an existing base PSM. Edit the base PSM to
     * convert to the new glycan
     * @param existingPSM
     * @param newGlycan
     * @param psmCols
     * @param params
     */
    public PSM(PSM existingPSM, Glycan newGlycan, HashMap<PSM_File.PSMattr, Integer> psmCols, Params params) {

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
                if (Math.abs(mass - glycanComp.mass) < 0.0001) {
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
        return String.format("%s_%d", protein, glycoSiteIndex);
    }

    public int getScanNum(HashMap<PSM_File.PSMattr, Integer> psmColumns) {
        String specStr = splits[psmColumns.get(PSM_File.PSMattr.spectrum)];
        String[] specSplits = specStr.split("\\.");
        return Integer.parseInt(specSplits[1]);
    }

}
