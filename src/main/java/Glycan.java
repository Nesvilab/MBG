import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Glycan {
    public TreeMap<GlycanResidue, Integer> composition;
    public double mass;

    private static final Pattern numberPattern = Pattern.compile("[0-9]");
    private static final Pattern letterPattern = Pattern.compile("[a-zA-Z]+");
    private static final Pattern byonicPattern = Pattern.compile("[A-Za-z]+\\([0-9]+\\)");      // e.g., HexNAc(1)Hex(1)

    public Glycan(String compString, Params analysis) {
        composition = parseGlycan(compString, analysis);
        mass = computeMonoisotopicMass();
    }


    /**
     * Parse glycans of varying formats given an input Regex Pattern.Matcher
     * @return
     */
    private TreeMap<GlycanResidue, Integer> parseGlycan(String line, Params mbgAnalysis) {
        TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();
        Matcher glycanMatch = byonicPattern.matcher(line);
        while(glycanMatch.find()) {
            String glycanToken = glycanMatch.group();
            Matcher glycanName = letterPattern.matcher(glycanToken);
            GlycanResidue residue;
            if (glycanName.find()) {
                String glycan = glycanName.group();
                residue = mbgAnalysis.findResidueName(glycan);
                if (residue == null) {
                    System.out.printf("Warning: glycan type %s in glycan %s not recognized. This glycan will not be included in search! Please add %s to the glycan_residues.tsv table.%n", glycanName, line, glycanName);
                    return new TreeMap<>();
                }
            } else {
                System.out.printf("Warning: glycan type %s in glycan %s not recognized. This glycan will not be included in search! Please add %s to the glycan_residues.tsv table.%n", glycanName, line, glycanName);
                return new TreeMap<>();
            }
            Matcher glycanCount = numberPattern.matcher(glycanToken);
            if (glycanCount.find()) {
                int count = Integer.parseInt(glycanCount.group());
                glycanComp.put(residue, count);
            } else {
                System.out.printf("Warning: glycan count not found in glycan %s and was ignored%n", line);
            }
        }
        return glycanComp;
    }

    /**
     * Compute exact mass of a given composition
     * @return monoisotopic mass
     */
    public double computeMonoisotopicMass() {
        double mass = 0;
        for (Map.Entry<GlycanResidue, Integer> glycanEntry : composition.entrySet()) {
            // mass = residue mass * residue count
            mass += (glycanEntry.getKey().mass * glycanEntry.getValue());
        }
        return mass;
    }

    public String toCompString() {
        StringBuilder stringBuilder = new StringBuilder();
        for (GlycanResidue residueKey : composition.keySet()) {
            if (composition.getOrDefault(residueKey, 0) > 0) {
                stringBuilder.append(String.format("%s(%d)", residueKey.name, composition.get(residueKey)));
            }
        }
        return stringBuilder.toString();
    }

}
