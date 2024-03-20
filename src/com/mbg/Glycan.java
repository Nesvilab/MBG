package com.mbg;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Glycan {
    public HashMap<GlycanResidue, Integer> composition;
    public double mass;
    public String name;

    private static final Pattern numberPattern = Pattern.compile("[0-9]+");
    private static final Pattern letterPattern = Pattern.compile("[a-zA-Z]+");
    private static final Pattern byonicPattern = Pattern.compile("[A-Za-z]+\\([0-9]+\\)");      // e.g., HexNAc(1)Hex(1)

    public Glycan(String compString, Params params) {
        compString = compString.replaceAll(" % [\\d.]+", "");
        composition = parseGlycan(compString, params);
        mass = computeMonoisotopicMass();
        name = this.toString();
    }


    /**
     * Parse glycans of varying formats given an input Regex Pattern.Matcher
     * @return
     */
    private HashMap<GlycanResidue, Integer> parseGlycan(String line, Params params) {
        HashMap<GlycanResidue, Integer> glycanComp = new HashMap<>();
        Matcher glycanMatch = byonicPattern.matcher(line);
        while(glycanMatch.find()) {
            String glycanToken = glycanMatch.group();
            Matcher glycanName = letterPattern.matcher(glycanToken);
            GlycanResidue residue;
            if (glycanName.find()) {
                String glycan = glycanName.group();
                residue = params.findResidueName(glycan);
                if (residue == null) {
                    System.out.printf("Warning: glycan type %s in glycan %s not recognized. This glycan will not be included in search! Please add %s to the glycan_residues.tsv table.\n", glycanName, line, glycanName);
                    return new HashMap<>();
                }
            } else {
                System.out.printf("Warning: glycan type %s in glycan %s not recognized. This glycan will not be included in search! Please add %s to the glycan_residues.tsv table.\n", glycanName, line, glycanName);
                return new HashMap<>();
            }
            Matcher glycanCount = numberPattern.matcher(glycanToken);
            if (glycanCount.find()) {
                int count = Integer.parseInt(glycanCount.group());
                glycanComp.put(residue, count);
            } else {
                System.out.printf("Warning: glycan count not found in glycan %s and was ignored\n", line);
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

    @Override
    public String toString() {
        StringBuilder output = new StringBuilder();
        ArrayList<GlycanResidue> sortedKeys = new ArrayList<>(composition.keySet());
        Collections.sort(sortedKeys);
        for (GlycanResidue residue: sortedKeys) {
            output.append(residue.toString());
            output.append("(");
            output.append(composition.get(residue));
            output.append(")");
        }
        return output.toString();
    }

    /**
     * Generates Byonic format composition for writing to psm.tsv
     * @return
     */
    public String toOutputString() {
        return String.format("%s %% %.4f", this, mass);
    }

}
