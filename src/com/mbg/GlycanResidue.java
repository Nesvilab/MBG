package com.mbg;

public class GlycanResidue implements Comparable<GlycanResidue> {
    public String name;
    public double mass;
    String[] alternateNames;
    private final int printIndex;

    public GlycanResidue(String name, double mass, String[] altNames, int printIndex) {
        this.name = name;
        this.mass = mass;
        this.alternateNames = altNames;
        this.printIndex = printIndex;
    }

    /**
     * Parse a residue from a FragPipe-formatted residues list. The file includes extra information for other
     * tools that doesn't need to be parsed here.
     * @param line
     * @return
     */
    public static GlycanResidue parseResidue(String line, int printIndex) {
        String[] splits = line.replace("\"", "").split("\t");
        double mass = Double.parseDouble(splits[1]);

        String[] altNames = splits[4].split(",");
        return new GlycanResidue(splits[0],
                mass,
                altNames,
                printIndex
        );
    }

    public String toString() {
        return name;
    }

    @Override
    public int compareTo(GlycanResidue o) {
        return Integer.compare(printIndex, o.printIndex);
    }
}
