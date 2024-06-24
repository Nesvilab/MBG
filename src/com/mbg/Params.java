package com.mbg;

import umich.ms.glyco.DatabaseType;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanParser;
import umich.ms.glyco.GlycanResidue;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import static java.lang.System.out;

public class Params {

    public ArrayList<Glycan> glycanDatabase;
    public int numGlycanResidues;
    public double maxGlycanQval = 0.01;
    public int minPSMs = 5;
    public int minGlycans = 2;
    private static final String defaultResiduePath = "glycan_residues.txt";
    public Path psmPath;
    public Path glycanDBpath;
    public Path residueDBpath;
    public HashMap<String, GlycanResidue> glycanResidues;

    /**
     * Parse command line args to parameters for matching (inference) mode
     * @param args command line args
     */
    public Params(String[] args) {
        for (int i=1; i < args.length; i+=2) {
            switch (args[i]) {
                case "--psm":
                    psmPath = checkFilePath(args[i + 1].trim(), "psm");
                    break;
                case "--glycandb":
                    glycanDBpath = checkFilePath(args[i + 1].trim(), "glycan database");
                    break;
                case "--residuedb":
                    residueDBpath = checkFilePath(args[i + 1].trim(), "residue database");
                    break;
                case "--maxq":
                    try {
                        maxGlycanQval = Double.parseDouble(args[i + 1].trim());
                    } catch (NumberFormatException ex) {
                        out.printf("Error parsing max glycan q-value argument %s: not a number", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--minpsms":
                    try {
                        minPSMs = Integer.parseInt(args[i + 1].trim());
                    } catch (NumberFormatException ex) {
                        out.printf("Error parsing min PSMs argument %s: not a number", args[i + 1]);
                        System.exit(1);
                    }
                    break;
            }
        }
        if (psmPath == null) {
            out.println("Error: PSM file path not specified!");
            System.exit(1);
        }
        if (glycanDBpath == null) {
            out.println("Error: glycan database not specified!");
            System.exit(1);
        }
        // use default residues if none provided
        if (residueDBpath == null) {
            glycanResidues = GlycanParser.parseGlycoResiduesDBStream(Params.class.getResourceAsStream(defaultResiduePath));
        } else {
            glycanResidues = GlycanParser.parseGlycoResiduesDB(residueDBpath.toString());
        }
        glycanDatabase = GlycanParser.loadGlycansFromText(glycanDBpath.toString(), DatabaseType.byonic, glycanResidues);
    }

    public static Path checkFilePath(String inputPath, String pathName) {
        Path path = null;
        try {
            path = Paths.get(inputPath.replaceAll("['\"]", ""));
        } catch (Exception e) {
            e.printStackTrace();
            out.printf("Error parsing %s file path: [%s]\n", pathName, inputPath);
            System.exit(1);
        }
        if (!Files.exists(path)) {
            out.printf("%s file path does not exist: [%s]\n", pathName, inputPath);
            System.exit(1);
        }
        return path;
    }
}
