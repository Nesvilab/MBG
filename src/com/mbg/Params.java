package com.mbg;

import umich.ms.glyco.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static java.lang.System.out;

public class Params {

    public double maxGlycanQval = 0.01;
    public int minPSMs = 5;
    public int minGlycans = 2;
    public float rtTol = 0.4F;
    public float imTol = 0.05F;
    public float fdr = 0.05F;
    public int mztol = 5;
    private static final String defaultResiduePath = "glycan_residues.txt";
    private static final String defaultModPath = "glycan_mods.txt";
    public boolean noPASEF = false;
    public boolean runTMT = false;
    public int expandDB = 0;
    public int minIsotopesIonQuant = 2;
    public int minScansIonQuant = 3;
    public int numThreads = 4;
    public Path psmPath;
    public Path residueDBpath;
    public Path glycanModPath;
    public HashMap<String, String> msFilePathsMap;
    public HashMap<String, GlycanResidue> glycanResidues;
    public HashMap<String, GlycanMod> glycanMods;
    public HashMap<String, GlycanResidue> toAddResidues;
    public ArrayList<Glycan> toAddGlycans;
    private String[] residualStrings;
    private Path manifestPath;

    /**
     * Parse command line args to parameters for matching (inference) mode
     *
     * @param args command line args
     */
    public Params(String[] args) {
        for (int i = 1; i < args.length; i += 2) {
            switch (args[i]) {
                case "--psm":
                    psmPath = checkFilePath(args[i + 1].trim(), "psm");
                    break;
                case "--manifest":
                    manifestPath = checkFilePath(args[i + 1].trim(), "manifest");
                    msFilePathsMap = readManifest();
                case "--residuedb":
                    residueDBpath = checkFilePath(args[i + 1].trim(), "residue database");
                    break;
                case "--glycanmoddb":
                    glycanModPath = checkFilePath(args[i + 1].trim(), "glycan mods database");
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
                case "--minglycans":
                    try {
                        minGlycans = Integer.parseInt(args[i + 1].trim());
                    } catch (NumberFormatException ex) {
                        out.printf("Error parsing min Glycans argument %s: not a number", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--rttol":
                    try {
                        rtTol = Float.parseFloat(args[i + 1].trim());
                    } catch (NumberFormatException ex) {
                        out.printf("Error parsing rt tolerance argument %s: not a number", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--imtol":
                    try {
                        imTol = Float.parseFloat(args[i + 1].trim());
                    } catch (NumberFormatException ex) {
                        out.printf("Error parsing im tolerance argument %s: not a number", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--fdr":
                    try {
                        fdr = Float.parseFloat(args[i + 1].trim());
                    } catch (NumberFormatException ex) {
                        out.printf("Error parsing fdr argument %s: not a number", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--nopasef":
                    try {
                        noPASEF = Boolean.parseBoolean(args[i + 1].trim());
                    } catch (Exception e) {
                        out.printf("Error parsing pasef argument %s: neither true nor false", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--runtmt":
                    try {
                        runTMT = Boolean.parseBoolean(args[i + 1].trim());
                    } catch (Exception e) {
                        out.printf("Error parsing runtmt argument %s: neither true nor false", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--expanddb":
                    try {
                        expandDB = Integer.parseInt(args[i + 1].trim());
                    } catch (Exception e) {
                        out.printf("Error parsing expanddb argument %s: neither true nor false", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--minscans":
                    try {
                        minScansIonQuant = Integer.parseInt(args[i + 1].trim());
                    } catch (Exception e) {
                        out.printf("Error parsing minscans argument %s: must be integer", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--minisotopes":
                    try {
                        minIsotopesIonQuant = Integer.parseInt(args[i + 1].trim());
                    } catch (Exception e) {
                        out.printf("Error parsing minisotopes argument %s: must be integer", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--numthreads":
                    try {
                        numThreads = Integer.parseInt(args[i + 1].trim());
                    } catch (Exception e) {
                        out.printf("Error parsing numthreads argument %s: must be integer", args[i + 1]);
                        System.exit(1);
                    }
                    break;
                case "--toaddresiduals":
                    try {
                        residualStrings = args[i + 1].trim().split(",");
                    } catch (Exception e) {
                        out.printf("Error parsing toaddresiduals argument %s", args[i + 1]);
                        System.exit(1);
                    }
                    break;
            }
        }
        if (psmPath == null) {
            out.println("Error: PSM file path not specified!");
            System.exit(1);
        }
        // use default residues if none provided
        if (residueDBpath == null) {
            glycanResidues = GlycanParser.parseGlycoResiduesDBStream(Params.class.getResourceAsStream(defaultResiduePath));
        } else {
            glycanResidues = GlycanParser.parseGlycoResiduesDB(residueDBpath.toString());
        }
        if (glycanModPath == null) {
            glycanMods = GlycanParser.parseGlycoModsDBStream(Params.class.getResourceAsStream(defaultModPath), glycanResidues.size(), glycanResidues);
        } else {
            glycanMods = GlycanParser.parseGlycoModsDB(glycanModPath.toString(), glycanResidues.size(), glycanResidues);
        }
        for (Map.Entry<String, GlycanMod> stringGlycanModEntry : glycanMods.entrySet()) {
            glycanResidues.put(stringGlycanModEntry.getKey(), stringGlycanModEntry.getValue());
        }
//        toAddResidues = getToAddResidues(residualStrings, glycanResidues);
        toAddGlycans = getToAddGlycans(residualStrings, glycanResidues);
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

    private HashMap<String, GlycanResidue> getToAddResidues(String[] residualStrings, HashMap<String, GlycanResidue> residues) {
        HashMap<String, GlycanResidue> toAddResidues = new HashMap<>();
        for (String residue : residualStrings) {
            toAddResidues.put(residue, residues.get(residue));
        }
        return toAddResidues;
    }

    private ArrayList<Glycan> getToAddGlycans(String[] residualStrings, HashMap<String, GlycanResidue> residues) {
        ArrayList<Glycan> toAddGlycans = new ArrayList<>();
        for (String residue : residualStrings) {
            if (residue.contains("(")) {
                // treat as a glycan
                toAddGlycans.add(GlycanParser.parseGlycanString(residue, residues));
            } else {
            HashMap<GlycanResidue, Integer> cComp = new HashMap<>();
            cComp.put(residues.get(residue), 1);
            toAddGlycans.add(new Glycan(cComp));
            }
        }
        return toAddGlycans;
    }

    private HashMap<String, String> readManifest() {
        HashMap<String, String> msFilesMap = new HashMap<>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(manifestPath.toFile()));
            String line;
            while ((line = reader.readLine()) != null) {
                String[] splits = line.split("\t");
                Path path = Paths.get(splits[0]);
                String fileName = path.getFileName().toString();
                int dotIndex = fileName.lastIndexOf('.');
                msFilesMap.put(fileName.substring(0, dotIndex), splits[0].replace("\\", "\\\\"));
            }
            reader.close();
        } catch (IOException ex) {
            out.printf("Error reading FragPipe manifest file %s\n", manifestPath);
            ex.printStackTrace();
            System.exit(1);
        }
        return msFilesMap;
    }
}