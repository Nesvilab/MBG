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
import java.util.StringTokenizer;

import static java.lang.System.out;

public class Params {

    public ArrayList<Glycan> glycanDatabase;
    public int numGlycanResidues;
    public double maxGlycanQval = 0.01;
    public int minPSMs = 5;
    public int minGlycans = 2;
    public static float rtTol = 0.4F;
    public static float imTol = 0.05F;
    public static float fdr = 0.01F;
    public static int mztol = 5;
    private static final String defaultResiduePath = "glycan_residues.txt";
    private static final String defaultModPath = "glycan_mods.txt";
    private static Path defaultPTMSconfPath;
    public static boolean noPASEF = false;
    public Integer threads;
    public Path psmPath;
    public Path glycanDBpath;
    public Path residueDBpath;
    public Path glycanModPath;
    public String msFilePath;
    public HashMap<String, GlycanResidue> glycanResidues;
    public HashMap<String, GlycanMod> glycanMods;
    public HashMap<String, GlycanResidue> toAddResidues;
    private String[] residualStrings;

    /**
     * Parse command line args to parameters for matching (inference) mode
     *
     * @param args command line args
     */
    public Params(String[] args) throws Exception {
        for (int i = 1; i < args.length; i += 2) {
            switch (args[i]) {
                case "--psm":
                    psmPath = checkFilePath(args[i + 1].trim(), "psm");
                    defaultPTMSconfPath = psmPath.getParent().resolve("shepherd.config");
                    break;
                case "--glycandb":
                    glycanDBpath = checkFilePath(args[i + 1].trim(), "glycan database");
                    break;
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
        if (glycanModPath == null) {
            glycanMods = GlycanParser.parseGlycoModsDBStream(Params.class.getResourceAsStream(defaultModPath), glycanResidues.size(), glycanResidues);
        } else {
            glycanMods = GlycanParser.parseGlycoModsDB(glycanModPath.toString(), glycanResidues.size(), glycanResidues);
        }
        for (Map.Entry<String, GlycanMod> stringGlycanModEntry : glycanMods.entrySet()) {
            glycanResidues.put(stringGlycanModEntry.getKey(), stringGlycanModEntry.getValue());
        }
        toAddResidues = getToAddResidues(residualStrings, glycanResidues);
        glycanDatabase = GlycanParser.loadGlycansFromText(glycanDBpath.toString(), DatabaseType.byonic, glycanResidues);
        // parse ptmsherped config file to get the path of spectra
        msFilePath = getMSFilePath(defaultPTMSconfPath);
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

    private String getMSFilePath(Path shepConfigPath) throws IOException {
        HashMap<String, String> pathMap = new HashMap<>();
        BufferedReader in = new BufferedReader(new FileReader(shepConfigPath.toFile()));
        String cline;
        while ((cline = in.readLine()) != null) {
            int comments = cline.indexOf("//");
            if (comments >= 0)
                cline = cline.substring(0, comments);
            cline = cline.trim();
            if (!cline.contains("="))
                continue;
            String key = cline.substring(0, cline.indexOf("=")).trim();
            String value = cline.substring(cline.indexOf("=") + 1).trim();
            if (key.equals("dataset")) {
                StringTokenizer st = new StringTokenizer(value);
                String dsName = st.nextToken();
                String tsvTxt = st.nextToken();
                String mzPath = st.nextToken();
                pathMap.put(tsvTxt, mzPath);
            } else if (key.equals("threads")) {
                threads = Integer.parseInt(value.trim());
            }
        }
        in.close();
        // get current ms file path
        return pathMap.get(psmPath.toString());
    }
}