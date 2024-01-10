import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

import static java.lang.System.out;

public class Params {

    public ArrayList<Glycan> glycanDatabase;
    private ArrayList<GlycanResidue> glycanResidues;
    public double maxGlycanQval = 0.01;
    public int minPSMs = 5;
    public int minGlycans = 2;
    private static final String defaultResiduePath = "glycan_residues.txt";
    private static final String defaultGlycanPath = "glycan_database.txt";
    public Path psmPath;
    public Path glycanDBpath;
    public Path residueDBpath;

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
        glycanResidues = parseGlycoResiduesDB(residueDBpath);
        glycanDatabase = parseGlycanDB(glycanDBpath);
    }

    private Path checkFilePath(String inputPath, String pathName) {
        Path path = null;
        try {
            path = Paths.get(inputPath.replaceAll("['\"]", ""));
        } catch (Exception e) {
            e.printStackTrace();
            out.printf("Error parsing %s file path: [%s]", pathName, inputPath);
            System.exit(1);
        }
        if (!Files.exists(path)) {
            out.printf("%s file path does not exist: [%s]", pathName, inputPath);
            System.exit(1);
        }
        return path;
    }

    private ArrayList<Glycan> parseGlycanDB(Path glycanDBpath) {
        ArrayList<Glycan> glycanDB = new ArrayList<>();
        try {
            BufferedReader in;
            if (glycanDBpath == null) {
                // no glycan database provided - fall back to default glycan list
                in = new BufferedReader(new InputStreamReader(Params.class.getResourceAsStream(defaultGlycanPath)));
            } else {
                in = new BufferedReader(new FileReader(glycanDBpath.toFile()));
            }

            HashMap<String, Boolean> glycansInDB = new HashMap<>();
            String readline;
            while ((readline = in.readLine()) != null) {
                Glycan glycan = new Glycan(readline, this);
                if (glycan.composition.size() == 0) {
                    continue;
                }
                String compStr = glycan.toString();
                // prevent duplicated glycans
                if (!glycansInDB.containsKey(compStr)) {
                    glycansInDB.put(compStr, true);
                    glycanDB.add(glycan);
                }
            }
        } catch (FileNotFoundException e) {
            System.out.printf("Glycan database file not found: [%s]%n", glycanDBpath);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("IO Exception while reading database file");
        }
        return glycanDB;
    }

    private ArrayList<GlycanResidue> parseGlycoResiduesDB(Path glycanResiduesPath) {
        ArrayList<GlycanResidue> residues = new ArrayList<>();
        BufferedReader in;
        try {
            // no glycan database provided - fall back to default glycan list in PeakAnnotator
            if (glycanResiduesPath == null) {
                in = new BufferedReader(new InputStreamReader(Params.class.getResourceAsStream(defaultResiduePath)));
            } else {
                in = new BufferedReader(new FileReader(glycanResiduesPath.toFile()));
            }
            String line;
            while ((line = in.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
                GlycanResidue residue = GlycanResidue.parseResidue(line);
                residues.add(residue);
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            System.out.printf("Error parsing input monosaccharide table %s%n", glycanResiduesPath);
        }
        return residues;
    }

    /**
     * Search all GlycanResidue names and alternate names to try to match this input to one of them.
     * @param inputName string
     * @return
     */
    public GlycanResidue findResidueName(String inputName) {
        for (GlycanResidue residue: glycanResidues) {
            if (residue.name.equalsIgnoreCase(inputName)) {
                return residue;
            }
        }
        // no exact match: check alternate names
        for (GlycanResidue residue: glycanResidues) {
            for (String altName : residue.alternateNames) {
                if (altName.equalsIgnoreCase(inputName)) {
                    return residue;
                }
            }
        }
        // no match - return null and warn user
        return null;
    }


}
