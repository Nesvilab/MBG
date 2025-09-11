package com.mbg;

import umich.ms.glyco.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static java.lang.System.out;

/**
 * Handles parsing and storage of command-line parameters for the MBG application.
 * Validates inputs, loads necessary resource files (like glycan definitions),
 * and provides access to configuration settings for the analysis workflow.
 */
public class Params {

    // --- Input Data and Filtering Parameters ---
    /** Path to the input Peptide Spectrum Match (PSM) file. Required. (--psm) */
    public Path psmPath;
    /** Path to the FragPipe manifest file to map raw file names to full paths. Optional. (--manifest) */
    private Path manifestPath;
    /** Maximum Glycan q-value allowed for a PSM to be considered as confidently identified input. Default: 0.01. (--maxq) */
    public double maxGlycanQval = 0.01;
    /** Minimum number of observed PSMs required for a specific glycan+peptide combination across all combinations to be used for RT/IM shift calculation and inference. Default: 2. (--minpsms) */
    public int minPSMs = 2;
    /** Minimum number of *different* glycans observed for a specific peptide sequence to consider it for inference. Default: 2. (--minglycans) */
    public int minGlycans = 2;

    // --- Inference Tolerances and FDR ---
    /** Retention Time tolerance (in minutes) for matching inferred glycoforms. Default: 0.4. (--rttol) */
    public float rtTol = 0.4F;
    /** Ion Mobility tolerance (unitless or V·s·cm^-2) for matching inferred glycoforms. Default: 0.05. (--imtol) */
    public float imTol = 0.05F;
    /** False Discovery Rate threshold for filtering inferred glycoforms after scoring. Default: 0.05 (5%). (--fdr) */
    public float fdr = 0.05F;
    /** Mass tolerance (in ppm) for matching inferred glycoforms in MS1 spectra via IonQuant. Default: 5. (--mztol) */
    public int mztol = 5; // Note: Parameter name might be misleading, it's ppm tolerance.

    // --- Glycan Databases and Definitions ---
    /** Path to the glycan residue definition file. Optional, uses internal default if not provided. (--residuedb) */
    public Path residueDBpath;
    /** Path to the glycan modification definition file. Optional, uses internal default if not provided. (--glycanmoddb) */
    public Path glycanModPath;
    /** Internal path to the default glycan residue definitions resource file. */
    private static final String defaultResiduePath = "glycan_residues.txt";
    /** Internal path to the default glycan modification definitions resource file. */
    private static final String defaultModPath = "glycan_mods.txt";
    /** Map storing parsed glycan residues and modifications, keyed by their string representation. */
    public HashMap<String, GlycanResidue> glycanResidues;
    /** Map storing parsed glycan modifications, keyed by their string representation. */
    public HashMap<String, GlycanMod> glycanMods;
    /** List of specific residues/glycans to generate "dummy" glycopeptides for, based on --toaddresiduals. */
    public ArrayList<Glycan> toAddGlycans;
    /** Raw strings provided via --toaddresiduals before parsing. */
    private String[] residualStrings;

    // --- IonQuant and Processing Parameters ---
    /** Flag to indicate if the data is NOT from PASEF (TIMS-TOF). Affects IonQuant parameters? Default: false. (--nopasef) */
    public boolean noPASEF = false; // TODO: Verify exact effect of this flag
    /** Flag to indicate if TMT quantification data is present/should be considered. Default: false. (--runtmt) */
    public boolean runTMT = false; // TODO: Verify exact effect of this flag
    /** Maximum difference (rounds of inference) between parent and inferred glycans. Continue rounds of inference on found peaks up to this many times */
    public int expandDB = 0;
    /** Number of missing glycoforms to allow while expanding database. I.e., skips=1 allows one round of continued inference on a missed peak */
    public int maxSkips = 0;
    /** Minimum isotopes required for a feature in IonQuant. Default: 2. (--minisotopes) */
    public int minIsotopesIonQuant = 2;
    /** Minimum scans required for a feature in IonQuant. Default: 3. (--minscans) */
    public int minScansIonQuant = 3;
    /** Number of threads to use for parallel processing steps. Default: 4. (--numthreads) */
    public int numThreads = 4;
    /** Flag to indicate if allowing chimeric spectra when MBG trying to find MS2 hits */
    public boolean allowChimeric = true;
    // --- Internal State ---
    /** ArrayList storing the matched scan numbers, if allowChimeric is false, then MBP will not find MS2 in this list. */
    public HashMap<String, HashSet<Integer>> excludeScanNumbersMAp = new HashMap<>();
    /** Map storing the mapping from base raw file name (e.g., "file1") to its full path from the manifest. */
    public HashMap<String, String> msFilePathsMap;
    /** Map storing residues specified by --toaddresiduals. Currently seems unused after refactoring to toAddGlycans. */
    public HashMap<String, GlycanResidue> toAddResidues; // TODO: Remove if truly unused?


    /**
     * Parses command line arguments to configure the MBG analysis parameters.
     * Handles file path validation and loading of glycan definitions.
     *
     * @param args Command line arguments passed to the application.
     */
    public Params(String[] args) {
        // Iterate through command line arguments, assuming pairs of flag and value
        for (int i = 1; i < args.length; i += 2) {
            String flag = args[i];
            String value = args[i + 1].trim(); // Get value and remove leading/trailing whitespace

            switch (flag) {
                case "--psm": // Path to the input PSM file
                    psmPath = checkFilePath(value, "psm");
                    break;
                case "--manifest": // Path to the FragPipe manifest file
                    manifestPath = checkFilePath(value, "manifest");
                    msFilePathsMap = readManifest(); // Read the manifest immediately
                    break; // Added missing break statement
                case "--residuedb": // Path to the glycan residue definitions
                    residueDBpath = checkFilePath(value, "residue database");
                    break;
                case "--glycanmoddb": // Path to the glycan modification definitions
                    glycanModPath = checkFilePath(value, "glycan mods database");
                    break;
                case "--maxq": // Maximum glycan q-value threshold
                    try {
                        maxGlycanQval = Double.parseDouble(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing max glycan q-value argument '%s': not a valid number%n", value);
                        System.exit(1);
                    }
                    break;
                case "--minpsms": // Minimum PSMs per parent glycopeptide (glycan + peptide combination)
                    try {
                        minPSMs = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing min PSMs argument '%s': not a valid integer%n", value);
                        System.exit(1);
                    }
                    break;
                case "--minglycans": // Minimum different glycans per peptide
                    try {
                        minGlycans = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing min Glycans argument '%s': not a valid integer%n", value);
                        System.exit(1);
                    }
                    break;
                case "--rttol": // Retention time tolerance
                    try {
                        rtTol = Float.parseFloat(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing RT tolerance argument '%s': not a valid number%n", value);
                        System.exit(1);
                    }
                    break;
                case "--imtol": // Ion mobility tolerance
                    try {
                        imTol = Float.parseFloat(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing IM tolerance argument '%s': not a valid number%n", value);
                        System.exit(1);
                    }
                    break;
                case "--fdr": // False Discovery Rate threshold
                    try {
                        fdr = Float.parseFloat(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing FDR argument '%s': not a valid number%n", value);
                        System.exit(1);
                    }
                    break;
                case "--mztol": // Mass tolerance (ppm)
                     try {
                        mztol = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing m/z tolerance (ppm) argument '%s': not a valid integer%n", value);
                        System.exit(1);
                    }
                    break; // Added missing break statement
                case "--nopasef": // Flag for non-PASEF data
                    try {
                        noPASEF = Boolean.parseBoolean(value);
                    } catch (Exception e) { // Catch potential issues beyond just format
                        out.printf("Error parsing --nopasef argument '%s': must be 'true' or 'false'%n", value);
                        System.exit(1);
                    }
                    break;
                case "--runtmt": // Flag for TMT experiment
                    try {
                        runTMT = Boolean.parseBoolean(value);
                    } catch (Exception e) {
                        out.printf("Error parsing --runtmt argument '%s': must be 'true' or 'false'%n", value);
                        System.exit(1);
                    }
                    break;
                case "--expanddb": // Database expansion level
                    try {
                        expandDB = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) { // More specific catch
                        out.printf("Error parsing --expanddb argument '%s': not a valid integer%n", value);
                        System.exit(1);
                    }
                    break;
                case "--maxskips": // Database expansion level
                    try {
                        maxSkips = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) { // More specific catch
                        out.printf("Error parsing --maxskips argument '%s': not a valid integer%n", value);
                        System.exit(1);
                    }
                    break;
                case "--minscans": // Minimum scans for IonQuant
                    try {
                        minScansIonQuant = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing --minscans argument '%s': must be an integer%n", value);
                        System.exit(1);
                    }
                    break;
                case "--minisotopes": // Minimum isotopes for IonQuant
                    try {
                        minIsotopesIonQuant = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing --minisotopes argument '%s': must be an integer%n", value);
                        System.exit(1);
                    }
                    break;
                case "--numthreads": // Number of processing threads
                    try {
                        numThreads = Integer.parseInt(value);
                    } catch (NumberFormatException numberFormatException) {
                        out.printf("Error parsing --numthreads argument '%s': must be an integer%n", value);
                        System.exit(1);
                    }
                    break;
                case "--toaddresiduals": // Comma-separated list of residues/glycans to add
                    try {
                        residualStrings = value.split(",");
                    } catch (Exception e) { // Catch potential issues during split
                        out.printf("Error parsing --toaddresiduals argument '%s'%n", value);
                        e.printStackTrace(); // Print stack trace for more details
                        System.exit(1);
                    }
                    break;
                case "--allowchimeric": // Comma-separated list of residues/glycans to add
                    try {
                        allowChimeric = Boolean.parseBoolean(value);
                    } catch (Exception e) { // Catch potential issues during split
                        out.printf("Error parsing --allowchimeric argument '%s': must be 'true' or 'false'%n", value);
                        System.exit(1);
                    }
                    break;
                default:
                    out.printf("Warning: Unknown command line argument '%s'. Ignoring.%n", flag);
                    i--; // Decrement i because this flag didn't have a matching value, so the next arg is potentially a flag
                    break;
            }
        }

        // --- Post-Parsing Validation and Setup ---

        // Ensure required PSM path is provided
        if (psmPath == null) {
            out.println("Error: PSM file path (--psm) not specified!");
            System.exit(1);
        }

        // Load glycan residue definitions (default or user-specified)
        out.println("Loading glycan residue definitions...");
        if (residueDBpath == null) {
            out.printf("No --residuedb specified, using internal default: %s%n", defaultResiduePath);
            glycanResidues = GlycanParser.parseGlycoResiduesDBStream(Params.class.getResourceAsStream(defaultResiduePath));
        } else {
            out.printf("Using residue definitions from: %s%n", residueDBpath);
            glycanResidues = GlycanParser.parseGlycoResiduesDB(residueDBpath.toString());
        }
        if (glycanResidues == null || glycanResidues.isEmpty()) {
             out.println("Error: Failed to load glycan residue definitions.");
             System.exit(1);
        }
         out.printf("Loaded %d residue definitions.%n", glycanResidues.size());


        // Load glycan modification definitions (default or user-specified)
        out.println("Loading glycan modification definitions...");
        if (glycanModPath == null) {
             out.printf("No --glycanmoddb specified, using internal default: %s%n", defaultModPath);
            glycanMods = GlycanParser.parseGlycoModsDBStream(Params.class.getResourceAsStream(defaultModPath), glycanResidues.size(), glycanResidues);
        } else {
            out.printf("Using modification definitions from: %s%n", glycanModPath);
            glycanMods = GlycanParser.parseGlycoModsDB(glycanModPath.toString(), glycanResidues.size(), glycanResidues);
        }
         if (glycanMods == null) { // Can be empty if no mods are defined
             out.println("Warning: Failed to load glycan modification definitions or none were defined.");
             glycanMods = new HashMap<>(); // Ensure it's not null
         } else {
            out.printf("Loaded %d modification definitions.%n", glycanMods.size());
         }


        // Merge modifications into the main residue map for unified lookup
        for (Map.Entry<String, GlycanMod> modEntry : glycanMods.entrySet()) {
            glycanResidues.put(modEntry.getKey(), modEntry.getValue());
        }
         out.printf("Total glycan residues/mods available: %d%n", glycanResidues.size());

        // Parse the residues/glycans specified for addition (used for generating dummy glycopeptides)
        if (residualStrings != null) {
            out.println("Parsing residues/glycans specified by --toaddresiduals...");
            toAddGlycans = getToAddGlycans(residualStrings, glycanResidues);
             out.printf("Parsed %d entries for glycan expansion.%n", toAddGlycans.size());
        } else {
             out.println("No --toaddresiduals specified, skipping specific glycan expansion list.");
            toAddGlycans = new ArrayList<>(); // Initialize to empty list
        }

        // TODO: Decide if this map is still needed or should be removed.
        // toAddResidues = getToAddResidues(residualStrings, glycanResidues);
    }

    /**
     * Checks if a given file path exists and is accessible. Exits application if not.
     * Removes surrounding quotes from the input path string.
     *
     * @param inputPath The file path string provided by the user.
     * @param pathName  A descriptive name for the file type (e.g., "psm", "manifest") used in error messages.
     * @return A validated Path object.
     */
    public static Path checkFilePath(String inputPath, String pathName) {
        Path path = null;
        String cleanedPath = inputPath.replaceAll("^['\"]|['\"]$", ""); // Remove leading/trailing quotes
        try {
            path = Paths.get(cleanedPath);
        } catch (Exception e) {
            e.printStackTrace();
            out.printf("Error parsing %s file path: [%s] (Original: '%s')%n", pathName, cleanedPath, inputPath);
            System.exit(1);
        }
        if (!Files.exists(path)) {
            out.printf("Error: %s file path does not exist: [%s]%n", pathName, path.toAbsolutePath());
            System.exit(1);
        }
         if (!Files.isReadable(path)) {
            out.printf("Error: %s file is not readable: [%s]%n", pathName, path.toAbsolutePath());
            System.exit(1);
        }
        return path;
    }

    /**
     * Creates a map of residue names to GlycanResidue objects based on input strings.
     * Note: This method seems superseded by getToAddGlycans and might be unused.
     *
     * @param residueNames Array of residue name strings.
     * @param allResidues  Map of all known residue/mod names to their GlycanResidue objects.
     * @return A map containing only the specified residues.
     */
    private HashMap<String, GlycanResidue> getToAddResidues(String[] residueNames, HashMap<String, GlycanResidue> allResidues) {
        // TODO: Verify if this method is still needed.
        if (residueNames == null) return new HashMap<>(); // Handle null input

        HashMap<String, GlycanResidue> selectedResidues = new HashMap<>();
        for (String name : residueNames) {
            GlycanResidue residue = allResidues.get(name);
            if (residue != null) {
                selectedResidues.put(name, residue);
            } else {
                out.printf("Warning: Residue '%s' specified in --toaddresiduals not found in definitions.%n", name);
            }
        }
        return selectedResidues;
    }

    /**
     * Parses an array of strings, converting each into a Glycan object.
     * Strings containing parentheses are treated as complex glycan compositions (e.g., "HexNAc(2)Hex(3)"),
     * while others are treated as single residues (e.g., "HexNAc").
     *
     * @param compositionStrings Array of strings representing residues or glycan compositions.
     * @param allResidues        Map of all known residue/mod names to their GlycanResidue objects for lookup.
     * @return An ArrayList of parsed Glycan objects.
     */
    private ArrayList<Glycan> getToAddGlycans(String[] compositionStrings, HashMap<String, GlycanResidue> allResidues) {
        if (compositionStrings == null) return new ArrayList<>(); // Handle null input

        ArrayList<Glycan> parsedGlycans = new ArrayList<>();
        for (String compString : compositionStrings) {
            compString = compString.trim(); // Trim whitespace
            if (compString.isEmpty()) continue; // Skip empty strings

            try {
                if (compString.contains("(")) {
                    // Treat as a complex glycan string (e.g., "HexNAc(2)Hex(3)")
                    Glycan parsed = GlycanParser.parseGlycanString(compString, allResidues);
                    if (parsed != null) {
                        parsedGlycans.add(parsed);
                    } else {
                         out.printf("Warning: Could not parse glycan string '%s' from --toaddresiduals.%n", compString);
                    }
                } else {
                    // Treat as a single residue name (e.g., "HexNAc")
                    GlycanResidue residue = allResidues.get(compString);
                    if (residue != null) {
                        TreeMap<GlycanResidue, Integer> compositionMap = new TreeMap<>();
                        compositionMap.put(residue, 1);
                        parsedGlycans.add(new Glycan(compositionMap));
                    } else {
                        out.printf("Warning: Residue '%s' from --toaddresiduals not found in definitions.%n", compString);
                    }
                }
            } catch (Exception e) {
                 out.printf("Error parsing entry '%s' from --toaddresiduals: %s%n", compString, e.getMessage());
                 // Optionally continue or exit depending on desired strictness
            }
        }
        return parsedGlycans;
    }

    /**
     * Reads a FragPipe manifest file (.fp-manifest) to create a map between
     * base raw file names (without extension) and their full file paths.
     * Handles potential path separator issues (Windows '\' vs universal '/').
     *
     * @return A HashMap where keys are base raw file names and values are full paths.
     */
    private HashMap<String, String> readManifest() {
        HashMap<String, String> msFilesMap = new HashMap<>();
        if (manifestPath == null) {
            out.println("Warning: Manifest file path not provided (--manifest). Cannot map raw file names.");
            return msFilesMap; // Return empty map
        }

        out.printf("Reading manifest file: %s%n", manifestPath);
        try (BufferedReader reader = new BufferedReader(new FileReader(manifestPath.toFile()))) {
            String line;
            int lineNum = 0;
            while ((line = reader.readLine()) != null) {
                lineNum++;
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) { // Skip empty lines and comments
                    continue;
                }

                String[] columns = line.split("\t"); // Manifest is typically tab-separated
                if (columns.length < 1) {
                     out.printf("Warning: Skipping malformed line %d in manifest: '%s'%n", lineNum, line);
                    continue;
                }

                String fullPathString = columns[0].trim();
                try {
                    Path fullPath = Paths.get(fullPathString);
                    String baseFileNameWithExt = fullPath.getFileName().toString();
                    int dotIndex = baseFileNameWithExt.lastIndexOf('.');
                    String baseFileName = (dotIndex == -1) ? baseFileNameWithExt : baseFileNameWithExt.substring(0, dotIndex);

                    // Store the original path string, potentially replacing backslashes for consistency if needed later,
                    // but often better to keep the original format from the manifest.
                    // String storedPath = fullPathString.replace("\\", "/"); // Example: Convert to forward slashes
                    String storedPath = fullPathString; // Keep original path from manifest

                    if (msFilesMap.containsKey(baseFileName)) {
                         out.printf("Warning: Duplicate base file name '%s' found in manifest (Line %d). Overwriting previous entry.%n", baseFileName, lineNum);
                    }
                    msFilesMap.put(baseFileName, storedPath);

                } catch (Exception e) {
                     out.printf("Warning: Could not process path '%s' from manifest line %d: %s%n", fullPathString, lineNum, e.getMessage());
                }
            }
        } catch (IOException ioException) {
            out.printf("Error reading FragPipe manifest file '%s': %s%n", manifestPath, ioException.getMessage());
            ioException.printStackTrace();
            // Decide whether to exit or continue without the map
            // System.exit(1);
             out.println("Proceeding without manifest file mapping.");
             return new HashMap<>(); // Return empty map on error
        }
         out.printf("Read %d file mappings from manifest.%n", msFilesMap.size());
        return msFilesMap;
    }
}
