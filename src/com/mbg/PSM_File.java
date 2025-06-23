package com.mbg;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.*;

import static java.nio.file.Files.createTempFile;

/**
 * Represents and manages a Peptide Spectrum Match (PSM) file (typically .psm.tsv).
 * Responsibilities include:
 * - Reading the PSM file header and content.
 * - Parsing each line into a PSM object.
 * - Mapping column headers to standardized PSMattr enums.
 * - Adding necessary MBG-specific columns if they don't exist.
 * - Grouping PSMs by their source raw file name.
 * - Writing the original and newly inferred PSMs back to the file system.
 * - Providing utility methods for generating mass lists required by external tools (IonQuant).
 */
public class PSM_File {

    /** The header line split into column names. */
    public String[] headers;
    /** A map grouping PSMs by their source raw file name (key: raw file name without extension, value: list of PSMs). */
    public HashMap<String, ArrayList<PSM>> rawfilePSMmap;
    /** A list containing all PSM objects read from the file. */
    public final ArrayList<PSM> allPSMs;
    /** The File object representing the PSM file being processed. */
    public File file;
    /** A map linking standardized PSM attributes (PSMattr enum) to their corresponding column index in the headers array. */
    public HashMap<PSMattr, Integer> PSMcolumns;
    /** Static map defining the expected column header string for each PSMattr enum. Populated in the static initializer block. */
    private static final HashMap<PSMattr, String> columnNames;
    /** The number of new columns added to the header/PSM data (0 if columns already existed). */
    public int addCols;
    /** The index at which new MBG-specific columns should be inserted into the header and data rows. Determined during initialization. */
    public static int insertIndex; // Static might be problematic if multiple PSM_Files are handled concurrently. TODO: Review if static is appropriate here.
    // holds the maximum rank of spectra in the PSM file, used to determine the rank for adding MBG PSMs
    public final HashMap<Integer, Integer> maxRankOfSpectra;

    // Static initializer block to populate the map of standard column names.
    // This map defines the expected header string for each PSM attribute (enum).
    static {
        columnNames = new HashMap<>();
        // Mapping PSMattr enums to their expected header names in the PSM file
        columnNames.put(PSMattr.RT, "Retention"); // Retention time
        columnNames.put(PSMattr.peptide, "Peptide"); // Peptide sequence
        columnNames.put(PSMattr.protein, "Protein"); // Protein accession
        columnNames.put(PSMattr.IM, "Ion Mobility"); // Ion mobility value
        columnNames.put(PSMattr.protein_start, "Protein Start"); // Start position in protein
        columnNames.put(PSMattr.assigned_mods, "Assigned Modifications"); // Modifications from search engine
        columnNames.put(PSMattr.glycan, "Total Glycan Composition"); // Glycan composition from PTM-Shepherd/search
        columnNames.put(PSMattr.spectrum, "Spectrum"); // Spectrum identifier
        columnNames.put(PSMattr.delta_mass, "Delta Mass"); // Mass difference (Obs - Theo)
        columnNames.put(PSMattr.mod_peptide, "Modified Peptide"); // Peptide sequence with mods
        columnNames.put(PSMattr.CALIBRATED_OBSERVED_MASS, "Calibrated Observed Mass"); // Observed precursor mass
        columnNames.put(PSMattr.CALIBRATED_OBSERVED_MZ, "Calibrated Observed M/Z"); // Observed precursor m/z
        columnNames.put(PSMattr.Calculated_MASS, "Calculated Peptide Mass"); // Observed precursor mass
        columnNames.put(PSMattr.Calculated_MZ, "Calculated M/Z"); // Observed precursor m/z
        columnNames.put(PSMattr.Observed_MASS, "Observed Mass"); // Observed precursor mass
        columnNames.put(PSMattr.Observed_Mz, "Observed M/Z"); // Observed precursor m/z
        columnNames.put(PSMattr.charge, "Charge"); // Precursor charge
        // MBG specific columns start here ("Spectrum File" and "Intensity" will be modified in MBG PSMs)
        columnNames.put(PSMattr.has_ms2_spec, "Has MS2 Spectrum?"); // Did the original PSM have an MS2?
        columnNames.put(PSMattr.has_Y1, "Has Y1 ion?"); // Was Y1 ion observed in MS2?
        columnNames.put(PSMattr.MBG_score, "MBG Score"); // Score from MBG LDA model
        columnNames.put(PSMattr.mass_error, "Mass error(ppm)"); // Mass error from IonQuant matching
        columnNames.put(PSMattr.KL, "KL divergence"); // KL divergence score from IonQuant
        columnNames.put(PSMattr.shift_frequency, "Shift Frequency"); // Frequency of the observed glycan shift
        columnNames.put(PSMattr.rt_diff, "RT Difference"); // Predicted vs. Observed RT difference
        columnNames.put(PSMattr.im_diff, "IM Difference"); // Predicted vs. Observed IM difference
        columnNames.put(PSMattr.if_new_glycan, "MBG only composition"); // Was this glycan composition novel?
        columnNames.put(PSMattr.confidence_level, "Confidence Level"); // Confidence category based on the existence of MS2 hit and Y1 ion
        columnNames.put(PSMattr.parent_glycan, "Parent Glycan"); // Glycan used as basis for inference
        columnNames.put(PSMattr.added_residual, "Difference from parent"); // Residue(s) added/removed for inference
        columnNames.put(PSMattr.spectrum_file, "Spectrum File"); // Source raw file name
        columnNames.put(PSMattr.intensity, "Intensity"); // Precursor intensity
        // columns that need to leave blank for MBG PSMs
        columnNames.put(PSMattr.BestPositions, "Best Positions");
        columnNames.put(PSMattr.expectation, "expectation");
        columnNames.put(PSMattr.glyco_q_val, "Glycan q-value"); // Glycan assignment confidence
        columnNames.put(PSMattr.Glycanscore, "Glycan score");
        columnNames.put(PSMattr.hyperscore, "hyperscore");
        columnNames.put(PSMattr.IonsAllUnshifted, "Ions All Unshifted");
        columnNames.put(PSMattr.IonsBestPosition, "Ions Best Position");
        columnNames.put(PSMattr.IonsSecondBestPosition, "Ions Second Best Position");
        columnNames.put(PSMattr.IonsShiftedAllPositions, "Ions Shifted All Positions");
        columnNames.put(PSMattr.MSFraggerLocalization, "MSFragger Localization");
        columnNames.put(PSMattr.IonsShiftedBestPosition, "Ions Shifted Best Position");
        columnNames.put(PSMattr.nextscore, "nextscore");
        columnNames.put(PSMattr.NumberBestPositions, "Number Best Positions");
        columnNames.put(PSMattr.PositionScores, "Position Scores");
        columnNames.put(PSMattr.probability, "probability");
        columnNames.put(PSMattr.ScoreAllUnshifted, "Score All Unshifted");
        columnNames.put(PSMattr.ScoreBestPosition, "Score Best Position");
        columnNames.put(PSMattr.ScoreSecondBestPosition, "Score Second Best Position");
        columnNames.put(PSMattr.ScoreShiftedAllPositions, "Score Shifted All Positions");
        columnNames.put(PSMattr.ScoreShiftedBestPosition, "Score Shifted Best Position");
        columnNames.put(PSMattr.ShiftedOnlyPositionIons, "Shifted Only Position Ions");
        columnNames.put(PSMattr.ShiftedOnlyPositionScores, "Shifted Only Position Scores");
    }
    // PSM table headers to add to the PSM file if they do not already exist.
    private static final String[] headersToAdd = new String[] {
            // Insert in reverse order of definition to maintain desired final order
            columnNames.get(PSMattr.added_residual),
            columnNames.get(PSMattr.parent_glycan),
            columnNames.get(PSMattr.confidence_level),
            columnNames.get(PSMattr.if_new_glycan),
            columnNames.get(PSMattr.im_diff),
            columnNames.get(PSMattr.rt_diff),
            columnNames.get(PSMattr.shift_frequency),
            columnNames.get(PSMattr.KL),
            columnNames.get(PSMattr.mass_error),
            columnNames.get(PSMattr.MBG_score),
            columnNames.get(PSMattr.has_Y1),
            columnNames.get(PSMattr.has_ms2_spec)
    };

    /**
     * Constructor for PSM_File. Reads the PSM file, initializes column mappings,
     * adds MBG-specific columns if needed, and parses all PSM entries.
     *
     * @param params The application parameters containing the PSM file path and other settings.
     * @throws Exception If there's an error reading the file or parsing data.
     */
    public PSM_File(Params params) throws Exception {
        BufferedReader psmReader = new BufferedReader(new FileReader(params.psmPath.toFile()), 1 << 22); // Use a large buffer
        this.file = params.psmPath.toFile(); // Store the File object
        this.headers = psmReader.readLine().split("\t"); // Read and split the header line
        this.PSMcolumns = initPSMcolumns(); // Map header names to column indices
        allPSMs = new ArrayList<>(); // Initialize the list to store all PSM objects
        maxRankOfSpectra = new HashMap<>();     // store max rank of all scans in the PSM file to know what rank to use when adding MBG PSMs

        // Check if MBG-specific columns already exist. If not, add them to the header.
        if (PSMcolumns.getOrDefault(PSMattr.MBG_score, -1) == -1) {
            System.out.println("MBG-specific columns not found in PSM file. Adding them...");
            // Determine the insertion point for new columns (typically after glycan q-value or glycan composition)
            if (PSMcolumns.containsKey(PSMattr.glyco_q_val)) {
                insertIndex = PSMcolumns.get(PSMattr.glyco_q_val) + 1;
            } else if (PSMcolumns.containsKey(PSMattr.glycan)) {
                insertIndex = PSMcolumns.get(PSMattr.glycan) + 1;
            } else {
                // Fallback: insert after peptide sequence if glycan columns are missing (less ideal)
                insertIndex = PSMcolumns.getOrDefault(PSMattr.peptide, 0) + 1;
                System.out.println("Warning: Glycan q-value and Glycan Composition columns not found. Inserting MBG columns after Peptide.");
            }

            // Create a new header list and insert the MBG column names at the determined index
            ArrayList<String> newHeadersList = new ArrayList<>(Arrays.asList(this.headers));
            addCols = headersToAdd.length; // Number of MBG columns to add
            for (String header : headersToAdd) {
                newHeadersList.add(insertIndex, header);
            }

            headers = newHeadersList.toArray(new String[0]); // Update the main headers array
            this.PSMcolumns = initPSMcolumns(); // Re-initialize column mappings with the new headers
            System.out.printf("Added %d MBG columns starting at index %d.%n", addCols, insertIndex);
        } else {
            addCols = 0; // No columns were added
            System.out.println("Found existing MBG-specific columns in PSM file.");
        }

        // Read PSM data lines, create PSM objects, and group them by raw file
        rawfilePSMmap = new HashMap<>();
        String line;
        int lineNum = 1; // Start from 1 for header
        while((line = psmReader.readLine()) != null) {
            lineNum++;
            if (!line.trim().isEmpty()) { // Skip empty lines
                try {
                    // Create a new PSM object, passing necessary context
                    // The boolean flags indicate this is an original (not inferred) PSM
                    PSM psm = new PSM(line, PSMcolumns, params, addCols);
                    allPSMs.add(psm); // Add to the master list
                    maxRankOfSpectra.put(psm.scanNum, Math.max(maxRankOfSpectra.getOrDefault(psm.scanNum, 1), psm.rank)); // Update max rank for this spectrum

                    // Add the PSM to the map grouped by raw file name
                    if (psm.rawfile != null && !psm.rawfile.isEmpty()) {
                        // computeIfAbsent gets the list if key exists, or creates a new ArrayList if not, then adds the psm
                        rawfilePSMmap.computeIfAbsent(psm.rawfile, k -> new ArrayList<>()).add(psm);
                    } else {
                         System.out.printf("Warning: PSM at line %d has missing or empty Spectrum File information. Cannot group by raw file.%n", lineNum);
                    }
                } catch (Exception e) {
                     System.out.printf("Error parsing PSM at line %d: %s%n", lineNum, e.getMessage());
                     e.printStackTrace(); // Log the error but continue processing other lines if possible
                }
            }
        }
        psmReader.close(); // Close the reader
        // update the max rank of spectra in each PSM to know what rank to use when adding MBG PSMs
        for (PSM psm: allPSMs) {
            psm.maxRank = maxRankOfSpectra.getOrDefault(psm.scanNum, 1); // Default to 1 if not found
        }
        System.out.printf("Read %d total PSMs from %d raw files.%n", allPSMs.size(), rawfilePSMmap.size());
    }

    /**
     * Initializes the PSMcolumns map by finding the index of each standard column
     * defined in columnNames within the actual file headers.
     * Issues warnings for missing standard columns (excluding MBG-specific ones initially).
     *
     * @return A HashMap mapping PSMattr enums to their column index. Returns an empty map if headers are null.
     */
    private HashMap<PSMattr, Integer> initPSMcolumns() {
        HashMap<PSMattr, Integer> psmColumnsMap = new HashMap<>();
        if (this.headers == null) {
            System.err.println("Error: Headers array is null in initPSMcolumns. Cannot map columns.");
            return psmColumnsMap; // Return empty map
        }

        // Iterate through the standard column definitions (PSMattr enum mapped to expected header string)
        for (Map.Entry<PSMattr, String> nameEntry: columnNames.entrySet()) {
            String expectedHeader = nameEntry.getValue();
            PSMattr attribute = nameEntry.getKey();
            int index = getColumnIndex(expectedHeader); // Find the index in the current file headers

            if (index != -1) {
                psmColumnsMap.put(attribute, index); // Store the mapping if the header is found
            } else {
                // Only issue a warning if a non-MBG standard column is missing.
                // MBG columns might be added later, so don't warn about them yet.
                if (!isMbgSpecificColumn(attribute)) {
                    // Example: Warn if "Peptide" or "Retention" is missing.
                    // System.out.printf("Warning: Expected standard column '%s' not found in PSM file header.%n", expectedHeader);
                    // Depending on strictness, could potentially throw an error for critical missing columns.
                }
                // If an MBG column is missing here, it's expected if 'addCols' > 0.
                // If 'addCols' == 0, it means we expected MBG columns but didn't find one, which might indicate an issue.
                // However, the check in the constructor handles the primary logic of adding them if needed.
            }
        }
        // TODO: Add error checking here? E.g., ensure critical columns like 'Peptide', 'Spectrum File' are present.
        if (!psmColumnsMap.containsKey(PSMattr.peptide)) {
             System.err.println("CRITICAL WARNING: 'Peptide' column not found in PSM file. Processing may fail.");
        }
         if (!psmColumnsMap.containsKey(PSMattr.spectrum_file)) {
             System.err.println("CRITICAL WARNING: 'Spectrum File' column not found in PSM file. Grouping by raw file will fail.");
        }
        return psmColumnsMap;
    }

     /**
     * Helper method to check if a PSMattr corresponds to an MBG-specific output column.
     * Used to suppress warnings during initial column mapping before MBG columns might be added.
     * @param attr The PSMattr to check.
     * @return true if it's an MBG-specific column, false otherwise.
     */
    private boolean isMbgSpecificColumn(PSMattr attr) {
        // Check against all MBG-specific enum values
        return attr == PSMattr.has_ms2_spec ||
               attr == PSMattr.has_Y1 || attr == PSMattr.MBG_score ||
               attr == PSMattr.mass_error || attr == PSMattr.KL ||
               attr == PSMattr.shift_frequency || attr == PSMattr.rt_diff ||
               attr == PSMattr.im_diff || attr == PSMattr.if_new_glycan ||
               attr == PSMattr.confidence_level || attr == PSMattr.parent_glycan ||
               attr == PSMattr.added_residual;
    }


    /**
     * Finds the index of a specific column name within the headers array.
     * Performs a case-insensitive comparison.
     *
     * @param headerName The name of the column header to find.
     * @return The 0-based index of the column, or -1 if not found.
     */
    public int getColumnIndex(String headerName) {
        if (headers == null || headerName == null) {
            return -1;
        }
        for(int i = 0; i < headers.length; i++) {
            // Use equalsIgnoreCase for robustness against case variations in header names
            if(headerName.equalsIgnoreCase(headers[i])) {
                return i;
            }
        }
        return -1; // Header not found
    }

    /**
     * Writes the original PSMs and the newly inferred PSMs to a temporary file,
     * sorts them (if possible), and then replaces the original PSM file with the updated content.
     * Includes a basic backup mechanism (.mbg.bak) and uses safe file replacement via move/rename.
     *
     * @param inferredPSMs An ArrayList containing the newly inferred PSM objects generated by MBG.
     * @throws IOException If there is an error writing to the file system during the process.
     */
    public void printNewPSMs(ArrayList<PSM> inferredPSMs) throws IOException {
        // Define temporary and backup file names based on the original file path
        File originalFile = this.file;
        // Use a more specific temporary file extension
        File tempFile = new File(originalFile.getAbsolutePath() + ".mbg.tmp");
        File backupFile = new File(originalFile.getAbsolutePath() + ".mbg.bak");

        System.out.printf("Writing updated PSM data to temporary file: %s%n", tempFile.getName());
        try (PrintWriter psmWriter = new PrintWriter(new FileWriter(tempFile))) {

            // Write the potentially updated header line
            psmWriter.println(String.join("\t", this.headers));

            // Combine the original PSMs with the newly inferred ones
            ArrayList<PSM> combinedPSMs = new ArrayList<>(allPSMs); // Start with original PSMs
            combinedPSMs.addAll(inferredPSMs); // Add the inferred PSMs
            System.out.printf("Combined %d original and %d inferred PSMs for output.%n", allPSMs.size(), inferredPSMs.size());

            // Sort the combined list. Assumes PSM implements Comparable (e.g., sorting by spectrum name or RT).
            // Sorting helps maintain a consistent order in the output file.
            try {
                 Collections.sort(combinedPSMs);
                 System.out.println("Sorted combined PSMs.");
            } catch (Exception e) {
                 // Log a warning if sorting fails, but continue writing in potentially unsorted order.
                 System.out.println("Warning: Could not sort combined PSMs. Output order might be inconsistent. Error: " + e.getMessage());
            }

            // Write each PSM object's data line (which is stored in its 'splits' array)
            for (PSM psm : combinedPSMs) {
                // Ensure psm.splits is not null before joining
                if (psm != null && psm.splits != null) {
                    psmWriter.println(String.join("\t", psm.splits));
                } else {
                    System.out.println("Warning: Skipping null PSM or PSM with null splits during output writing.");
                }
            }
        } // PrintWriter is automatically closed here (try-with-resources)

        System.out.println("Finished writing to temporary file.");

        // --- Safely replace the original file with the temporary file ---
        try {
            // 1. Create a backup of the original file (if it exists)
            if (originalFile.exists()) {
                // Use Files.move for atomic operation where possible, replacing existing backup
                Files.move(originalFile.toPath(), backupFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                System.out.printf("Created backup of original file: %s%n", backupFile.getName());
            } else {
                 System.out.println("Original file does not exist, no backup created.");
            }

            // 2. Move the temporary file to the original file name
            // This is generally safer than delete + rename, especially on some OSes.
            Files.move(tempFile.toPath(), originalFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
            System.out.printf("Successfully replaced original PSM file with updated version: %s%n", originalFile.getName());

            // 3. Optional: Delete the backup file if the replacement was successful
            // Commented out by default for safety - allows manual recovery from backup if needed.
            // Files.deleteIfExists(backupFile.toPath());
            // System.out.printf("Deleted backup file: %s%n", backupFile.getName());

        } catch (IOException e) {
            // Log detailed error message if replacement fails
            System.err.printf("Error replacing original PSM file '%s' with temporary file '%s'.%n", originalFile.getName(), tempFile.getName());
            System.err.println("The updated data might be found in the temporary file: " + tempFile.getAbsolutePath());
            if (backupFile.exists()) {
                System.err.println("The original data might be found in the backup file: " + backupFile.getAbsolutePath());
            }
            e.printStackTrace();
            // Rethrow the exception so the caller knows the operation failed
            throw e;
        } finally {
             // 4. Ensure the temporary file is deleted if it still exists (e.g., if the move failed)
             try {
                 Files.deleteIfExists(tempFile.toPath());
                 Files.deleteIfExists(backupFile.toPath());
             } catch (IOException ignored) {
                 // Ignore errors during cleanup deletion
             }
        }
    }

    /**
     * Writes a list of calculated masses to the 'MBG_masses_list.txt' file.
     * This file is likely used as input for IonQuant or a similar tool for MS1 feature finding.
     * It also triggers appending the unique masses from this list to 'modmasses_ionquant.txt'.
     *
     * @param massesToWrite A List of Double values representing the masses to write.
     * @throws IOException If an error occurs during file writing.
     */
    public static void massWriter(List<Double> massesToWrite) throws IOException {
//        String mbgMassListPath = "MBG_masses_list.txt"; // File specifically for MBG masses
        String ionQuantMassPath = "modmasses_ionquant.txt"; // File potentially used by IonQuant
//
//        System.out.printf("Writing %d masses to %s%n", massesToWrite.size(), mbgMassListPath);
////         Using try-with-resources to ensure the writer is closed automatically
//        try (BufferedWriter massListWriter = new BufferedWriter(new FileWriter(mbgMassListPath))) {
//            // Write each mass to the MBG-specific file
//            for (Double mass : massesToWrite) {
//                if (mass != null) { // Avoid writing null values
//                    // Format mass to a fixed number of decimal places for consistency
//                    massListWriter.write(String.format(Locale.US, "%.6f", mass) + "\n");
//                }
//            }
//            System.out.printf("Finished writing masses to %s.%n", mbgMassListPath);
//        } catch (IOException e) {
//            System.err.printf("An error occurred while writing to %s: %s%n", mbgMassListPath, e.getMessage());
//            throw e; // Rethrow the exception so the caller is aware
//        }

        // Append the unique masses from the MBG list to the IonQuant list
//        System.out.printf("Appending unique masses from %s to %s%n", mbgMassListPath, ionQuantMassPath);
        try {
            appendUniqueMasses(ionQuantMassPath, massesToWrite);
            System.out.printf("Finished appending masses to %s.%n", ionQuantMassPath);
        } catch (IOException e) {
             System.err.printf("An error occurred while appending masses to %s: %s%n", ionQuantMassPath, e.getMessage());
             // Decide if this error is critical or can be ignored. Currently just prints error.
             // throw e; // Uncomment to make this error critical
        }
    }

    /**
     * Appends unique masses from a given list to a target file, ensuring no duplicate masses are added.
     * Existing masses in the target file are compared using a scaled integer representation (e.g., multiplied by 10,000)
     * to account for floating-point precision issues. This allows matching masses up to a fixed decimal precision.
     *
     * If the target file does not exist, a new file is created and all valid, unique masses are written to it.
     *
     * @param targetFilePath Path to the file to append masses to (e.g., "modmasses_ionquant.txt").
     * @param massesToWrite List of masses to be appended if not already present in the target file.
     * @throws IOException If an error occurs during file reading or writing.
     */

    private static void appendUniqueMasses(String targetFilePath, List<Double> massesToWrite) throws IOException {
        File targetFile = new File(targetFilePath);

        // --- Create target file if it does not exist ---
        if (!targetFile.exists()) {
            System.out.printf("Target file %s does not exist. Creating a new file.%n", targetFilePath);
            targetFile.createNewFile(); // Empty file created
        }

        // --- Prepare to detect duplicates in the target file ---
        Set<Long> existingMassesScaled = new HashSet<>();
        final long scaleFactor = 10000; // e.g., for 4 decimal places

        // 1. Read existing masses from the target file
        System.out.printf("Reading existing masses from target file: %s%n", targetFilePath);
        try (BufferedReader targetReader = new BufferedReader(new FileReader(targetFile))) {
            String line;
            while ((line = targetReader.readLine()) != null) {
                try {
                    double mass = Double.parseDouble(line.trim());
                    long scaledMass = Math.round(mass * scaleFactor);
                    existingMassesScaled.add(scaledMass);
                } catch (NumberFormatException e) {
                    // Skip invalid lines
                }
            }
        }
        System.out.printf("Found %d unique existing masses (scaled by %d) in target file.%n", existingMassesScaled.size(), scaleFactor);

        // 2. Append only new, unique masses from the list
        int appendedCount = 0;
        System.out.println("Appending unique masses from provided list...");
        try (BufferedWriter targetWriter = new BufferedWriter(new FileWriter(targetFile, true))) { // append mode
            for (Double mass : massesToWrite) {
                if (mass == null || Double.isNaN(mass) || Double.isInfinite(mass)) continue; // skip invalid values
                long scaledMass = Math.round(mass * scaleFactor);
                if (existingMassesScaled.add(scaledMass)) {
                    targetWriter.write(mass.toString());
                    targetWriter.newLine();
                    appendedCount++;
                }
            }
        }
        System.out.printf("Appended %d unique masses to %s.%n", appendedCount, targetFilePath);
    }



    /**
     * Enum defining standardized attributes (columns) expected or used within a PSM file.
     * This provides a type-safe way to refer to specific columns when parsing or writing.
     * Each enum constant maps to a standard column header name.
     */
    public enum PSMattr {
        // --- Core PSM Identification Info ---
        /** Peptide sequence without modifications. Header: "Peptide" */
        peptide,
        /** Peptide sequence with modification annotations (e.g., N(1)). Header: "Modified Peptide" */
        mod_peptide,
        /** Spectrum identifier (e.g., scan number, index, or title). Header: "Spectrum" */
        spectrum,
        /** Base name of the raw data file where the spectrum originated. Header: "Spectrum File" */
        spectrum_file,
        /** Precursor charge state. Header: "Charge" */
        charge,
        /** Retention time (typically in minutes). Header: "Retention" */
        RT,
        /** Ion mobility value (units can vary, e.g., 1/K0 or V.s.cm-2). Header: "Ion Mobility" */
        IM,
        /** Intensity of the precursor ion or MS1 feature. Header: "Intensity" */
        intensity,

        // --- Mass Information ---
        /** Observed precursor mass, potentially calibrated by the search engine. Header: "Calibrated Observed Mass" */
        CALIBRATED_OBSERVED_MASS,
        /** Observed precursor m/z, potentially calibrated. Header: "Calibrated Observed M/Z" */
        CALIBRATED_OBSERVED_MZ,
        Calculated_MASS,
        Calculated_MZ,
        Observed_MASS,
        Observed_Mz,
        /** Mass difference between observed and theoretical peptide mass (Da). Header: "Delta Mass" */
        delta_mass,

        // --- Protein Information ---
        /** Protein accession(s) the peptide maps to. Header: "Protein" */
        protein,
        /** Start position of the peptide sequence within the protein sequence(s). Header: "Protein Start" */
        protein_start,

        // --- Modification and Glycan Information ---
        /** Modifications identified by the primary search engine (excluding the main glycan). Header: "Assigned Modifications" */
        assigned_mods,
        /** Glycan composition string identified by PTM-Shepherd or glycan search. Header: "Total Glycan Composition" */
        glycan,
        /** Confidence score (e.g., q-value or probability) for the glycan assignment. Header: "Glycan q-value" */
        glyco_q_val,

        // --- MBG Inference Output Columns (Added by this tool if not present) ---
        /** Boolean flag indicating if the original identified PSM file contained an MS2 spectrum for this precursor. Header: "Has MS2 Spectrum?" */
        has_ms2_spec,
        /** Boolean flag indicating if a Y1 ion (peptide + core HexNAc) was observed in the MS2 spectrum (if available). Header: "Has Y1 ion?" */
        has_Y1,
        /** Score assigned by the MBG Linear Discriminant Analysis (LDA) model. Header: "MBG Score" */
        MBG_score,
        /** The glycan composition of the identified PSM that was used as the basis for inferring this glycoform. Header: "Parent Glycan" */
        parent_glycan,
        /** The glycan residue(s) or composition difference between the parent and the inferred glycan. Header: "Residual Difference" */
        added_residual,
        /** Mass error (ppm) observed by IonQuant when matching the inferred precursor m/z in MS1 data. Header: "Mass error(ppm)" */
        mass_error,
        /** KL divergence score from IonQuant's MS1 feature matching for the inferred glycoform. Header: "KL (IonQuant Score)" */
        KL,
        /** How often the specific glycan shift (parent -> inferred) was observed among confidently identified PSMs in the dataset. Header: "Shift Frequency" */
        shift_frequency,
        /** Difference between the predicted retention time (based on parent + shift) and the observed RT from IonQuant. Header: "RT Difference" */
        rt_diff,
        /** Difference between the predicted ion mobility (based on parent + shift) and the observed IM from IonQuant. Header: "IM Difference" */
        im_diff,
        /** Boolean flag indicating if the inferred glycan composition was not present in the set of originally identified glycans. Header: "New Glycan?" */
        if_new_glycan,
        /** Categorical confidence level (e.g., High, Medium) assigned based on the MBG score and FDR filtering. Header: "Confidence Level" */
        confidence_level,

        // columns to ignore
        hyperscore,
        nextscore,
        probability,
        expectation,
        MSFraggerLocalization,
        PositionScores,
        BestPositions,
        NumberBestPositions,
        ShiftedOnlyPositionScores,
        ShiftedOnlyPositionIons,
        ScoreBestPosition,
        IonsBestPosition,
        ScoreSecondBestPosition,
        IonsSecondBestPosition,
        ScoreAllUnshifted,
        IonsAllUnshifted,
        ScoreShiftedBestPosition,
        IonsShiftedBestPosition,
        ScoreShiftedAllPositions,
        IonsShiftedAllPositions,
        Glycanscore;

        // --- Potentially Deprecated/Unused Enums ---
        // ionMobilityDiff, // Seems redundant with im_diff? Check usage.

        // Note: The mapping to actual column names is done in the static 'columnNames' HashMap.
        // This enum provides the standardized internal representation.
    }

    public static final EnumSet<PSMattr> PSMattr_blank = EnumSet.of(
//            PSMattr.hyperscore,
//            PSMattr.nextscore,
//            PSMattr.probability,
//            PSMattr.expectation,
            PSMattr.MSFraggerLocalization,
            PSMattr.PositionScores,
            PSMattr.BestPositions,
            PSMattr.NumberBestPositions,
            PSMattr.ShiftedOnlyPositionScores,
            PSMattr.ShiftedOnlyPositionIons,
            PSMattr.ScoreBestPosition,
            PSMattr.IonsBestPosition,
            PSMattr.ScoreSecondBestPosition,
            PSMattr.IonsSecondBestPosition,
            PSMattr.ScoreAllUnshifted,
            PSMattr.IonsAllUnshifted,
            PSMattr.ScoreShiftedBestPosition,
            PSMattr.IonsShiftedBestPosition,
            PSMattr.ScoreShiftedAllPositions,
            PSMattr.IonsShiftedAllPositions,
            PSMattr.Glycanscore
//            PSMattr.glyco_q_val
    );

} // End of PSM_File class
