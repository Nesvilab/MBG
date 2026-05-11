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
    private final Params params;
    public static boolean runOPair;

    /**
     * Confident-identification index used by the four-class FDR framework.
     * Outer key: raw file name. Inner: set of "peptide|glycanComposition" strings derived
     * from PSMs whose glycanState == CONFIDENT. Used by type -2 decoy generation (per-file
     * MS2 scan lookup) and by {@link #confidentPsmAtScan}.
     */
    public final HashMap<String, HashSet<String>> confidentIDs = new HashMap<>();

    /**
     * Cross-file confident-identification set: union of all per-file confidentIDs.
     * Used by {@link #isConfidentlyIdentified(String, String)} to classify MBG candidates
     * as type 2 (known true positive anywhere in the experiment) vs type 1 (novel).
     * A glycan identified in any replicate is considered a known true positive even if it
     * was not detected in the specific raw file currently being inferred.
     */
    public final HashSet<String> globalConfidentIDs = new HashSet<>();

    /**
     * Confident-PSM-by-MS2-scan index. Outer key: raw file name. Inner: scan number ->
     * the confident PSM at that scan. Used by type -2 decoy generation to verify that
     * the +k Da shifted m/z lands on an MS2 spectrum that was confidently identified as
     * a different glycopeptide.
     */
    public final HashMap<String, HashMap<Integer, PSM>> ms2ScanToPSM = new HashMap<>();

    /**
     * All-identification scan index: set of every scan number that appears in the PSM
     * table for a given raw file, regardless of glycan q-value or glycan presence.
     * Includes non-glyco identifications.  Used by the relaxed type -2 decoy gate
     * ({@link #isScanIdentified}) to accept any identified precursor at the decoy m/z,
     * not only confident glycopeptides.
     */
    public final HashMap<String, HashSet<Integer>> allIdentifiedScans = new HashMap<>();

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
        columnNames.put(PSMattr.has_Y1, "Y1 ion intensity"); // Was Y1 ion observed in MS2?
        columnNames.put(PSMattr.MBG_score, "MBG Score"); // Score from MBG LDA model
        columnNames.put(PSMattr.mass_error, "Mass error(ppm)"); // Mass error from IonQuant matching
        columnNames.put(PSMattr.KL, "KL divergence"); // KL divergence score from IonQuant
        columnNames.put(PSMattr.shift_frequency, "Shift Frequency"); // Frequency of the observed glycan shift
        columnNames.put(PSMattr.rt_diff, "RT Difference"); // Predicted vs. Observed RT difference
        columnNames.put(PSMattr.im_diff, "IM Difference"); // Predicted vs. Observed IM difference
        columnNames.put(PSMattr.if_new_glycan, "MBG only composition"); // Was this glycan composition novel?
        columnNames.put(PSMattr.MBG_confidence_level, "MBG Confidence Level"); // Confidence category based on the existence of MS2 hit and Y1 ion
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
        // O-pair columns
        columnNames.put(PSMattr.o_pair_score, "O-Pair Score");
        columnNames.put(PSMattr.glycan_site_composition, "Glycan Site Composition(s)");
        columnNames.put(PSMattr.composition_delta_score, "Composition Delta Score");
        columnNames.put(PSMattr.next_best_composition, "Next Best Composition");
        columnNames.put(PSMattr.confidence_level, "Confidence Level");
        columnNames.put(PSMattr.localization_delta_score, "Localization Delta Score");
        columnNames.put(PSMattr.site_probabilities, "Site Probabilities");
        columnNames.put(PSMattr.ratio_138_144, "138/144 Ratio");
        columnNames.put(PSMattr.paired_scan_num, "Paired Scan Num");
        columnNames.put(PSMattr.msfragger_localization, "MSFragger Localization");
        columnNames.put(PSMattr.position_scores, "Position Scores");
        columnNames.put(PSMattr.best_positions, "Best Positions");
        columnNames.put(PSMattr.num_best_positions, "Number Best Positions");
        columnNames.put(PSMattr.shifted_only_position_scores, "Shifted Only Position Scores");
        columnNames.put(PSMattr.shifted_only_position_ions, "Shifted Only Position Ions");
        columnNames.put(PSMattr.score_best_position, "Score Best Position");
        columnNames.put(PSMattr.ions_best_position, "Ions Best Position");
        columnNames.put(PSMattr.score_second_best_position, "Score Second Best Position");
        columnNames.put(PSMattr.ions_second_best_position, "Ions Second Best Position");
        columnNames.put(PSMattr.score_all_unshifted, "Score All Unshifted");
        columnNames.put(PSMattr.ions_all_unshifted, "Ions All Unshifted");
        columnNames.put(PSMattr.score_shifted_best_position, "Score Shifted Best Position");
        columnNames.put(PSMattr.ions_shifted_best_position, "Ions Shifted Best Position");
        columnNames.put(PSMattr.score_shifted_all_positions, "Score Shifted All Positions");
        columnNames.put(PSMattr.ions_shifted_all_positions, "Ions Shifted All Positions");
    }
    // PSM table headers to add to the PSM file if they do not already exist.
    private static final String[] headersToAdd = new String[] {
            // Insert in reverse order of definition to maintain desired final order
            columnNames.get(PSMattr.added_residual),
            columnNames.get(PSMattr.parent_glycan),
            columnNames.get(PSMattr.MBG_confidence_level),
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
        this.params = params;
        this.file = params.psmPath.toFile(); // Store the File object
        this.headers = psmReader.readLine().split("\t"); // Read and split the header line
        this.PSMcolumns = initPSMcolumns(); // Map header names to column indices
        allPSMs = new ArrayList<>(); // Initialize the list to store all PSM objects
        maxRankOfSpectra = new HashMap<>();     // store max rank of all scans in the PSM file to know what rank to use when adding MBG PSMs

        this.runOPair = Arrays.asList(this.headers).contains("O-Pair Score");

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
                        allIdentifiedScans.computeIfAbsent(psm.rawfile, k -> new HashSet<>()).add(psm.scanNum);
                        // Only confidently-identified glyco PSMs contribute to chimeric-MS2 exclusion.
                        // Above-maxq glyco rows have a parsed glycan but should not exclude their scans
                        // (preserving pre-refactor semantics, since they were previously nulled out).
                        if (psm.glycanState == PSM.GlycanState.CONFIDENT) {
                            params.excludeScanNumbersMAp.computeIfAbsent(psm.rawfile, k -> new HashSet<>()).add(psm.scanNum);
                            // Populate the confidentIDs and ms2ScanToPSM indexes used by the
                            // four-class FDR framework (commits 2 and 4 read these).
                            String idKey = psm.peptide + "|" + psm.glycan.toString();
                            confidentIDs
                                .computeIfAbsent(psm.rawfile, k -> new HashSet<>())
                                .add(idKey);
                            globalConfidentIDs.add(idKey);
                            ms2ScanToPSM
                                .computeIfAbsent(psm.rawfile, k -> new HashMap<>())
                                .put(psm.scanNum, psm);
                        }
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
        int confidentCount = confidentIDs.values().stream().mapToInt(Set::size).sum();
        System.out.printf("Indexed %d confident (peptide, glycan) pairs across %d raw files for four-class FDR lookup.%n",
                confidentCount, confidentIDs.size());
    }

    /**
     * Returns true if the original PSM table contains a confident identification of
     * (peptide, glycanComposition) on the given raw file. Used by type -2 decoy
     * generation where a per-file check is required.
     */
    public boolean isConfidentlyIdentified(String rawFile, String peptide, String glycanComposition) {
        HashSet<String> set = confidentIDs.get(rawFile);
        if (set == null) {
            return false;
        }
        return set.contains(peptide + "|" + glycanComposition);
    }

    /**
     * Cross-file overload: returns true if (peptide, glycanComposition) was confidently
     * identified in ANY raw file in this experiment. Used by the four-class FDR type
     * classifier — a glycan observed in another replicate is a known true positive even
     * when it hasn't yet been detected in the file currently being inferred.
     */
    public boolean isConfidentlyIdentified(String peptide, String glycanComposition) {
        return globalConfidentIDs.contains(peptide + "|" + glycanComposition);
    }

    /**
     * Returns true if the given MS2 scan was identified by MSFragger in the original
     * PSM table — regardless of glycan q-value or whether a glycan was assigned.
     * Includes non-glyco peptide identifications.  Used by the relaxed type -2 decoy
     * gate: any identification at the decoy scan is sufficient evidence that the shifted
     * m/z lands on a real fragmented precursor.
     */
    public boolean isScanIdentified(String rawFile, int scanNum) {
        HashSet<Integer> scans = allIdentifiedScans.get(rawFile);
        return scans != null && scans.contains(scanNum);
    }

    /**
     * Returns the confident PSM at the given (raw file, MS2 scan number), or null if
     * no confident PSM was identified at that scan. Used by type -2 decoy generation
     * to verify that a +k Da shifted m/z lands on an MS2 spectrum confidently
     * identified as a different glycopeptide.
     */
    public PSM confidentPsmAtScan(String rawFile, int scanNum) {
        HashMap<Integer, PSM> map = ms2ScanToPSM.get(rawFile);
        if (map == null) {
            return null;
        }
        return map.get(scanNum);
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
                    if (attribute == PSMattr.IM && params.noPASEF) {
                        continue;
                    }
//                    System.err.println("Warning: Could not find column '" + expectedHeader + "' in PSM file. " +
//                            "Processing may fail. Please check the formatting of the PSM file.");
                }
            }
        }
        if (!psmColumnsMap.containsKey(PSMattr.peptide)) {
            System.err.println("WARNING: 'Peptide' column not found in PSM file. Processing may fail.");
        }
        if (!psmColumnsMap.containsKey(PSMattr.spectrum_file)) {
            System.err.println("WARNING: 'Spectrum File' column not found in PSM file. Grouping by raw file will fail.");
        }
        if (!params.noPASEF) {
            // make sure IM column is present if data is set to IM type to avoid a crash later
            if (!psmColumnsMap.containsKey(PSMattr.IM)) {
                MBG.exit("Error: 'Ion Mobility' column not found in PSM file but IM-based processing was requested." +
                        " Please check that the input data type (MS or IM-MS) matches the data", 1);
            }
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
               attr == PSMattr.MBG_confidence_level || attr == PSMattr.parent_glycan ||
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

    private void addOPairCols() {

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
        MBG_confidence_level,

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
        Glycanscore,
        // O-pair
        o_pair_score,
        /** Glycan site composition(s). Header: "Glycan Site Composition(s)" */
        glycan_site_composition,
        /** Composition delta score. Header: "Composition Delta Score" */
        composition_delta_score,
        /** Next best composition. Header: "Next Best Composition" */
        next_best_composition,
        /** Confidence level. Header: "Confidence Level" */
        confidence_level,
        /** Localization delta score. Header: "Localization Delta Score" */
        localization_delta_score,
        /** Site probabilities. Header: "Site Probabilities" */
        site_probabilities,
        /** 138/144 Ratio. Header: "138/144 Ratio" */
        ratio_138_144,
        /** Paired scan number. Header: "Paired Scan Num" */
        paired_scan_num,
        /** MSFragger localization. Header: "MSFragger Localization" */
        msfragger_localization,
        /** Position scores. Header: "Position Scores" */
        position_scores,
        /** Best positions. Header: "Best Positions" */
        best_positions,
        /** Number of best positions. Header: "Number Best Positions" */
        num_best_positions,
        /** Shifted only position scores. Header: "Shifted Only Position Scores" */
        shifted_only_position_scores,
        /** Shifted only position ions. Header: "Shifted Only Position Ions" */
        shifted_only_position_ions,
        /** Score of the best position. Header: "Score Best Position" */
        score_best_position,
        /** Ions of the best position. Header: "Ions Best Position" */
        ions_best_position,
        /** Score of the second best position. Header: "Score Second Best Position" */
        score_second_best_position,
        /** Ions of the second best position. Header: "Ions Second Best Position" */
        ions_second_best_position,
        /** Score of all unshifted positions. Header: "Score All Unshifted" */
        score_all_unshifted,
        /** Ions of all unshifted positions. Header: "Ions All Unshifted" */
        ions_all_unshifted,
        /** Score of the shifted best position. Header: "Score Shifted Best Position" */
        score_shifted_best_position,
        /** Ions of the shifted best position. Header: "Ions Shifted Best Position" */
        ions_shifted_best_position,
        /** Score of all shifted positions. Header: "Score Shifted All Positions" */
        score_shifted_all_positions,
        /** Ions of all shifted positions. Header: "Ions Shifted All Positions" */
        ions_shifted_all_positions;
        ;

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

    public static final EnumSet<PSMattr> PSMattr_blank_OPair = EnumSet.of(
            PSMattr.glycan_site_composition,
            PSMattr.composition_delta_score,
            PSMattr.next_best_composition,
            PSMattr.confidence_level,
            PSMattr.localization_delta_score,
            PSMattr.site_probabilities,
            PSMattr.ratio_138_144,
            PSMattr.paired_scan_num,
            PSMattr.msfragger_localization,
            PSMattr.position_scores,
            PSMattr.best_positions,
            PSMattr.num_best_positions,
            PSMattr.shifted_only_position_scores,
            PSMattr.shifted_only_position_ions,
            PSMattr.score_best_position,
            PSMattr.ions_best_position,
            PSMattr.score_second_best_position,
            PSMattr.ions_second_best_position,
            PSMattr.score_all_unshifted,
            PSMattr.ions_all_unshifted,
            PSMattr.score_shifted_best_position,
            PSMattr.ions_shifted_best_position,
            PSMattr.score_shifted_all_positions,
            PSMattr.ions_shifted_all_positions,
            PSMattr.o_pair_score
    );

} // End of PSM_File class
