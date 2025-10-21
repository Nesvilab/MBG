package com.mbg;

import umich.ms.glyco.Glycan;

import java.io.*;
import java.util.*;

public class Utils {

    public static final String ionQuantMassPath = "modmasses_ionquant.txt"; // File potentially used by IonQuant
    public static final String newGlycanPath = "mbg_added_glycans.txt"; // File potentially used by IonQuant
    public static final String rtshiftsPath = "rt_shifts.csv";
    public static final String imshiftsPath = "im_shifts.csv";

    public static void writeSkylineMods(HashSet<String> addedGlycans) {
        if (addedGlycans.isEmpty()) {
            return;
        }
        String output = String.join(",", addedGlycans);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(newGlycanPath))) {
            writer.write(output);
            System.out.printf("Wrote %d new glycans to %s%n", addedGlycans.size(), newGlycanPath);
        } catch (IOException e) {
            System.err.printf("An error occurred while writing new glycans to %s: %s%n", newGlycanPath, e.getMessage());
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
                    targetWriter.write(String.format("%.4f", mass));
                    targetWriter.newLine();
                    appendedCount++;
                }
            }
        }
        System.out.printf("Appended %d unique masses to %s.%n", appendedCount, targetFilePath);
    }
    public static void writeRtImShifts(HashMap<String, List<Float>> rtShifts, HashMap<String, List<Float>> imShifts, Boolean noPASEF, ArrayList<Glycan> toAddResidues) {
        // write rt shifts
        List<String> toAddResiduesStrings = new ArrayList<>();

        for (Glycan g : toAddResidues) {
            toAddResiduesStrings.add(g.toString());
        }
        // filter rt/im shifts
        HashMap<String, List<Float>> rtShifts_filtered = new HashMap<>();
        for (String key : toAddResiduesStrings) {
            if (rtShifts.containsKey(key)) {
                rtShifts_filtered.put(key, rtShifts.get(key));
            }
        }

//        writeHashMapToCSV(rtShifts_filtered, rtshiftsPath);
        writeHashMapToCSV(rtShifts, rtshiftsPath);
        if (!noPASEF) {
            HashMap<String, List<Float>> imShifts_filtered = new HashMap<>();
            for (String key : toAddResiduesStrings) {
                if (imShifts.containsKey(key)) {
                    imShifts_filtered.put(key, imShifts.get(key));
                }
            }
//            writeHashMapToCSV(imShifts_filtered, imshiftsPath);
            writeHashMapToCSV(imShifts, imshiftsPath);
        }
    }

    private static void writeHashMapToCSV(HashMap<String, List<Float>> data, String filePath) {
        try (FileWriter writer = new FileWriter(filePath)) {
            // Write header
            writer.write("glycan_differences,shift\n");

            // Sort keys for consistent output
            List<String> sortedKeys = new ArrayList<>(data.keySet());
            Collections.sort(sortedKeys);

            // Write each key-value pair as multiple rows
            for (String key : sortedKeys) {
                List<Float> shifts = data.get(key);
                for (Float value : shifts) {
                    writer.write(key + "," + value + "\n");
                }
            }

            System.out.println("CSV file written in long format: " + filePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
