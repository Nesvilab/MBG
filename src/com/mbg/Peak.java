package com.mbg;

import ionquant.api.Entry;
import ionquant.api.IonQuantAPI;

import java.io.File;
import java.util.Objects;

/**
 * Peak handles MS1 feature detection and quantification using IonQuant.
 * 
 * This class provides functionality to:
 * 1. Initialize and build IonQuant indices for MS data files
 * 2. Match and quantify features in MS1 scans
 * 3. Calculate average masses from isotope patterns
 * 4. Find MS data files with prioritized extensions
 */
public final class Peak {
    /** IonQuant API instance for MS1 feature matching */
    public static IonQuantAPI api;
    
    /** Average mass calculated from isotope pattern */
    public float averageMass;
    
    /** Quantified feature entry from IonQuant */
    public Entry quantifiedEntry;

    /**
     * Creates a new Peak instance.
     */
    public Peak() {
        // Default constructor
    }

    /**
     * Initializes IonQuant and builds the feature index for a mass spectrometry file.
     * 
     * @param filePath Path to the MS data file
     * @param params Parameters controlling the feature detection
     */
    public static void indexBuilder(String filePath, Params params) {
        api = new IonQuantAPI(
            filePath, 
            params.numThreads, 
            params.mztol, 
            params.rtTol, 
            params.imTol, 
            params.minIsotopesIonQuant, 
            params.minScansIonQuant, 
            params.noPASEF
        );
        api.buildIndex();
    }

    /**
     * Matches an MS1 feature based on m/z, retention time, ion mobility, and charge.
     * Also calculates the average mass from the isotope pattern if detected.
     * 
     * @param mz Precursor m/z
     * @param rt Retention time (in minutes)
     * @param im Ion mobility
     * @param charge Charge state
     * @param cv Collision voltage (for FAIMS data)
     */
    public void matchMS1(float mz, float rt, float im, int charge, int cv) {
        quantifiedEntry = null;
        
        // Check for M-1 isotope to filter out potential interferences
        Entry isotopeMinus1 = api.quantXIC((float) (mz - 1.002), rt, im, charge, cv);
        Entry monoisotopicPeak = api.quantXIC(mz, rt, im, charge, cv);
        Entry isotopePlus1 = api.quantXIC((float) (mz + 1.002), rt, im, charge, cv);
        Entry isotopePlus2 = api.quantXIC((float) (mz + 2.004), rt, im, charge, cv);
        
        if (isotopeMinus1 == null) {
            // No M-1 peak detected, use regular quantification
            quantifiedEntry = api.quantXIC123(mz, rt, im, charge, cv);
        }
        
        // Calculate average mass if isotope pattern is complete
        if (monoisotopicPeak == null || isotopePlus1 == null || isotopePlus2 == null) {
            averageMass = 0;
        } else {
            calculateAverageMass(monoisotopicPeak, isotopePlus1, isotopePlus2, mz);
        }
    }

    /**
     * Calculates the intensity-weighted average mass from an isotope pattern.
     * 
     * @param monoisotopicPeak The monoisotopic peak (M)
     * @param isotopePlus1 The M+1 isotope peak
     * @param isotopePlus2 The M+2 isotope peak
     * @param baseMz The monoisotopic m/z value
     */
    private void calculateAverageMass(Entry monoisotopicPeak, Entry isotopePlus1, Entry isotopePlus2, float baseMz) {
        float intensitySum = monoisotopicPeak.intensity + isotopePlus1.intensity + isotopePlus2.intensity;
        
        // Calculate weighted average mass
        averageMass = (float) (
            (monoisotopicPeak.intensity * baseMz) + 
            (isotopePlus1.intensity * (baseMz + 1.002)) + 
            (isotopePlus2.intensity * (baseMz + 2.004))
        ) / intensitySum;
    }

    /**
     * Searches for a file with the specified prefix and one of the provided extensions,
     * prioritizing extensions in the order they appear in the array.
     * 
     * @param directoryPath Directory to search
     * @param filePrefix Prefix of the file to find
     * @param extensions Array of extensions to search for, in priority order
     * @return Full path to the first matching file, or null if none found
     */
    public static String findFirstPriorityFile(String directoryPath, String filePrefix, String[] extensions) {
        File directory = new File(directoryPath);
        File[] files = directory.listFiles();
        
        if (files == null) {
            return null;
        }

        // Search for each extension in priority order
        for (String extension : extensions) {
            for (File file : files) {
                String fileName = file.getName();
                
                // Check if this file matches our criteria
                if (fileName.endsWith(extension) && 
                    fileName.substring(0, fileName.length() - extension.length()).equals(filePrefix)) {
                    return file.getAbsolutePath();
                }
            }
        }

        return null; // No matching file found
    }

    public static void main(String[] args) {
        IonQuantAPI api = new IonQuantAPI("C:\\Users\\jishen\\Downloads\\cwq_mix2-1_726.raw", 11, 10, 0.4f, 0.05f, 2, 3, true);
//        IonQuantAPI api = new IonQuantAPI("Z:\\dpolasky\\data\\glyco\\pGlyco2\\brain_PXD005411\\MouseBrain-Z-T-1.mzML", 11, 10, 0.4f, 0.05f, 2, 3, true);

        api.buildIndex();
        Entry quantifiedEntry123 = api.quantXIC123(410.7354f, 1636.7035f / 60, 0, 2, 0);
//        Entry quantifiedEntry = api.quantXIC(410.7354f, 1636.7035f / 60, 0, 2, 0);
//        String directoryPath = "Z:\\jishen\\MBG\\data";
//        String filePrefix = "cwq_mix2-1_726";
//
//        String foundFile = findFirstPriorityFile(directoryPath, filePrefix);
//
//        if (foundFile != null) {
//            System.out.println("Found file: " + foundFile);
//        } else {
//            System.out.println("No file found with the specified prefix and extensions.");
//        }
    }
}


