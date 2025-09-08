package com.mbg;

import com.google.common.base.Joiner;
import com.mbg.msreader.MXMLReader;
import com.mbg.msreader.Spectrum;
import umich.ms.fileio.filetypes.bruker.FrameInfo;
import umich.ms.fileio.filetypes.mzbin.MZBINFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * MSMS handles reading and processing mass spectrometry data files.
 * 
 * This class provides functionality to:
 * 1. Read MS/MS spectra from various file formats
 * 2. Search for specific spectra based on RT, m/z, and IM criteria
 * 3. Extract ion mobility (IM) information from PASEF data
 * 4. Match fragment ions in MS/MS spectra
 */
public class MSMS {
    /** Reader for MS data files */
    public static MXMLReader mReader;
    
    /** Raw file name */
    public String rawFileName;
    
    /** MS data file */
    public File msDataFile;
    
    /** Map of scan numbers to ion mobility values */
    public static HashMap<Integer, Double> ms2IMs;
    
    /** Maximum length of scan numbers for formatting */
    public static int maxScanNumLength;

    public static boolean ifMzBin;

    /**
     * Creates a new MSMS instance for the specified run.
     * 
     * @param runName Name of the MS run (raw file name)
     * @param msFilePath Path to the directory containing MS data files
     * @throws Exception If file cannot be found or accessed
     */
    public MSMS(String runName, String msFilePath) throws Exception {
        rawFileName = runName;
        
        // Attempt to find the MS file with various supported extensions
        File basePath = new File(msFilePath);
        String[] supportedExtensions = {
            ".mzBIN_cache", "_calibrated.mzML", "_uncalibrated.mzML", 
            "_calibrated.mgf", "_uncalibrated.mgf", ".mzML", ".mzXML", 
            ".mzBIN", ".mgf", ".raw"
        };
        
        String filePath = Peak.findFirstPriorityFile(basePath.getParent(), runName, supportedExtensions);
        if (filePath == null) {
            throw new FileNotFoundException(
                String.format("Could not locate MS file %s in directory %s\n", 
                runName, basePath.getParent())
            );
        }
        getIfMzBin(filePath);
        msDataFile = new File(filePath);
    }

    public static void getIfMzBin(String path) {
        ifMzBin = false;
        int lastDot = path.lastIndexOf('.');
        if (lastDot != -1 && lastDot != path.length() - 1) {
            if (path.substring(lastDot + 1).equals("mzBIN")) {
                ifMzBin = true;
            }
        }
    }

    /**
     * Reads the MS data file using PTMShepherd.
     * 
     * @param numThreads Number of threads to use for reading
     * @throws Exception If file cannot be read
     */
    public void msFileReader(int numThreads) throws Exception {
        mReader = new MXMLReader(msDataFile, numThreads);
        mReader.readFully();
        
        // Calculate maximum scan number length for formatting
        maxScanNumLength = String.valueOf(mReader.specs[mReader.specs.length - 1].scanNum).length();
    }

    /**
     * Searches for spectra within specified retention time and m/z ranges.
     * 
     * @param minRT Minimum retention time
     * @param maxRT Maximum retention time
     * @param tMZ target precursor m/z
     * @return List of spectra matching the criteria
     */
    public static List<Spectrum> spectraSearcher(double minRT, double maxRT, double tMZ) {
        // Step 1: Get the initial list of spectra based on RT and mass/mz range
        List<Spectrum> matchedSpectra;
        if (ifMzBin) {
            matchedSpectra = Arrays.stream(mReader.specs)
                    .filter(spectrum -> spectrum.rt >= minRT && spectrum.rt <= maxRT)
                    .filter(spectrum -> spectrum.precursorMass >= (tMZ-1) && spectrum.precursorMass <= (tMZ+1))
                    .collect(Collectors.toList());
        } else {
            matchedSpectra = Arrays.stream(mReader.specs)
                    .filter(spectrum -> spectrum.rt >= minRT && spectrum.rt <= maxRT)
                    .filter(spectrum -> spectrum.mzRangeEnd >= tMZ && spectrum.mzRangeStart <= tMZ)
                    .filter(spectrum -> spectrum.precursorMass >= (tMZ-1) && spectrum.precursorMass <= (tMZ+1))
                    .collect(Collectors.toList());
        }

        // Step 2: Check for the presence of specific scan types in the matched list
        boolean hasHCD = matchedSpectra.stream()
                .anyMatch(spectrum -> Objects.equals(spectrum.scanType, "hcd"));

        boolean hasETDOrETHCD = matchedSpectra.stream()
                .anyMatch(spectrum -> Objects.equals(spectrum.scanType, "etd") || Objects.equals(spectrum.scanType, "ethcd"));

        // Step 3: Apply the conditional filter
        if (hasHCD && hasETDOrETHCD) {
            return matchedSpectra.stream()
                    .filter(spectrum -> Objects.equals(spectrum.scanType, "hcd"))
                    .collect(Collectors.toList());
        } else {
            return matchedSpectra;
        }
    }

    /**
     * Searches for spectra within specified retention time, m/z, and ion mobility ranges.
     * 
     * @param minRT Minimum retention time
     * @param maxRT Maximum retention time
     * @param tMZ target precursor m/z
     * @param minIM Minimum ion mobility
     * @param maxIM Maximum ion mobility
     * @return List of spectra matching the criteria
     */
    public static List<Spectrum> spectraSearcher(double minRT, double maxRT, double tMZ, float minIM, float maxIM) {
        // First filter by ion mobility
        Set<Integer> qualifiedScanNumbers = MSMS.ms2IMs.entrySet().stream()
                .filter(entry -> entry.getValue() >= minIM && entry.getValue() <= maxIM)
                .map(Map.Entry::getKey)
                .collect(Collectors.toSet());

        // Then filter by RT, m/z, and scan numbers
        List<Spectrum> matchedSpectra;
        if (ifMzBin) {
            matchedSpectra = Arrays.stream(mReader.specs)
                    .filter(spectrum -> spectrum.rt >= minRT && spectrum.rt <= maxRT)
                    .filter(spectrum -> spectrum.precursorMass >= (tMZ - 1) && spectrum.precursorMass <= (tMZ + 1))
                    .filter(spectrum -> qualifiedScanNumbers.contains(spectrum.scanNum))
                    .collect(Collectors.toList());
        } else {
            // matched spectrum must within isolation window and 1 mz range
            matchedSpectra = Arrays.stream(mReader.specs)
                    .filter(spectrum -> spectrum.rt >= minRT && spectrum.rt <= maxRT)
                    .filter(spectrum -> spectrum.mzRangeEnd >= tMZ && spectrum.mzRangeStart <= tMZ)
                    .filter(spectrum -> spectrum.precursorMass >= (tMZ-1) && spectrum.precursorMass <= (tMZ+1))
                    .filter(spectrum -> qualifiedScanNumbers.contains(spectrum.scanNum))
                    .collect(Collectors.toList());
        }
        // Step 2: Check for the presence of specific scan types in the matched list
        boolean hasHCD = matchedSpectra.stream()
                .anyMatch(spectrum -> Objects.equals(spectrum.scanType, "hcd"));

        boolean hasETDOrETHCD = matchedSpectra.stream()
                .anyMatch(spectrum -> Objects.equals(spectrum.scanType, "etd") || Objects.equals(spectrum.scanType, "ethcd"));

        // Step 3: Apply the conditional filter
        if (hasHCD && hasETDOrETHCD) {
            return matchedSpectra.stream()
                    .filter(spectrum -> Objects.equals(spectrum.scanType, "hcd"))
                    .collect(Collectors.toList());
        } else {
            return matchedSpectra;
        }
    }

    /**
     * Matches Y0 and Y1 ions in spectra and calculates their relative intensities.
     * 
     * @param qualifiedSpectra List of spectra to search
     * @param basePsm Base PSM providing mass information
     * @param params Parameter object containing glycan residue definitions
     * @return List of relative intensities for Y0+Y1 ions in each spectrum
     */
    public static List<Double> y1Matcher(List<Spectrum> qualifiedSpectra, PSM basePsm, Params params) {
        List<Double> yIonIntensities = new ArrayList<>();
        
        // Calculate theoretical masses for Y0 and Y1 ions
        double peptideMass = basePsm.calcMass - basePsm.glycan.mass;
        double y1Mass = peptideMass + params.glycanResidues.get("HexNAc").mass;
        
        // Process each spectrum
        for (Spectrum spectrum : qualifiedSpectra) {
            // Find base peak intensity for normalization
            double basePeakIntensity = getBasePeakIntensity(spectrum);
            
            if (basePeakIntensity < 1) {
                yIonIntensities.add(0.0);
                continue; // Avoid division by zero
            }
            
            // Find Y0 and Y1 ions
            double y1Intensity = spectrum.findIonNeutral(y1Mass, 20, spectrum.charge);
            double y0Intensity = spectrum.findIonNeutral(peptideMass, 20, spectrum.charge);
            
            // Calculate relative intensity
            double relativeIntensity = (y1Intensity + y0Intensity) / basePeakIntensity;
            yIonIntensities.add(relativeIntensity);
        }
        
        return yIonIntensities;
    }

    /**
     * Helper method to find the base peak intensity in a spectrum.
     */
    private static double getBasePeakIntensity(Spectrum spectrum) {
        double maxIntensity = spectrum.peakInt[0];
        for (double intensity : spectrum.peakInt) {
            if (intensity > maxIntensity) {
                maxIntensity = intensity;
            }
        }
        return maxIntensity;
    }

    /**
     * Extracts ion mobility values for MS/MS spectra from an mzBIN file.
     * 
     * @param mzBinFile mzBIN file containing IM data
     * @param numThreads Number of threads to use for processing
     * @throws Exception If file cannot be read
     */
    public static void getMS2IMs(File mzBinFile, int numThreads) throws Exception {
        ms2IMs = new HashMap<>();
        MZBINFile mzbinSource = new MZBINFile(numThreads, mzBinFile, true);
        
        for (MZBINFile.MZBINSpectrum spectrum : mzbinSource.specs) {
            ms2IMs.put(spectrum.scanNum, (double) spectrum.precursorIM);
        }
    }

    /**
     * Determines the scan mode from a list of frame information objects.
     * 
     * @param frameInfoList List of frame information objects
     * @return Scan mode code
     */
    public static int getScanMode(List<FrameInfo> frameInfoList) {
        // Scan mode codes:
        // 0 = MS, 1 = AutoMSMS, 2 = MRM, 3 = in-source CID, 4 = broadband CID,
        // 8 = PASEF, 9 = DIA, 10 = PRM, 20 = MALDI
        
        // Collect unique scan modes
        Set<Integer> scanModes = new HashSet<>();
        for (FrameInfo frameInfo : frameInfoList) {
            scanModes.add(frameInfo.getScanMode());
        }

        // Handle different cases
        if (scanModes.size() == 1) {
            // Only one scan mode
            return scanModes.iterator().next();
        } else if (scanModes.size() > 1) {
            // Multiple scan modes
            if (scanModes.contains(8) && !scanModes.contains(10)) {
                System.out.println("There are multiple scan modes: " + 
                                  Joiner.on(",").join(scanModes) + 
                                  ". Will only use the scans with mode 8 (PASEF).");
                return 8;
            } else if (!scanModes.contains(8) && scanModes.contains(10)) {
                System.out.println("There are multiple scan modes: " + 
                                  Joiner.on(",").join(scanModes) + 
                                  ". Will only use the scans with mode 10 (PRM).");
                return 10;
            } else if (scanModes.contains(8) && scanModes.contains(10)) {
                System.out.println("There are multiple scan modes: " + 
                                  Joiner.on(",").join(scanModes) + 
                                  ". MSFragger cannot process the data with both 8 (PASEF) and 10 (PRM).");
                System.exit(1);
            } else {
                System.out.println("There are multiple scan modes: " + 
                                  Joiner.on(",").join(scanModes) + 
                                  ". None of them is 8 (PASEF) or 10 (PRM).");
                System.exit(1);
            }
        }

        return -1; // No valid scan mode found
    }

    public static void main(String [] args) throws Exception {
//        MXMLReader mr = new MXMLReader(new File("Z:\\jishen\\MBG\\data\\cwq_mix2-1_726_calibrated.mzML"), 7);
//        mr.readFully();
//        for(int i = 0; i < mr.specs.length; i++)
//        	System.out.println(mr.specs[i].scanName);
        String dfile = "Z:\\dpolasky\\data\\glyco\\MBG\\20230728_TIMSTOF_UM6_5858593_SA_EXT00_newdigest_H2_T1_GB12_1_1450.d";
        getMS2IMs(new File(dfile), 1);
    }
}
