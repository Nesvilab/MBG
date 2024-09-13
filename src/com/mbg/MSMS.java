package com.mbg;

import com.google.common.base.Joiner;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import umich.ms.fileio.filetypes.bruker.FrameInfo;
import umich.ms.fileio.filetypes.mzbin.MZBINFile;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.Collectors;

public class MSMS {
    public static MXMLReader mReader;
    public String rawFileName;
    public File msDir;
    public static HashMap<Integer, Double> ms2IMs;
    public static int maxScanNumLength;

    /**
     * @param runName
     * @param msFilePath string path like "Z:\jishen\MBG\data" that storage ms data
     * @throws Exception
     */
    public MSMS(String runName, String msFilePath) throws Exception {
        // init arguments
        rawFileName = runName;
        File path = new File(msFilePath);
        String[] extensions = {".mzBIN_cache", "_calibrated.mzML", "_uncalibrated.mzML", "_calibrated.mgf", "_uncalibrated.mgf", ".mzML", ".mzXML", ".mzBIN", ".mgf", ".raw"};
        String filePath = Peak.findFirstPriorityFile(path.getParent(), runName, extensions);
        if (filePath == null) {
            throw new FileNotFoundException(String.format("Could not locate MS file %s in directory %s\n", runName, path.getParent()));
        }
        msDir = new File(filePath);
    }

    /** read ms file using PTMShepherd
     * @throws Exception
     */
    public void msFileReader(int numThreads) throws Exception {
        mReader = new MXMLReader(msDir, numThreads);
        mReader.readFully();
        // get the length of spectrum number of the last spectrum
        maxScanNumLength = String.valueOf(mReader.specs[mReader.specs.length - 1].scanNum).length();
    }

    /** generate all candidates of given range of RT and MZ
     * @param minRT
     * @param maxRT
     * @param minMZ
     * @param maxMZ
     */
    public static List<Spectrum> spectraSearcher(double minRT, double maxRT, float minMZ, float maxMZ) {
        List<Spectrum> candidates;
        candidates = Arrays.stream(mReader.specs)
                .filter(c -> c.rt >= minRT && c.rt <= maxRT)
                .filter(c -> c.precursorMass >= minMZ && c.precursorMass <= maxMZ)
                .collect(Collectors.toList());
        return candidates;
    }

    public static List<Spectrum> spectraSearcher(double minRT, double maxRT, float minMZ, float maxMZ, float minIM, float maxIM) {
        List<Spectrum> candidates;
        Set<Integer> qualifiedScanNumber = MSMS.ms2IMs.entrySet().stream()
                .filter(entry -> entry.getValue() >= minIM && entry.getValue() <= maxIM)
                .map(Map.Entry::getKey)
                .collect(Collectors.toSet());

        candidates = Arrays.stream(mReader.specs)
                .filter(c -> c.rt >= minRT && c.rt <= maxRT)
                .filter(c -> c.precursorMass >= minMZ && c.precursorMass <= maxMZ)
                .filter(c -> qualifiedScanNumber.contains(c.scanNum))
                .collect(Collectors.toList());
        return candidates;
    }

    public static List<Boolean> y1Matcher(List<Spectrum> qualifiedSpectra, PSM basepsm, Params params) {
        List<Boolean> hasY1 = new ArrayList<>();
        // match Y1 ion
        double pepMass = basepsm.calcMass - basepsm.glycan.mass;
        double y1Mass = pepMass + params.glycanResidues.get("HexNAc").mass;
        for (Spectrum spec :qualifiedSpectra) {
            if (spec.findIonNeutral(y1Mass, 20, spec.charge) > 0) {
                hasY1.add(true);
            } else {
                hasY1.add(false);
            }
        }
        return hasY1;
    }

    public static void getMS2IMs(File pathDir, int nThreads) throws Exception {
        ms2IMs = new HashMap<>();
        MZBINFile mzbinSource = new MZBINFile(nThreads, pathDir, true);
        for (MZBINFile.MZBINSpectrum spec: mzbinSource.specs) {
            ms2IMs.put(spec.scanNum, (double) spec.precursorIM);
        }
    }

    public static int getScanMode(List<FrameInfo> frameInfoList) { // 0 = MS, 1 = AutoMSMS, 2 = MRM, 3 = in-source CID, 4 = broadband CID, 8 = PASEF, 9 = DIA, 10 = PRM, 20 = MALDI
        Set<Integer> scanModes = new HashSet<>();
        for (FrameInfo frameInfo : frameInfoList) {
            scanModes.add(frameInfo.getScanMode());
        }

        if (scanModes.size() == 1) {
            return scanModes.iterator().next();
        } else if (scanModes.size() > 1) {
            if (scanModes.contains(8) && !scanModes.contains(10)) {
                System.out.println("There are multiple scan modes: " + Joiner.on(",").join(scanModes) + ". Will only use the scans with mode 8 (PASEF).");
                return 8;
            } else if (!scanModes.contains(8) && scanModes.contains(10)) {
                System.out.println("There are multiple scan modes: " + Joiner.on(",").join(scanModes) + ". Will only use the scans with mode 10 (PRM).");
                return 10;
            } else if (scanModes.contains(8) && scanModes.contains(10)) {
                System.out.println("There are multiple scan modes: " + Joiner.on(",").join(scanModes) + ". MSFragger cannot process the data with both 8 (PASEF) and 10 (PRM).");
                System.exit(1);
            } else {
                System.out.println("There are multiple scan modes: " + Joiner.on(",").join(scanModes) + ". None of them is 8 (PASEF) or 10 (PRM).");
                System.exit(1);
            }
        }

        return -1;
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
