package com.mbg;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

public class MSMS {
    public Path msPath;
    public HashMap<String, File> mappings = new HashMap<>();  // cwq_mix2-1_726 -> {File@993} "Z:\jishen\MBG\data\cwq_mix2-1_726_calibrated.mzML"
    public static MXMLReader mReader;
    public String rawFileName;
    public File msDir;

    /**
     * @param runName
     * @param msFilePath string path like "Z:\jishen\MBG\data" that storage ms data
     * @throws Exception
     */
    public MSMS(String runName, String msFilePath) throws Exception {
        // init arguments
        rawFileName = runName;
        HashSet<String> runNames = new HashSet<>();
        runNames.add(runName);
        mappings.put(runName, null);
        File path = new File(msFilePath);
        // put all ms file names into runName and mappings
//        runNames = (HashSet<String>) psmFile.rawfilePSMmap.keySet();
//        for (String fs: psmFile.rawfilePSMmap.keySet()) {
//            mappings.put(fs, null);
//        }
        // fill in mappings
        PSMFile.getMappings(path, mappings, runNames);
    }

    /** read ms file using PTMShepherd
     * @throws Exception
     */
    public void msFileReader() throws Exception {
        //
        msDir = mappings.get(rawFileName);
        mReader = new MXMLReader(msDir, 7);
        mReader.readFully();
    }

    /** generate all candidates of given range of RT and MZ
     * @param minRT
     * @param maxRT
     * @param minMZ
     * @param maxMZ
     * @return
     */
    public static List<Spectrum> spectraSearcher(double minRT, double maxRT, float minMZ, float maxMZ) {
        List<Spectrum> candidates;
        candidates = Arrays.stream(mReader.specs)
                .filter(c -> c.rt >= minRT && c.rt <= maxRT)
                .filter(c -> c.precursorMass >= minMZ && c.precursorMass <= maxMZ)
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


    public static void main(String [] args) throws Exception {
        MXMLReader mr = new MXMLReader(new File("Z:\\jishen\\MBG\\data\\cwq_mix2-1_726_calibrated.mzML"), 7);
        mr.readFully();
        for(int i = 0; i < mr.specs.length; i++)
        	System.out.println(mr.specs[i].scanName);
    }
}
