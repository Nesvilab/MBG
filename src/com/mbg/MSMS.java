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
    public MXMLReader mReader;

    /**
     * @param psmFile
     * @param msFilePath string path like "Z:\jishen\MBG\data" that storage ms data
     * @throws Exception
     */
    public MSMS(PSM_File psmFile, String msFilePath) throws Exception {
        // init arguments
        HashSet<String> runNames;
        File path = new File(msFilePath);
        // put all ms file names into runName and mappings
        HashSet<String> runNames = new HashSet<>();
        runNames = (HashSet<String>) psmFile.rawfilePSMmap.keySet();
        for (String fs: psmFile.rawfilePSMmap.keySet()) {
            mappings.put(fs, null);
        }
        // fill in mappings
        PSMFile.getMappings(path, mappings, runNames);
    }

    /** read ms file using PTMShepherd
     * @param runName
     * @throws Exception
     */
    public void msFileReader(String runName) throws Exception {
        //
        File cp = mappings.get(runName);
        mReader = new MXMLReader(cp, 7);
        mReader.readFully();
    }

    /** generate all candidates of given range of RT and MZ
     * @param minRT
     * @param maxRT
     * @param minMZ
     * @param maxMZ
     * @return
     */
    public List<Spectrum> spectraSearcher(PSM basepsm, Params params, double minRT, double maxRT, float minMZ, float maxMZ) {
        List<Spectrum> candidates;
        candidates = Arrays.stream(mReader.specs)
                .filter(c -> c.charge >= minRT && c.charge <= maxRT)
                .filter(c -> c.getPrecursorMass() >= minMZ && c.getPrecursorMass() <= maxMZ)
                .collect(Collectors.toList());
        return candidates;
    }
    public void y1Matcher(List<Spectrum> candidates, PSM basepsm, Params params) {
        if (candidates.isEmpty()) {
            return;
        } else {
            // match Y1 ion
            List<Spectrum> qualifiedSpectra = new ArrayList<>();
            double pepMass = basepsm.calcMass - basepsm.glycan.mass;
            double y1Mass = pepMass + params.glycanResidues.get("HexNAc").mass;
            for (Spectrum spec :mReader.specs) {
                if (spec.findIonNeutral(y1Mass, 20, spec.charge) > 0) {
                    qualifiedSpectra.add(spec);
                }
            // if none qualified spectrum

            //
            }
        }


    }


    public static void main(String [] args) throws Exception {
        MXMLReader mr = new MXMLReader(new File("Z:\\jishen\\MBG\\data\\cwq_mix2-1_726_calibrated.mzML"), 7);
        mr.readFully();
        for(int i = 0; i < mr.specs.length; i++)
        	System.out.println(mr.specs[i].scanName);
    }
}
