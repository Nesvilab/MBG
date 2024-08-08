package com.mbg;

import ionquant.api.Entry;
import ionquant.api.IonQuantAPI;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public final class Peak {
    public static IonQuantAPI api;
    public Entry quantifiedEntry;
    public Entry quantifiedEntry123;

    public Peak() {
    }

    public static void indexBuilder(String path) {
        api = new IonQuantAPI(path, 11, 10, 0.4f, 0.05f, 2, 3, true);
        api.buildIndex();
    }

    public void matchMS1(float mz, float rt, float im, int charge, int cv) {
        quantifiedEntry = api.quantXIC(mz, rt, im, charge, cv);
        quantifiedEntry123 = api.quantXIC123(mz, rt, im, charge, cv);
    }

    public void matchMS1forPSMs(HashMap<String, GlycoSite> observedGlycoSites) {
        for (Map.Entry<String, GlycoSite> siteEntry : observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            for (Map.Entry<String, ArrayList<PSM>> glycanEntry: site.observedGlycans.entrySet()) {
                ArrayList<PSM> glycan_psm = glycanEntry.getValue();
                String glycanName = glycanEntry.getKey();

                Map<String, Float> glycanRTs = new HashMap<>();
                float topApexRT = 0;  // get top intensity apex RT (in minute) of this glycopeptide
                float topIntensity = 0;
                for (PSM psm : glycan_psm) {
                    matchMS1((float) psm.calcMz, (float) psm.rt / 60, 0, psm.charge, 0);
                    if (quantifiedEntry123 != null) {
                        psm.apexRT = quantifiedEntry123.apexRt;
                        psm.intensity = quantifiedEntry123.intensity;
                        if (psm.intensity > topIntensity) {
                            topApexRT = (float) psm.apexRT;
                            topIntensity = (float) psm.intensity;

                        }
                    }
                }
                if (topApexRT > 0) {
                    glycanRTs.put(glycanName, topApexRT);
                }
            // storage into RT adjuster for calculating RT median shift

            }
        }
    }



    public static String findFirstPriorityFile(String directoryPath, String filePrefix) {
        File directory = new File(directoryPath);

        // Define priority order for extensions
        String[] extensions = {".mzML", ".raw", ".d"};

        // Iterate through files in the directory
        for (File file : Objects.requireNonNull(directory.listFiles())) {
            if (file.isFile() && file.getName().startsWith(filePrefix)) {
                // Check for each extension in priority order
                for (String extension : extensions) {
                    if (file.getName().endsWith(extension)) {
                        return file.getAbsolutePath(); // Return the first match
                    }
                }
            }
        }

        return null; // No file found with the specified criteria
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


