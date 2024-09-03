package com.mbg;

import ionquant.api.Entry;
import ionquant.api.IonQuantAPI;

import java.io.File;
import java.util.Objects;

public final class Peak {
    public static IonQuantAPI api;
    public float averageMass;
    public Entry quantifiedEntry123;

    public Peak() {
    }

    public static void indexBuilder(String path) {
        api = new IonQuantAPI(path, 11, Params.mztol, Params.rtTol, Params.imTol, 2, 3, Params.noPASEF);
        api.buildIndex();
    }

    public void matchMS1(float mz, float rt, float im, int charge, int cv) {
        quantifiedEntry123 = null;
        Entry quantifiedEntry0 = api.quantXIC((float) (mz-1.002), rt, im, charge, cv);
        Entry quantifiedEntry1 = api.quantXIC(mz, rt, im, charge, cv);
        Entry quantifiedEntry2 = api.quantXIC((float) (mz+1.002), rt, im, charge, cv);
        Entry quantifiedEntry3 = api.quantXIC((float) (mz+2.004), rt, im, charge, cv);
        if (quantifiedEntry0 == null) {
            // remove peaks if detect M-1 peak
            quantifiedEntry123 = api.quantXIC123(mz, rt, im, charge, cv);
        }
        if (quantifiedEntry1 == null || quantifiedEntry2 == null || quantifiedEntry3 == null) {
            averageMass = 0;
        } else {
            getAverageMass(quantifiedEntry1, quantifiedEntry2, quantifiedEntry3, mz);
        }
    }

    private void getAverageMass(Entry entry1, Entry entry2, Entry entry3, float mz) {
        float intensitySum = entry1.intensity + entry2.intensity + entry3.intensity;
        averageMass = (float) (((entry1.intensity*mz)+(entry2.intensity*(mz+1.002))+(entry3.intensity*(mz+2.004))) / intensitySum);
    }

    public static String findFirstPriorityFile(String directoryPath, String filePrefix) {
        File directory = new File(directoryPath);

        // Define priority order for extensions
        String[] extensions = {".mzML", ".raw", ".d"};

        // Iterate through files in the directory
        for (String extension : extensions) {
            for (File file : Objects.requireNonNull(directory.listFiles())) {
                String fileNameWithoutExtension = file.getName();
                int lastIndexOfDot = fileNameWithoutExtension.lastIndexOf('.');
                if (lastIndexOfDot != -1) {
                    fileNameWithoutExtension = fileNameWithoutExtension.substring(0, lastIndexOfDot);
                } else {
                    continue;
                }
                if (fileNameWithoutExtension.equals(filePrefix)) {
                    // Check for each extension in priority order
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


