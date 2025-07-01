package com.mbg.msreader;

import com.mbg.msreader.Spectrum;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class MSFMGFFile {
    private static final int MSFMGFVersion = 20200824;
    public final File f;
    public final String runName;
    public ArrayList<com.mbg.msreader.Spectrum> specs;
    public ArrayList<Integer> scanNums;
    public int maxScan;
    private int runStartTime;
    private String mzMLHeader;
    private String rawSHA1;
    private String activationMethods;
    private String batmassVersion;
    private String timsLibVersion;

    public MSFMGFFile(ExecutorService executorService, int nThread, String filePath, boolean loadScans) throws Exception {
        this(executorService, nThread, new File(filePath), loadScans);
    }

    public MSFMGFFile(ExecutorService executorService, int nThread, File f, boolean loadScans) throws Exception {
        this.runStartTime = 0;
        this.mzMLHeader = "";
        this.rawSHA1 = "";
        this.activationMethods = "";
        this.batmassVersion = "";
        this.timsLibVersion = "";
        this.f = f;
        this.runName = this.getBasename(f.getName());
        this.specs = new ArrayList();
        this.scanNums = new ArrayList();
        this.readMSFMGF(executorService, nThread, loadScans);
    }

    public static String getBasenameStat(String f) {
        String baseName;
        if (f.contains("_calibrated")) {
            baseName = f.substring(0, f.indexOf("_calibrated"));
        } else if (f.contains("_uncalibrated")) {
            baseName = f.substring(0, f.indexOf("_uncalibrated"));
        } else {
            baseName = f.substring(0, f.lastIndexOf("."));
        }

        return baseName;
    }

    public String getBasename(String f) {
        String baseName;
        if (f.contains("_calibrated")) {
            baseName = f.substring(0, f.indexOf("_calibrated"));
        } else if (f.contains("_uncalibrated")) {
            baseName = f.substring(0, f.indexOf("_uncalibrated"));
        } else {
            baseName = f.substring(0, f.lastIndexOf("."));
        }

        return baseName;
    }

    private void readMSFMGF(ExecutorService executorService, int nThread, boolean loadScans) throws Exception {
        this.loadStructure();
        if (loadScans) {
            this.loadMSFMGFScans(executorService, nThread);
        }

    }

    private void loadStructure() throws Exception {
        this.scanNums = new ArrayList();
        BufferedReader br = new BufferedReader(new FileReader(this.f));
        boolean readFlag = false;
        int bScan = 0;
        int eScan = 0;

        String cline;
        while((cline = br.readLine()) != null) {
            cline = cline.trim();
            if (cline.equals("BEGIN IONS")) {
                ++bScan;
                readFlag = true;
            } else if (cline.equals("END IONS")) {
                ++eScan;
            }

            if (readFlag && cline.substring(0, 4).equals("TITLE")) {
                this.scanNums.add(this.getScanNum(cline));
                readFlag = false;
            }
        }

        this.maxScan = Math.max(bScan, eScan);
    }

    private int getScanNum(String line) {
        String scanStr = line.trim().split("=")[1].trim();
        int scanNum = Integer.parseInt(scanStr.split(".")[1]);
        return scanNum;
    }

    private void loadMSFMGFScans(ExecutorService executorService, int nThread) throws Exception {
        this.scanNums = new ArrayList();
        this.specs = new ArrayList();
        BufferedReader br = new BufferedReader(new FileReader(this.f));
        int nBlocks = this.maxScan / 500;
        if (this.maxScan % 500 != 0) {
            ++nBlocks;
        }

        ArrayList<String> fileLines = new ArrayList();
        int[] scanStarts = new int[this.maxScan];
        int scanC = 0;

        String cline;
        for(int lineC = 0; (cline = br.readLine()) != null; ++lineC) {
            fileLines.add(cline);
            if (cline.equals("BEGIN IONS")) {
                scanStarts[scanC] = lineC;
            } else if (cline.equals("END IONS")) {
                ++scanC;
            }
        }

        ArrayList<Future> futureList = new ArrayList(nBlocks);

        for(int i = 0; i < nBlocks; ++i) {
            int istart = i * 500;
            int iend = (i + 1) * 500;
            int ilinestart = scanStarts[istart];
            int ilineend;
            if (iend >= scanC) {
                ilineend = fileLines.size();
            } else {
                ilineend = scanStarts[iend];
            }

            futureList.add(executorService.submit(() -> {
                this.processMSFMGFSpectraBlock(fileLines, ilinestart, ilineend);
            }));
        }

        Iterator var17 = futureList.iterator();

        while(var17.hasNext()) {
            Future future = (Future)var17.next();
            future.get();
        }

    }

    private void processMSFMGFSpectraBlock(ArrayList<String> lines, int istart, int iend) {
        ArrayList<com.mbg.msreader.Spectrum> spectra = new ArrayList();
        ArrayList<String> cspecLines = new ArrayList();

        for(int i = istart; i < iend; ++i) {
            String line = (String)lines.get(i);
            cspecLines.add(line);
            if (line.equals("END IONS")) {
                try {
                    spectra.add(this.loadMSFMGFSpectrum(cspecLines));
                } catch (Exception var10) {
                    Exception e = var10;
                    System.out.printf("Malformed scan: %s\n", cspecLines.get(1));
                    e.printStackTrace();
                }

                cspecLines.clear();
            }
        }

        this.addSpecs(spectra);
    }

    private synchronized void addSpecs(ArrayList<com.mbg.msreader.Spectrum> specs) {
        this.specs.addAll(specs);
        Iterator var2 = specs.iterator();

        while(var2.hasNext()) {
            com.mbg.msreader.Spectrum spec = (com.mbg.msreader.Spectrum)var2.next();
            this.scanNums.add(spec.scanNum);
        }

    }

    private int getScanNumFromTitle(String title) {
        String[] bits = title.trim().split("\\.");
        int scanNum;
        if (bits.length == 4) {
            scanNum = Integer.parseInt(bits[2]);
        } else {
            scanNum = Integer.parseInt(bits[1]);
        }

        return scanNum;
    }

    private String getScanNameFromTitle(int scanNum) {
        String scanname = String.format("%s.%d.%d", this.runName, scanNum, scanNum);
        return scanname;
    }

    private Spectrum loadMSFMGFSpectrum(ArrayList<String> block) throws Exception {
        String scanName = "";
        int scanNum = -1;
        int charge = 0;
        int msLevel = 2;
        double precursorMZ = -1.0;
        double retTime = -1.0;
        ArrayList<Float> peakmz = new ArrayList();
        ArrayList<Float> peakint = new ArrayList();

        try {
            Iterator var19 = block.iterator();

            while(var19.hasNext()) {
                String line = (String)var19.next();
                if (!line.equals("BEGIN IONS") && !line.equals("END IONS")) {
                    String[] cinfo;
                    if (line.contains("=")) {
                        cinfo = line.trim().split("=");
                        String att = cinfo[0];
                        String val = cinfo[1];
                        if (att.equals("TITLE")) {
                            scanNum = this.getScanNumFromTitle(val);
                            scanName = this.getScanNameFromTitle(scanNum);
                        } else if (att.equals("CHARGE")) {
                            charge = Integer.parseInt(val.substring(0, val.indexOf("+")));
                        } else if (att.equals("RTINSECONDS")) {
                            retTime = Double.parseDouble(val);
                        } else if (att.equals("PEPMASS")) {
                            precursorMZ = Double.parseDouble(val.split(" ")[0]);
                        }
                    } else {
                        cinfo = line.trim().split(" ");
                        peakmz.add(Float.parseFloat(cinfo[0]));
                        peakint.add(Float.parseFloat(cinfo[1]));
                    }
                }
            }
        } catch (Exception var18) {
            Exception e = var18;
            System.out.printf("Malformed scan: %s\n", block.get(1));
            e.printStackTrace();
            System.exit(1);
        }

        float[] peakmzArr = new float[peakmz.size()];

        for(int i = 0; i < peakmz.size(); ++i) {
            peakmzArr[i] = (Float)peakmz.get(i);
        }

        float[] peakintArr = new float[peakint.size()];

        for(int i = 0; i < peakmz.size(); ++i) {
            peakintArr[i] = (Float)peakint.get(i);
        }

        return new Spectrum(scanName, scanNum, charge, msLevel, precursorMZ, retTime, peakmzArr, peakintArr);
    }
}
