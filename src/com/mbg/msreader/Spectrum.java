package com.mbg.msreader;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;

import java.util.*;



import umich.ms.fileio.filetypes.mzbin.MZBINFile;



public class Spectrum {
    public int scanNum;
    public int charge;
    public double precursorMass;
    public double rt;
    double monoMass;
    double targetMass;
    public float[] peakMZ;
    public float[] peakInt;
    double norm;
    public String scanName;
    public int msLevel;
    public double[] averageFragMass;
    public double averageIonMass;
    public double basePeakInt;
    static final float[] fact = new float[128];

    public Spectrum(int nfrags) {
        this.peakMZ = new float[nfrags];
        this.peakInt = new float[nfrags];
        this.norm = -1.0;
        this.basePeakInt = this.findBasePeakInt();
    }

    public Spectrum(String scanname, int scannum, int z, int mslevel, double precursormass, double rettime, float[] peakmz, float[] peakint) {
        this.scanName = MXMLReader.stripChargeState(scanname);
        this.scanNum = scannum;
        this.charge = z;
        this.precursorMass = precursormass;
        this.rt = rettime;
        this.peakMZ = peakmz;
        this.peakInt = peakint;
        this.msLevel = mslevel;
        this.norm = -1.0;
        this.basePeakInt = this.findBasePeakInt();
    }

    public Spectrum(MZBINFile.MZBINSpectrum s, String runName) {
        this.scanNum = s.scanNum;
        this.scanName = runName + "." + this.scanNum + "." + this.scanNum;
        this.charge = s.charge;
        this.msLevel = s.msLevel;
        this.rt = (double)s.retentionTime;
        this.precursorMass = s.selectedMz;
        this.monoMass = s.selectedMz;
        this.targetMass = s.isolatedMz;
        this.peakMZ = Floats.toArray(Doubles.asList(s.peakMZ));
        this.peakInt = s.peakInt;
        this.basePeakInt = this.findBasePeakInt();
        this.norm = -1.0;
    }


    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(this.scanNum + " -");

        for(int i = 0; i < this.peakMZ.length; ++i) {
            sb.append(String.format(" [%.2f %.2f]", this.peakMZ[i], this.peakInt[i]));
        }

        return sb.toString();
    }

    public int compareTo(com.mbg.msreader.Spectrum o) {
        return this.scanNum - o.scanNum;
    }

    public double getSimilarity(com.mbg.msreader.Spectrum s, double ppmTol) {
        double sum = 0.0;
        double ppm = ppmTol * 1.0E-6;
        int pos = 0;

        for(int i = 0; i < this.peakMZ.length; ++i) {
            for(int j = pos; j < s.peakMZ.length; ++j) {
                if ((double)s.peakMZ[j] < (double)this.peakMZ[i] * (1.0 - ppm)) {
                    pos = j;
                } else {
                    if ((double)s.peakMZ[j] > (double)this.peakMZ[i] * (1.0 + ppm)) {
                        break;
                    }

                    sum += (double)(this.peakInt[i] * s.peakInt[j]);
                }
            }
        }

        return sum;
    }

    public void logTransform() {
        for(int i = 0; i < this.peakInt.length; ++i) {
            this.peakInt[i] = (float)Math.log10((double)this.peakInt[i]);
        }

    }

    public double findBasePeakInt() {
        double bpInt = 0.0;

        for(int i = 0; i < this.peakInt.length; ++i) {
            if ((double)this.peakInt[i] > bpInt) {
                bpInt = (double)this.peakInt[i];
            }
        }

        return bpInt;
    }

    public double findIonNeutral(double neutralIonMass, double ppmTol, int maxCharge) {
        double ionIntensity = 0.0;
        int charge = this.charge == 0 ? maxCharge : Math.min(this.charge, maxCharge);

        for(int z = 1; z <= charge; ++z) {
            double ion = (neutralIonMass + 1.00727 * (double)z) / (double)z;
            double ppmRange = ppmTol * 1.0E-6 * ion;
            double max = ion + ppmRange;

            for(int i = 0; i < this.peakMZ.length && !((double)this.peakMZ[i] > max); ++i) {
                if (Math.abs((double)this.peakMZ[i] - ion) < ppmRange) {
                    ionIntensity += (double)this.peakInt[i];
                }
            }
        }

        return ionIntensity;
    }


}
