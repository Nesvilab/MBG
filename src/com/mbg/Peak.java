package com.mbg;

import ionquant.api.Entry;
import ionquant.api.IonQuantAPI;

public final class Peak {
    public static IonQuantAPI api;

    private Peak() {
    }

    public static void indexBuilder(String path) {
        api = new IonQuantAPI(path, 11, 10, 0.4f, 0.05f, 2, 3, true);
        api.buildIndex();
    }

    public static void matchMS1(float mz, float rt, float im, int charge, int cv) {
        Entry quantifiedEntry;
        quantifiedEntry = api.quantXIC(mz, rt, im, charge, cv);
    }


    public static void main(String[] args) {
        IonQuantAPI api = new IonQuantAPI("C:\\Users\\jishen\\Downloads\\cwq_mix2-1_726.raw", 11, 10, 0.4f, 0.05f, 2, 3, true);
//        IonQuantAPI api = new IonQuantAPI("Z:\\dpolasky\\data\\glyco\\pGlyco2\\brain_PXD005411\\MouseBrain-Z-T-1.mzML", 11, 10, 0.4f, 0.05f, 2, 3, true);

        api.buildIndex();
        Entry quantifiedEntry123 = api.quantXIC123(410.7354f, 1636.7035f / 60, 0, 2, 0);
        Entry quantifiedEntry = api.quantXIC(410.7354f, 1636.7035f / 60, 0, 2, 0);
    }
}


