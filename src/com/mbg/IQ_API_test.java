package com.mbg;

import ionquant.api.Entry;
import ionquant.api.IonQuantAPI;

public class IQ_API_test {

    public static void main(String[] args) {
        IonQuantAPI api = new IonQuantAPI("C:\\_Local\\_software-tests\\skyline\\cwq_mix2-1_726.mzML", 11, 10, 0.4f, 0.05f, 2, 3, true);

        // build index
        api.buildIndex();

        Entry quantifiedEntry123 = api.quantXIC123(725.3695f, 186.2914f / 60, 0, 2, 0);
        Entry quantifiedEntry = api.quantXIC(725.3695f, 186.2914f / 60, 0, 2, 0);

        IonQuantAPI api2 = new IonQuantAPI("C:\\_Local\\misc_testing\\Wright_glyco-turboID\\UMICH\\raw\\20240727_bkc_25-0700_Wright_Sample-12_sceHCD.raw", 11, 10, 0.4f, 0.05f, 2, 3, true);
        api2.buildIndex();
        Entry quantifiedEntry1232 = api2.quantXIC123(1398.6183f, 5928.4805f / 60, 0, 2, 0);
        Entry quantifiedEntry2 = api2.quantXIC(1398.6183f, 5928.4805f / 60, 0, 2, 0);

    }
}
