package com.mbg;

import java.util.ArrayList;
import java.util.HashMap;

public class GlycoSite {

    ArrayList<PSM> psms;
    HashMap<String, ArrayList<PSM>> observedGlycans;
    String protein;
    int siteIndex;
    String id;
    PSM basePSM;

    public GlycoSite(PSM firstPSM, String id) {
        this.protein = firstPSM.protein;
        this.id = id;
        psms = new ArrayList<>();
        psms.add(firstPSM);
        basePSM = firstPSM;
        observedGlycans = new HashMap<>();
        ArrayList<PSM> glycoPSMs = new ArrayList<>();
        glycoPSMs.add(firstPSM);
        observedGlycans.put(firstPSM.glycan.name, glycoPSMs);
    }

    public void addPSM(PSM psm) {
        psms.add(psm);
        String glycanComp = psm.glycan.name;
        ArrayList<PSM> glycoPSMs;
        if (observedGlycans.containsKey(glycanComp)) {
            glycoPSMs = observedGlycans.get(glycanComp);
        } else {
            glycoPSMs = new ArrayList<>();
        }
        glycoPSMs.add(psm);
        observedGlycans.put(glycanComp, glycoPSMs);
    }

}
