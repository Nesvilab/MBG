package com.mbg;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import umich.ms.glyco.Glycan;

public class GlycanMatching {
    public PSM_File psmFile;
    public ArrayList<PSM> psms;
    private final Params params;
    private HashMap<String, GlycoSite> observedGlycoSites;

    public GlycanMatching(Params params, PSM_File psmFile, ArrayList<PSM> psms) throws Exception {
        this.psmFile = psmFile;
        this.psms = psms;
        this.params = params;
    }

    public void groupGlycosites() {
        observedGlycoSites = new HashMap<>();

        for (PSM psm: psms) {
            if (psm.glycan == null) {
                continue;
            }
            String siteID = psm.getGlycoSiteID();
            if (observedGlycoSites.containsKey(siteID)) {
                observedGlycoSites.get(siteID).addPSM(psm);
            } else {
                // new site
                observedGlycoSites.put(siteID, new GlycoSite(psm, siteID));
            }
        }
    }

    /**
     * Given a glycan database and rules from Params, fill in the unobserved glycans for each glycosite
     * and generate dummy PSMs for them.
     * @return
     */
    public ArrayList<PSM> fillInGlycoSites() {
        ArrayList<PSM> newPSMs = new ArrayList<>();

        for (Map.Entry<String, GlycoSite> siteEntry: observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            int rankCounter = 1;
            if (site.psms.size() >= params.minPSMs && site.observedGlycans.size() >= params.minGlycans) {
                for (Glycan dbGlycan : params.glycanDatabase) {
                    String dbGlycanStr = dbGlycan.name;
                    if (site.observedGlycans.containsKey(dbGlycanStr)) {
                        // this glycan was observed for this site, no action required
                        continue;
                    }
                    // glycan not observed. Make a PSM for it
                    newPSMs.add(PSM.makeDummyPSM(site.basePSM, dbGlycan, psmFile.PSMcolumns, params, rankCounter));
                    rankCounter++;
                }
            }
        }

        return newPSMs;
    }


}
