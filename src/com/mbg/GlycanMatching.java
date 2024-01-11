package com.mbg;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class GlycanMatching {
    public PSM_File psmFile;
    private Params params;
    private HashMap<String, GlycoSite> observedGlycoSites;

    public GlycanMatching(Params params) throws Exception {
        psmFile = new PSM_File(params);
        this.params = params;
    }

    public void groupGlycosites() {
        observedGlycoSites = new HashMap<>();

        for (PSM psm: psmFile.psms) {
            if (psm.glycan == null || psm.inferredGlycoform) {
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
                    String dbGlycanStr = dbGlycan.toString();
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
