package com.mbg;

import java.util.*;

import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanMod;

import static java.lang.System.out;

public class GlycanMatching {
    public PSM_File psmFile;
    public ArrayList<PSM> psms;
    private final Params params;
    private HashMap<String, GlycoSite> observedGlycoSites;
    private String[] qualifiedGlycans;  // glycans in psm.tsv with a specific minimum count

    public GlycanMatching(Params params, PSM_File psmFile, ArrayList<PSM> psms) throws Exception {
        this.psmFile = psmFile;
        this.psms = psms;
        this.params = params;
    }

    public void groupGlycosites() {
        observedGlycoSites = new HashMap<>();
        HashMap<String, Integer> allObservedGlycans = new HashMap<>();

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
            // counts observed glycans
            if (allObservedGlycans.containsKey(psm.glycan.toString())) {
                allObservedGlycans.put(psm.glycan.toString(), allObservedGlycans.get(psm.glycan.toString()) + 1);
            } else {
                allObservedGlycans.put(psm.glycan.toString(), 1);
            }
        }
        // filter to keep entries where the glycan psm counts >5
        qualifiedGlycans = allObservedGlycans.entrySet()
                .stream()
                .filter(entry -> entry.getValue() > 5)
                .map(Map.Entry::getKey).toArray(String[]::new);
    }

    /**
     * Given a glycan database and rules from Params, fill in the unobserved glycans for each glycosite
     * and generate dummy PSMs for them.
     * @return
     */
    public ArrayList<PSM> fillInGlycoSites(String rawFileName, String msFilePath) throws Exception {
        ArrayList<PSM> newPSMs = new ArrayList<>();
        GlycanMod mod_nh3 = new GlycanMod(params.glycanResidues.get("NH4"));

        // read msFile
        MSMS ms2 = new MSMS(rawFileName, msFilePath);
        ms2.msFileReader();
        for (Map.Entry<String, GlycoSite> siteEntry: observedGlycoSites.entrySet()) {
            GlycoSite site = siteEntry.getValue();
            if (site.psms.size() >= params.minPSMs && site.observedGlycans.size() >= params.minGlycans) {
                for (PSM base_PSM : site.psms) {
                    int rankCounter = 1;
                    ArrayList<Glycan> dummy_glycans = PSM.getDummyGlycans(base_PSM, params.toAddResidues);
                    for (Glycan new_glycan: dummy_glycans) {
                        if (site.observedGlycans.containsKey(new_glycan.toString()) || site.basePSM == null) {
                            // this glycan was observed for this site, no action required
                            continue;
                        }
                        if (!Arrays.asList(qualifiedGlycans).contains(new_glycan.toString())) {
                            // this glycan not detected in psm.tsv
                            continue;
                        }

                        // glycan not observed. Make a PSM for it
                        if (base_PSM.scanNum == 10876 || base_PSM.scanNum == 10907) {
                            out.print("sadsad");
                        }
                        List<PSM> addedPSMs = PSM.makeDummyPSM(base_PSM, new_glycan, psmFile.PSMcolumns, params, rankCounter);
                        newPSMs.addAll(addedPSMs);
                        // update observed glycans
                        site.observedGlycans.put(new_glycan.toString(), new ArrayList<>());
                        rankCounter++;
                    }

                }
                // make NH3 glycan
                int rankCounter = 1;
                for (PSM base_psm : site.psms) {
                    Glycan new_glycan = PSM.getNH3PSM(base_psm.glycan, mod_nh3);
                    List<PSM> addedPSMs = PSM.makeDummyPSM(base_psm, new_glycan, psmFile.PSMcolumns, params, rankCounter);
                    newPSMs.addAll(addedPSMs);
                }
            }
        }
        return newPSMs;
    }
}
