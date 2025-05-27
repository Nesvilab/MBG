package com.mbg;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * GlycoSite represents a glycosylation site on a peptide, tracking all observed
 * glycans at this site along with their corresponding PSMs (Peptide-Spectrum Matches).
 * 
 * This class is used to group PSMs based on peptide sequence and manage collections
 * of glycans observed at each site during the glycoform inference process.
 */
public class GlycoSite {

    /** All PSMs associated with this glycosite */
    ArrayList<PSM> psms;
    
    /** Map of glycan compositions to their corresponding PSMs */
    HashMap<String, ArrayList<PSM>> observedGlycans;
    
    /** Protein identifier */
    String protein;
    
    /** Peptide sequence identifier for this glycosite */
    String id;
    
    /** Reference to the primary PSM used as a template for inference */
    PSM basePSM;
    
    /** Collection of the most intense PSMs for each glycopeptide at this site */
    ArrayList<PSM> basePSMs;
    
    /** List of observed glycan masses at this site (used for mass filtering) */
    ArrayList<Double> glycanMasses;

    /**
     * Creates a new GlycoSite instance with the first observed PSM.
     * 
     * @param firstPSM The first PSM observed at this site
     * @param peptideSeq The peptide sequence defining this glycosite
     */
    public GlycoSite(PSM firstPSM, String peptideSeq) {
        this.protein = firstPSM.protein;
        this.id = peptideSeq;
        
        // Initialize collections
        psms = new ArrayList<>();
        psms.add(firstPSM);
        basePSM = firstPSM;
        
        // Initialize glycan tracking
        observedGlycans = new HashMap<>();
        ArrayList<PSM> glycoPSMs = new ArrayList<>();
        glycoPSMs.add(firstPSM);
        observedGlycans.put(firstPSM.glycan.name, glycoPSMs);
        
        basePSMs = new ArrayList<>();
        glycanMasses = new ArrayList<>();
        glycanMasses.add(firstPSM.glycan.mass);
    }

    /**
     * Adds a PSM to this glycosite and updates the observed glycans collection.
     * If the glycan has not been observed before, it's added to the collection.
     * 
     * @param psm The PSM to add to this glycosite
     */
    public void addPSM(PSM psm) {
        // Add to master list of PSMs
        psms.add(psm);
        
        // Get glycan composition name
        String glycanComposition = psm.glycan.name;
        
        // Get or create list of PSMs for this glycan
        ArrayList<PSM> glycoPSMs;
        if (observedGlycans.containsKey(glycanComposition)) {
            glycoPSMs = observedGlycans.get(glycanComposition);
        } else {
            glycoPSMs = new ArrayList<>();
            glycanMasses.add(psm.glycan.mass);
        }
        
        // Add PSM to the glycan-specific list
        glycoPSMs.add(psm);
        observedGlycans.put(glycanComposition, glycoPSMs);
    }
}
