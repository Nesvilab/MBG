package com.mbg;

import java.util.HashMap;

/**
 * AA_masses provides monoisotopic masses for amino acid residues used in peptide mass calculations.
 * 
 * The class contains a static map of single-letter amino acid codes to their corresponding
 * monoisotopic residue masses (in Daltons). These masses represent the mass of each amino acid
 * residue within a peptide sequence (without the mass of H2O that would be added for a free amino acid).
 */
public class AA_masses {

    /**
     * Map of amino acid single-letter codes to their monoisotopic residue masses in Daltons.
     * Standard amino acids have their accurate masses, while non-standard codes (J, U, X, Z) 
     * are assigned zero mass.
     */
    public static final HashMap<Character, Double> aa_masses;
    
    static {
        aa_masses = new HashMap<>();
        
        // Standard amino acids
        aa_masses.put('A', 71.03711);  // Alanine
        aa_masses.put('C', 103.00919); // Cysteine
        aa_masses.put('D', 115.02694); // Aspartic acid
        aa_masses.put('E', 129.04259); // Glutamic acid
        aa_masses.put('F', 147.06841); // Phenylalanine
        aa_masses.put('G', 57.02146);  // Glycine
        aa_masses.put('H', 137.05891); // Histidine
        aa_masses.put('I', 113.08406); // Isoleucine
        aa_masses.put('K', 128.09496); // Lysine
        aa_masses.put('L', 113.08406); // Leucine
        aa_masses.put('M', 131.04049); // Methionine
        aa_masses.put('N', 114.04293); // Asparagine
        aa_masses.put('P', 97.05276);  // Proline
        aa_masses.put('Q', 128.05858); // Glutamine
        aa_masses.put('R', 156.10111); // Arginine
        aa_masses.put('S', 87.03203);  // Serine
        aa_masses.put('T', 101.04768); // Threonine
        aa_masses.put('V', 99.06841);  // Valine
        aa_masses.put('W', 186.07931); // Tryptophan
        aa_masses.put('Y', 163.06333); // Tyrosine
        
        // Special cases
        aa_masses.put('O', 114.07931); // Pyrrolysine
        
        // Non-standard codes (assigned zero mass)
        aa_masses.put('J', 0.0); // Not a standard amino acid (sometimes used as I/L ambiguity)
        aa_masses.put('U', 0.0); // Selenocysteine (not commonly used in standard analyses)
        aa_masses.put('X', 0.0); // Unknown amino acid
        aa_masses.put('Z', 0.0); // Not a standard amino acid (sometimes used as E/Q ambiguity)
    }
    
    /**
     * Calculates the total mass of a peptide sequence from its constituent amino acids.
     * 
     * @param sequence The peptide sequence as a string of single-letter amino acid codes
     * @return The calculated mass of the peptide sequence in Daltons
     */
    public static double calculatePeptideMass(String sequence) {
        // Start with mass of H2O to account for the extra H on N-terminus and OH on C-terminus
        double mass = 18.01056;
        
        // Add mass of each amino acid residue
        for (char aminoAcid : sequence.toCharArray()) {
            Double residueMass = aa_masses.get(aminoAcid);
            if (residueMass != null) {
                mass += residueMass;
            } else {
                // Handle unknown amino acid codes by using a default mass or logging
                System.err.println("Unknown amino acid code: " + aminoAcid);
            }
        }
        
        return mass;
    }
}
