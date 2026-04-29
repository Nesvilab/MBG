package com.mbg;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Accumulates per-candidate diagnostic records and writes mbg_diagnostics.tsv.
 *
 * Called from GlycanMatching.fillInGlycoSites after scoring each candidate.
 * Decoy records and posterior/q-value columns are left blank until the
 * four-class scoring refactor in Commits 5–6.
 */
public class Diagnostics {

    /** One row per candidate in the diagnostics output. */
    public static class CandidateRecord {
        String rawFile;
        String peptide;
        String glycanComposition;
        double precursorMz;
        int charge;
        /** type1_target | type2_target (decoy variants added in Commit 5) */
        String candidateClass;
        float rtShift;
        Float imShift;   // null when --nopasef true
        double massErrorPpm;
        double klDivergence;
        double precursorIntensity;
        double shiftFrequency;
        double featureIntensityY0Y1;
        double ldaScore;
        /** blank until Commit 6 EM step */
        String posteriorProbability = "";
        /** blank until Commit 6 q-value step */
        String qValue = "";
    }

    /** Accumulated records across all raw files processed so far. */
    public static List<CandidateRecord> records = new ArrayList<>();

    /** Reset between experiments (called alongside Scoring.init()). */
    public static void reset() {
        records = new ArrayList<>();
    }

    /**
     * Add a record for a scored target candidate.
     *
     * @param psm             the inferred PSM
     * @param ldaScore        the LDA score from Scoring.getMBGScore
     * @param posteriorProb   posterior probability from mixture-model EM; NaN for TYPE_2
     * @param qValue          q-value from sorted posteriors; NaN for TYPE_2
     * @param noPASEF         true if IM data is absent (leaves im_shift blank)
     */
    public static void addTargetRecord(PSM psm, double ldaScore, double posteriorProb, double qValue, boolean noPASEF) {
        CandidateRecord r = new CandidateRecord();
        r.rawFile = psm.rawfile;
        r.peptide = psm.peptide;
        r.glycanComposition = psm.glycan != null ? psm.glycan.toString() : "";
        r.precursorMz = psm.calcMz;
        r.charge = psm.charge;
        r.candidateClass = (psm.candidateType == PSM.CandidateType.TYPE_2)
                ? "type2_target" : "type1_target";
        r.rtShift = psm.rtDiff;
        r.imShift = noPASEF ? null : psm.imDiff;
        r.massErrorPpm = psm.avePPM;
        r.klDivergence = psm.kl;
        r.precursorIntensity = psm.intensity;
        r.shiftFrequency = psm.shiftFrequency;
        r.featureIntensityY0Y1 = psm.featureIntensity;
        r.ldaScore = ldaScore;
        r.posteriorProbability = Double.isNaN(posteriorProb) ? "" : String.format("%.6f", posteriorProb);
        r.qValue = Double.isNaN(qValue) ? "" : String.format("%.6f", qValue);
        records.add(r);
    }

    private static final String TSV_HEADER =
            "raw_file\tpeptide\tglycan_composition\tprecursor_mz\tcharge\t" +
            "class\trt_shift\tim_shift\tmass_error_ppm\tkl_divergence\t" +
            "precursor_intensity\tshift_frequency\tfeature_intensity_y0_y1\t" +
            "lda_score\tposterior_probability\tq_value";

    /**
     * Writes mbg_diagnostics.tsv next to the psm.tsv output file.
     *
     * @param psmFile  the psm.tsv File object (output directory is its parent)
     * @param noPASEF  true if IM column should be blank
     */
    public static void write(File psmFile, boolean noPASEF) {
        File outputFile = new File(psmFile.getParent(), "mbg_diagnostics.tsv");
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            writer.write(TSV_HEADER);
            writer.newLine();
            for (CandidateRecord r : records) {
                StringBuilder sb = new StringBuilder();
                sb.append(r.rawFile).append('\t');
                sb.append(r.peptide).append('\t');
                sb.append(r.glycanComposition).append('\t');
                sb.append(String.format("%.4f", r.precursorMz)).append('\t');
                sb.append(r.charge).append('\t');
                sb.append(r.candidateClass).append('\t');
                sb.append(String.format("%.4f", r.rtShift)).append('\t');
                sb.append(r.imShift != null ? String.format("%.4f", r.imShift) : "").append('\t');
                sb.append(String.format("%.4f", r.massErrorPpm)).append('\t');
                sb.append(String.format("%.6f", r.klDivergence)).append('\t');
                sb.append(String.format("%.4f", r.precursorIntensity)).append('\t');
                sb.append(String.format("%.6f", r.shiftFrequency)).append('\t');
                sb.append(String.format("%.6f", r.featureIntensityY0Y1)).append('\t');
                sb.append(String.format("%.6f", r.ldaScore)).append('\t');
                sb.append(r.posteriorProbability).append('\t');
                sb.append(r.qValue);
                writer.write(sb.toString());
                writer.newLine();
            }
            System.out.printf("Wrote %d candidate records to %s%n", records.size(), outputFile.getPath());
        } catch (IOException e) {
            System.err.printf("Warning: could not write mbg_diagnostics.tsv: %s%n", e.getMessage());
        }
    }
}
