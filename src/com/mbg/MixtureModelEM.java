package com.mbg;

import java.util.List;

/**
 * Two-component 1-D Gaussian mixture model fitted by EM for the four-class FDR
 * framework (Commit 6).  Fitted on the LDA scores of TYPE_1 targets and decoys;
 * the "true" and "false" components are initialized from the respective score
 * statistics so the EM converges to the right labelling.
 */
public class MixtureModelEM {

    /** Minimum combined (target + decoy) sample size required to run EM. */
    public static final int N_MIN = 50;

    private static final int MAX_ITER = 200;
    private static final double CONV_TOL = 1e-6;
    private static final double VAR_FLOOR = 1e-6;

    private double muTrue, muFalse;
    private double sigmaTrue, sigmaFalse;
    private double pi;  // mixing weight for "true" component

    private boolean converged = false;

    /**
     * Initialize parameters from the score statistics of TYPE_1 targets and decoys.
     * The component with the higher mean is labelled "true".
     */
    public MixtureModelEM(List<Double> targetScores, List<Double> decoyScores) {
        double mt = mean(targetScores);
        double md = mean(decoyScores);
        double st = std(targetScores, mt);
        double sd = std(decoyScores, md);
        int nt = targetScores.size();
        int nd = decoyScores.size();

        // Label the component with the larger mean as "true"
        if (mt >= md) {
            muTrue = mt; sigmaTrue = st;
            muFalse = md; sigmaFalse = sd;
            pi = (double) nt / (nt + nd);
        } else {
            muTrue = md; sigmaTrue = sd;
            muFalse = mt; sigmaFalse = st;
            pi = (double) nd / (nt + nd);
        }
    }

    /**
     * Fit the mixture model by EM on the combined score list (targets + decoys, unlabeled).
     * Populates muTrue, muFalse, sigmaTrue, sigmaFalse, and pi.
     */
    public void fit(List<Double> scores) {
        double[] data = scores.stream().mapToDouble(Double::doubleValue).toArray();
        double prevLogLik = Double.NEGATIVE_INFINITY;

        for (int iter = 0; iter < MAX_ITER; iter++) {
            // E-step: compute responsibility of "true" component for each datum
            double[] r = new double[data.length];
            double logLik = 0.0;
            for (int i = 0; i < data.length; i++) {
                double pt = pi * gaussian(data[i], muTrue, sigmaTrue);
                double pf = (1.0 - pi) * gaussian(data[i], muFalse, sigmaFalse);
                double total = pt + pf;
                r[i] = (total > 0) ? pt / total : 0.5;
                logLik += Math.log(total + 1e-300);
            }

            // M-step: update parameters from responsibilities
            double sumR = 0, sumRx = 0, sumRx2 = 0;
            double sum1R = 0, sum1Rx = 0, sum1Rx2 = 0;
            for (int i = 0; i < data.length; i++) {
                sumR  += r[i];
                sumRx  += r[i] * data[i];
                sumRx2 += r[i] * data[i] * data[i];
                sum1R  += (1.0 - r[i]);
                sum1Rx  += (1.0 - r[i]) * data[i];
                sum1Rx2 += (1.0 - r[i]) * data[i] * data[i];
            }

            pi = sumR / data.length;
            if (sumR > 0) {
                muTrue = sumRx / sumR;
                double varTrue = sumRx2 / sumR - muTrue * muTrue;
                sigmaTrue = Math.sqrt(Math.max(varTrue, VAR_FLOOR));
            }
            if (sum1R > 0) {
                muFalse = sum1Rx / sum1R;
                double varFalse = sum1Rx2 / sum1R - muFalse * muFalse;
                sigmaFalse = Math.sqrt(Math.max(varFalse, VAR_FLOOR));
            }

            // Enforce muTrue >= muFalse (relabel if converged to swapped state)
            if (muTrue < muFalse) {
                double tmp;
                tmp = muTrue;    muTrue    = muFalse;    muFalse    = tmp;
                tmp = sigmaTrue; sigmaTrue = sigmaFalse; sigmaFalse = tmp;
                pi = 1.0 - pi;
            }

            // Convergence check
            double relDelta = Math.abs(logLik - prevLogLik) / (Math.abs(prevLogLik) + 1e-300);
            if (iter > 0 && relDelta < CONV_TOL) {
                converged = true;
                break;
            }
            prevLogLik = logLik;
        }
        if (!converged) {
            System.out.println("MixtureModelEM: did not converge within " + MAX_ITER + " iterations.");
        }
    }

    /**
     * Returns the posterior probability that the given LDA score belongs to the "true" component.
     */
    public double posterior(double score) {
        double pt = pi * gaussian(score, muTrue, sigmaTrue);
        double pf = (1.0 - pi) * gaussian(score, muFalse, sigmaFalse);
        double total = pt + pf;
        return (total > 0) ? pt / total : 0.0;
    }

    private static double gaussian(double x, double mu, double sigma) {
        double z = (x - mu) / sigma;
        return Math.exp(-0.5 * z * z) / (sigma * Math.sqrt(2.0 * Math.PI));
    }

    private static double mean(List<Double> values) {
        if (values.isEmpty()) return 0.0;
        return values.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
    }

    private static double std(List<Double> values, double m) {
        if (values.size() < 2) return 1.0;
        double var = values.stream().mapToDouble(v -> (v - m) * (v - m)).average().orElse(1.0);
        return Math.sqrt(Math.max(var, VAR_FLOOR));
    }
}
