package com.mbg;

import java.util.*;

import org.apache.commons.math3.linear.*;

/**
 * The Scoring class implements Linear Discriminant Analysis (LDA) to score and 
 * filter inferred glycopeptides based on multiple features.
 * 
 * This class calculates discriminant functions to separate target PSMs from decoy PSMs,
 * and applies an FDR-based threshold to filter glycopeptides.
 */
public class Scoring {
    /** Collection of feature vectors for decoy PSMs */
    public static List<double[]> decoyData;
    
    /** Collection of feature vectors for target PSMs */
    public static List<double[]> targetData;
    
    /** Score threshold for accepting glycopeptides */
    public static double thresholdScore;
    
    /** Coefficients of the LDA discriminant function */
    private static double[] coefficients;

    /** Intercept term of the LDA discriminant function */
    private static double intercept;

    /** Per-feature means used for z-score standardization, fit on combined target+decoy training data (per raw file). */
    private static double[] featureMeans;

    /** Per-feature standard deviations used for z-score standardization (per raw file). */
    private static double[] featureStds;

    /** Minimum std to avoid divide-by-zero; features with std below this are treated as inert (scaled value = 0). */
    private static final double STD_EPS = 1e-9;

    /** Index of the intensity feature in the per-PSM feature vector. */
    private static final int INTENSITY_IDX = 0;

    /** Natural-log denominator used to compute log base 2 of the intensity feature. */
    private static final double LOG2 = Math.log(2.0);

    /**
     * Initializes the scoring system by creating new data collections.
     * Must be called before any scoring operations are performed.
     */
    public static void init() {
        decoyData = new ArrayList<>();
        targetData = new ArrayList<>();
        featureMeans = null;
        featureStds = null;
    }

    /**
     * Calculates the LDA model and determines the score threshold based on the 
     * specified false discovery rate (FDR).
     * 
     * @param fdrCutoff The maximum acceptable FDR
     * @throws Exception If matrix operations fail
     */
    public static void getCutOffScore(double fdrCutoff) throws Exception {
        // Z-score standardize features using combined target+decoy training data (per raw file).
        // This brings features onto comparable scales so the pooled covariance matrix is well
        // conditioned and the printed LDA coefficients reflect relative feature importance.
        standardizeData();

        // Calculate mean vectors for target and decoy datasets
        double[] decoyMean = calculateMean(decoyData);
        double[] targetMean = calculateMean(targetData);

        // Calculate covariance matrices
        RealMatrix decoyCov = calculateCovariance(decoyData, decoyMean);
        RealMatrix targetCov = calculateCovariance(targetData, targetMean);

        // Calculate pooled covariance matrix
        RealMatrix pooledCov = calculatePooledCovariance(targetCov, decoyCov, 
                                                        targetData.size(), decoyData.size());

        // Check for singularity and regularize if necessary
        double det = new LUDecomposition(pooledCov).getDeterminant();
        if (Double.isNaN(det) || det == 0) {
            System.out.println("Warning: Pooled covariance matrix is singular or invalid. Regularizing...");
            pooledCov = regularizeCovarianceMatrix(pooledCov, 0.001);
        }
//        if (new LUDecomposition(pooledCov).getDeterminant() == 0) {
//            pooledCov = regularizeCovarianceMatrix(pooledCov, 0.001);
//        }
        
        // Invert the covariance matrix
        LUDecomposition lu = new LUDecomposition(pooledCov);
        RealMatrix invPooledCov = lu.getSolver().getInverse();

        // Calculate LDA model parameters
        coefficients = calculateCoefficients(targetMean, decoyMean, invPooledCov);
        intercept = calculateIntercept(targetMean, decoyMean, coefficients);

        // Output model coefficients
        System.out.println("Feature means (z-score): " + java.util.Arrays.toString(featureMeans));
        System.out.println("Feature stds (z-score): " + java.util.Arrays.toString(featureStds));
        System.out.println("Coefficients: " + java.util.Arrays.toString(coefficients));

        // Calculate scores and determine threshold
        List<Double> targetScores = calculateScores(targetData, coefficients, intercept);
        List<Double> decoyScores = calculateScores(decoyData, coefficients, intercept);
        thresholdScore = getCutOff(fdrCutoff, targetScores, decoyScores);
        System.out.println("Threshold Score: " + thresholdScore);
    }

    /**
     * Checks if a score exceeds the calculated threshold.
     * 
     * @param score The score to evaluate
     * @return true if the score is above the threshold
     */
    public static boolean overThreshold(double score) {
        return score > thresholdScore;
    }

    /**
     * Calculates the MBG score for a PSM and updates the PSM's score.
     * 
     * @param psm The PSM to score
     * @param psmCols Column mapping for PSM fields
     * @param noPASEF Whether to exclude ion mobility features
     * @return The calculated MBG score
     */
    public static double getMBGScore(PSM psm, HashMap<PSM_File.PSMattr, Integer> psmCols, boolean noPASEF) {
        // Create feature vector with or without ion mobility data
        double[] featureVector;
        if (noPASEF) {
            featureVector = new double[]{
                psm.intensity, 
                psm.avePPM, 
                psm.kl, 
                psm.shiftFrequency, 
                psm.featureIntensity, 
                (double) psm.rtDiff
            };
        } else {
            featureVector = new double[]{
                psm.intensity, 
                psm.avePPM, 
                psm.kl, 
                psm.shiftFrequency, 
                psm.featureIntensity, 
                (double) psm.rtDiff, 
                (double) psm.imDiff
            };
        }
        
        // Apply the same z-score standardization used during LDA fitting
        double[] scaled = standardize(featureVector);

        // Calculate score using LDA model
        double score = calculateScore(scaled, coefficients, intercept);

        // Update PSM with calculated score
        psm.add_MBG_score(psmCols, score);
        return score;
    }

    /**
     * Computes per-feature mean and std from the combined target+decoy training set
     * and rescales both lists in place to z-scores. Features with near-zero variance
     * are collapsed to 0 to avoid numerical instability.
     */
    private static void standardizeData() {
        if (targetData.isEmpty() && decoyData.isEmpty()) {
            return;
        }
        int nFeatures = targetData.isEmpty() ? decoyData.get(0).length : targetData.get(0).length;
        long total = (long) targetData.size() + (long) decoyData.size();

        // Compress intensity's 6+ order-of-magnitude dynamic range before z-scoring so the
        // resulting distribution is closer to Gaussian and less outlier-dominated.
        for (double[] v : targetData) logTransformInPlace(v);
        for (double[] v : decoyData) logTransformInPlace(v);

        double[] means = new double[nFeatures];
        for (double[] v : targetData) {
            for (int i = 0; i < nFeatures; i++) means[i] += v[i];
        }
        for (double[] v : decoyData) {
            for (int i = 0; i < nFeatures; i++) means[i] += v[i];
        }
        for (int i = 0; i < nFeatures; i++) means[i] /= total;

        double[] vars = new double[nFeatures];
        for (double[] v : targetData) {
            for (int i = 0; i < nFeatures; i++) {
                double d = v[i] - means[i];
                vars[i] += d * d;
            }
        }
        for (double[] v : decoyData) {
            for (int i = 0; i < nFeatures; i++) {
                double d = v[i] - means[i];
                vars[i] += d * d;
            }
        }
        double[] stds = new double[nFeatures];
        long dof = Math.max(total - 1, 1);
        for (int i = 0; i < nFeatures; i++) {
            stds[i] = Math.sqrt(vars[i] / dof);
        }

        featureMeans = means;
        featureStds = stds;

        for (double[] v : targetData) scaleInPlace(v);
        for (double[] v : decoyData) scaleInPlace(v);
    }

    /**
     * Returns a z-scored copy of a feature vector using the means/stds captured during training.
     * Falls through unchanged if standardization was not initialized (e.g. scorer bypassed).
     */
    private static double[] standardize(double[] vec) {
        if (featureMeans == null || featureStds == null) {
            return vec;
        }
        double[] out = vec.clone();
        logTransformInPlace(out);
        for (int i = 0; i < out.length; i++) {
            out[i] = featureStds[i] > STD_EPS ? (out[i] - featureMeans[i]) / featureStds[i] : 0.0;
        }
        return out;
    }

    /**
     * Applies a log2(x + 1) transform to the intensity slot (index 0) of the given vector
     * in place. Intensity spans roughly 10^3-10^11 in practice; compressing it before
     * z-scoring avoids the tail dominating the mean/std and keeps the covariance matrix
     * well behaved.
     */
    private static void logTransformInPlace(double[] vec) {
        vec[INTENSITY_IDX] = Math.log(vec[INTENSITY_IDX] + 1.0) / LOG2;
    }

    /** In-place variant of {@link #standardize(double[])} used during training-set rescaling. */
    private static void scaleInPlace(double[] vec) {
        for (int i = 0; i < vec.length; i++) {
            vec[i] = featureStds[i] > STD_EPS ? (vec[i] - featureMeans[i]) / featureStds[i] : 0.0;
        }
    }

    //--------------------------------
    // SCORING CALCULATION METHODS
    //--------------------------------

    /**
     * Calculates the LDA score for a data point.
     * 
     * @param dataPoint Feature vector to score
     * @param coefficients LDA model coefficients
     * @param intercept LDA model intercept
     * @return The calculated score
     */
    private static double calculateScore(double[] dataPoint, double[] coefficients, double intercept) {
        return dotProduct(dataPoint, coefficients) + intercept;
    }

    /**
     * Calculates scores for a collection of data points.
     * 
     * @param data Collection of feature vectors
     * @param coefficients LDA model coefficients
     * @param intercept LDA model intercept
     * @return List of calculated scores
     */
    public static List<Double> calculateScores(List<double[]> data, double[] coefficients, double intercept) {
        List<Double> scores = new ArrayList<>();
        for (double[] dataPoint : data) {
            double score = calculateScore(dataPoint, coefficients, intercept);
            scores.add(score);
        }
        return scores;
    }

    /**
     * Determines the score threshold based on the specified FDR.
     * 
     * @param fdrCutOff Maximum acceptable FDR
     * @param targetScores Scores for target PSMs
     * @param decoyScores Scores for decoy PSMs
     * @return The score threshold
     */
    public static double getCutOff(double fdrCutOff, List<Double> targetScores, List<Double> decoyScores) {
        // Sort target scores in descending order
        targetScores.sort(Comparator.reverseOrder());
        
        // Get decoy indices in the combined sorted list
        List<Integer> decoyIndexes = getDecoyIndexes(targetScores, decoyScores);
        
        // Calculate FDRs at each point
        List<Double> fdrs = new ArrayList<>();
        for (int i = 0; i < decoyIndexes.size(); i++) {
            double targetCount = decoyIndexes.get(i) - i;   // Number of total PSMs before this decoy minus number of decoys
            double decoyCount = i + 1;      // Number of decoys up to this point
            double fdr = (decoyCount + 1) / targetCount;        // q = (d+1)/t recommended per 10.1021/acs.jproteome.6b00144
            fdrs.add(fdr);
        }
        
        // Find the last index where FDR is below threshold
        int decoyIndex = getLastValueBelowThreshold(fdrs, fdrCutOff);
        int targetIndex = decoyIndexes.get(decoyIndex);
        
        // Handle edge case
        if (targetIndex - decoyIndex <= 0) {
            return targetScores.get(0);
        }
        
        return targetScores.get(targetIndex - decoyIndex - 1);
    }

    //--------------------------------
    // HELPER METHODS
    //--------------------------------

    /**
     * Gets the indices of decoy scores in the combined sorted score list.
     */
    private static List<Integer> getDecoyIndexes(List<Double> targetScores, List<Double> decoyScores) {
        // Combine scores
        List<Double> combinedScores = new ArrayList<>(targetScores);
        combinedScores.addAll(decoyScores);

        // Track original indices
        List<Integer> originalIndexes = new ArrayList<>();
        for (int i = 0; i < combinedScores.size(); i++) {
            originalIndexes.add(i);
        }

        // Sort by score in descending order
        originalIndexes.sort(Comparator.comparingDouble(combinedScores::get).reversed());

        // Extract decoy indices
        List<Integer> decoyIndexes = new ArrayList<>();
        for (int i = 0; i < originalIndexes.size(); i++) {
            if (originalIndexes.get(i) >= targetScores.size()) {
                decoyIndexes.add(i);
            }
        }

        return decoyIndexes;
    }

    /**
     * Finds the index of the last value in a list that is below a threshold.
     */
    public static Integer getLastValueBelowThreshold(List<Double> list, double threshold) {
        for (int i = list.size() - 1; i >= 0; i--) {
            if (list.get(i) <= threshold) {
                return i;
            }
        }
        return 0; // Default if no value found
    }

    /**
     * Calculates the mean vector for a dataset.
     */
    private static double[] calculateMean(List<double[]> data) {
        int featureCount = data.get(0).length;
        double[] mean = new double[featureCount];
        
        // Sum values for each feature
        for (double[] dataPoint : data) {
            for (int i = 0; i < featureCount; i++) {
                mean[i] += dataPoint[i];
            }
        }
        
        // Divide by count to get mean
        for (int i = 0; i < featureCount; i++) {
            mean[i] /= data.size();
        }
        
        return mean;
    }

    /**
     * Calculates the covariance matrix for a dataset.
     */
    private static RealMatrix calculateCovariance(List<double[]> data, double[] mean) {
        int featureCount = data.get(0).length;
        RealMatrix covariance = new Array2DRowRealMatrix(featureCount, featureCount);
        
        // Calculate covariance sums
        for (double[] dataPoint : data) {
            for (int i = 0; i < featureCount; i++) {
                for (int j = 0; j < featureCount; j++) {
                    double currentValue = covariance.getEntry(i, j);
                    double update = (dataPoint[i] - mean[i]) * (dataPoint[j] - mean[j]);
                    covariance.setEntry(i, j, currentValue + update);
                }
            }
        }
        
        // Normalize by degrees of freedom
        for (int i = 0; i < featureCount; i++) {
            for (int j = 0; j < featureCount; j++) {
                double normalizedValue = covariance.getEntry(i, j) / (data.size() - 1);
                covariance.setEntry(i, j, normalizedValue);
            }
        }
        
        return covariance;
    }

    /**
     * Calculates the pooled covariance matrix from two class-specific matrices.
     */
    private static RealMatrix calculatePooledCovariance(RealMatrix cov1, RealMatrix cov2, int n1, int n2) {
        double weight1 = (double) (n1 - 1) / (n1 + n2 - 2);
        double weight2 = (double) (n2 - 1) / (n1 + n2 - 2);
        return cov1.scalarMultiply(weight1).add(cov2.scalarMultiply(weight2));
    }

    /**
     * Regularizes a covariance matrix to ensure it's invertible.
     */
    private static RealMatrix regularizeCovarianceMatrix(RealMatrix covarianceMatrix, double alpha) {
        int size = covarianceMatrix.getRowDimension();
        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(size);
        return covarianceMatrix.scalarMultiply(1 - alpha).add(identity.scalarMultiply(alpha));
    }

    /**
     * Calculates the coefficients for the LDA discriminant function.
     */
    private static double[] calculateCoefficients(double[] mean1, double[] mean2, RealMatrix invPooledCov) {
        double[] coefficients = new double[mean1.length];
        double[] difference = new double[mean1.length]; 
        
        // Calculate mean difference vector
        for (int i = 0; i < mean1.length; i++) {
            difference[i] = mean1[i] - mean2[i]; 
        }
        
        // Calculate coefficients using matrix operations
        ArrayRealVector diffRealVector = new ArrayRealVector(difference);
        for (int i = 0; i < mean1.length; i++) {
            coefficients[i] = invPooledCov.getRowVector(i).dotProduct(diffRealVector);
        }
        
        return coefficients;
    }

    /**
     * Calculates the intercept term for the LDA discriminant function.
     */
    private static double calculateIntercept(double[] mean1, double[] mean2, double[] coefficients) {
        return -0.5 * (dotProduct(mean1, coefficients) + dotProduct(mean2, coefficients));
    }

    /**
     * Calculates the dot product of two vectors.
     */
    private static double dotProduct(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("Arrays must have the same length");
        }
        
        double result = 0;
        for (int i = 0; i < a.length; i++) {
            result += a[i] * b[i];
        }
        
        return result;
    }
}


