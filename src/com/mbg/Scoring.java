package com.mbg;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.ArrayRealVector;


public class Scoring {
    public static List<double[]> decoyData;
    public static List<double[]> targetData;
    public static double thresholdScore;
    private static double[] coefficients;
    private static double intercept;


    public static void init() {
        decoyData = new ArrayList<>();
        targetData = new ArrayList<>();
    }

    public static void getCutOffScore(double fdrCurOff) throws Exception {
        // Calculate mean vectors
        double[] decoyMean = calculateMean(decoyData);
        double[] targetMean = calculateMean(targetData);

        // Calculate covariance matrices
        RealMatrix decoyCov = calculateCovariance(decoyData, decoyMean);
        RealMatrix targetCov = calculateCovariance(targetData, targetMean);

        // Calculate pooled covariance matrix
        RealMatrix pooledCov = calculatePooledCovariance(targetCov, decoyCov, decoyData.size(), targetData.size());

        // Check for singularity and handle it
        LUDecomposition lu = new LUDecomposition(pooledCov);
        if (lu.getDeterminant() == 0) {
            System.err.println("Pooled covariance matrix is singular. Cannot calculate coefficients.");
            // Handle the error appropriately (e.g., exit the program, use a different method)
            return;
        }

        // Calculate inverse of pooled covariance matrix
        RealMatrix invPooledCov = lu.getSolver().getInverse();

        // Calculate coefficients of the discriminant function
        coefficients = calculateCoefficients(targetMean, decoyMean, invPooledCov);

        // Calculate intercept of the discriminant function
        intercept = calculateIntercept(targetMean, decoyMean, coefficients);

        // Print results
        System.out.println("Coefficients: " + java.util.Arrays.toString(coefficients));
        System.out.println("Intercept: " + intercept);

        // calculate scores for targets and decoys
        List<Double> targetScores = calculateScores(targetData, coefficients, intercept);
        List<Double> decoyScores = calculateScores(decoyData, coefficients, intercept);
        thresholdScore = getCutOff(fdrCurOff, targetScores, decoyScores);
        // Classify a new data point (replace with your actual data)
//        double[] newDataPoint = new double[]{3.5, 4.5, 5.5};
//        double score = calculateScore(newDataPoint, coefficients, intercept);
//        System.out.println("Score: " + score);
//        System.out.println("Classification: " + (score > 0 ? "Target" : "Decoy"));
    }

    // Helper functions
    private static double[] calculateMean(List<double[]> data) {
        double[] mean = new double[data.get(0).length];
        for (double[] row : data) {
            for (int i = 0; i < row.length; i++) {
                mean[i] += row[i];
            }
        }
        for (int i = 0; i < mean.length; i++) {
            mean[i] /= data.size();
        }
        return mean;
    }

    private static RealMatrix calculateCovariance(List<double[]> data, double[] mean) {
        int numFeatures = data.get(0).length;
        RealMatrix covariance = new Array2DRowRealMatrix(numFeatures, numFeatures);
        for (double[] row : data) {
            for (int i = 0; i < numFeatures; i++) {
                for (int j = 0; j < numFeatures; j++) {
                    covariance.setEntry(i, j, covariance.getEntry(i, j) + (row[i] - mean[i]) * (row[j] - mean[j]));
                }
            }
        }
        for (int i = 0; i < numFeatures; i++) {
            for (int j = 0; j < numFeatures; j++) {
                covariance.setEntry(i, j, covariance.getEntry(i, j) / (data.size() - 1));
            }
        }
        return covariance;
    }

    private static RealMatrix calculatePooledCovariance(RealMatrix cov1, RealMatrix cov2, int n1, int n2) {
        return cov1.scalarMultiply((double) (n1 - 1) / (n1 + n2 - 2)).add(cov2.scalarMultiply((double) (n2 - 1) / (n1 + n2 - 2)));
    }

    private static double[] calculateCoefficients(double[] mean1, double[] mean2, RealMatrix invPooledCov) {
        double[] coefficients = new double[mean1.length];
        double[] difference = new double[mean1.length]; // Array to store element-wise difference
        for (int i = 0; i < mean1.length; i++) {
            difference[i] = mean1[i] - mean2[i]; // Calculate element-wise difference
        }
        ArrayRealVector diffRealVector = new ArrayRealVector(difference);
        for (int i = 0; i < mean1.length; i++) {
            coefficients[i] = invPooledCov.getRowVector(i).dotProduct(diffRealVector);
        }
        return coefficients;
    }

    private static double calculateIntercept(double[] mean1, double[] mean2, double[] coefficients) {
        return -0.5 * (dotProduct(mean1, coefficients) + dotProduct(mean2, coefficients));
    }

    private static double calculateScore(double[] dataPoint, double[] coefficients, double intercept) {
        return dotProduct(dataPoint, coefficients) + intercept;
    }

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

    public static double getCutOff(double fdrCutOff, List<Double> targetScores, List<Double> decoyScores) {
        targetScores.sort(Comparator.reverseOrder());
        List<Integer> decoyIndexes = getDecoyIndexes(targetScores, decoyScores);
        // calculate FDRs for each point
        List<Double> fdrs = new ArrayList<>();
        for (int i = 0; i < decoyIndexes.size(); i++) {
            fdrs.add((double) (i+1)/ (decoyIndexes.get(i)+1));
        }
        int iOfDecoy = getLastValueBelowThreshold(fdrs, fdrCutOff);
        int n = decoyIndexes.get(iOfDecoy);
        if (n-iOfDecoy<=0) {
            return targetScores.get(0);
        }
        return targetScores.get(n-iOfDecoy-1);
    }

    private static List<Integer> getDecoyIndexes(List<Double> targetScores, List<Double> decoyScores) {
        // Combine the two lists into a single list
        List<Double> combinedScores = new ArrayList<>(targetScores);
        combinedScores.addAll(decoyScores);

        // Create a list to store the original indexes of the scores
        List<Integer> originalIndexes = new ArrayList<>();
        for (int i = 0; i < combinedScores.size(); i++) {
            originalIndexes.add(i);
        }

        // Sort the combined scores in descending order, maintaining the original indexes
        originalIndexes.sort(Comparator.comparingDouble(combinedScores::get).reversed());

        // Extract the indexes of the decoy scores in the ranked list
        List<Integer> decoyIndexes = new ArrayList<>();
        for (int i = 0; i < originalIndexes.size(); i++) {
            if (originalIndexes.get(i) >= targetScores.size()) {
                decoyIndexes.add(i);
            }
        }

        return decoyIndexes;
    }

    public static Integer getLastValueBelowThreshold(List<Double> list, double threshold) {
        for (int i = list.size() - 1; i >= 0; i--) {
            if (list.get(i) <= threshold) {
                return i;
            }
        }
        return 0; // No value found
    }

    public static List<Double> calculateScores(List<double[]> data, double[] coefficients, double intercept) {
        // Calculate scores for each data point
        List<Double> scores = new ArrayList<>();
        for (double[] dataPoint : data) {
            double score = calculateScore(dataPoint, coefficients, intercept);
            scores.add(score);
        }
        return scores;
    }

    public static boolean overThreshold(PSM psm, boolean noPASEF) {
        double[] data;
        if (noPASEF) {
            data = new double[]{psm.intensity, psm.avePPM, psm.kl, (double) psm.rtDiff};
        } else {
            data = new double[]{psm.intensity, psm.avePPM, psm.kl, (double) psm.rtDiff, (double) psm.imDiff};
        }
        double c_score = calculateScore(data, coefficients, intercept);
        return c_score > thresholdScore;
    }
}


