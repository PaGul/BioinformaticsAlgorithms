/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bioinformaticsalgorithms;

/**
 *
 * @author pavelgulaev Forward-Backward Alignment Algorithm
 */
public class FBAlignment {

    static double delta = 0.1;
    static double tau = 0.1;
    static double epsilon = 0.5;
    static double gamma = 0;
    static double qX = 1;
    static double qY = 1;
    static double pXY = 0.9;

    public static double[][] run(String seq1, String seq2) {
        int n = seq1.length();
        int m = seq2.length();
        double[][] fM = new double[n + 1][m + 1];
        double[][] fX = new double[n + 1][m + 1];
        double[][] fY = new double[n + 1][m + 1];
        double[][] bM = new double[n + 1][m + 1];
        double[][] bX = new double[n + 1][m + 1];
        double[][] bY = new double[n + 1][m + 1];
        double[][] res = new double[n][m];
        double fE;
        fM[0][0] = 1;
        fX[0][0] = 0;
        fY[0][0] = 0;
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                if (i == 0 && j == 0) {
                    continue;
                }
                double p;
                if (i == 0) {
                    fM[i][j] = 0;
                    fX[i][j] = 0;
                    fY[i][j] = qY * (delta * fM[i][j - 1] + epsilon * fY[i][j - 1] + gamma * fX[i][j - 1]);
                } else {
                    if (j == 0) {
                        fM[i][j] = 0;
                        fX[i][j] = qX * (delta * fM[i - 1][j] + epsilon * fX[i - 1][j] + gamma * fY[i - 1][j]);
                        fY[i][j] = 0;
                    } else {
                        if (seq1.charAt(i - 1) != seq2.charAt(j - 1)) {
                            p = 1 - pXY;
                        } else {
                            p = pXY;
                        }
                        fM[i][j] = p * ((1 - 2 * delta - tau) * fM[i - 1][j - 1]
                                + (1 - epsilon - tau - gamma) * (fX[i - 1][j - 1] + fY[i - 1][j - 1]));
                        fX[i][j] = qX * (delta * fM[i - 1][j] + epsilon * fX[i - 1][j] + gamma * fY[i - 1][j]);
                        fY[i][j] = qY * (delta * fM[i][j - 1] + epsilon * fY[i][j - 1] + gamma * fX[i][j - 1]);
                    }
                }
            }
        }
        fE = tau * (fM[n][m] + fX[n][m] + fY[n][m]);
        bM[n][m] = bX[n][m] = bY[n][m] = tau;
        for (int i = n; i >= 1; i--) {
            for (int j = m; j >= 1; j--) {
                if (i == n && j == m) {
                    continue;
                }
                if (i == n) {
                    bM[i][j] = delta * qY * bY[i][j + 1];
                    bX[i][j] = gamma * qY * bY[i][j + 1];
                    bY[i][j] = epsilon * qY * bY[i][j + 1];
                } else {
                    if (j == m) {
                        bM[i][j] = delta * qX * bX[i + 1][j];
                        bX[i][j] = epsilon * bX[i + 1][j];
                        bY[i][j] = gamma * qX * bX[i + 1][j];
                    } else {
                        // i+1 j+1 символы
                        double p = (seq1.charAt(i) == seq2.charAt(j)) ? pXY : (1 - pXY);
                        bM[i][j] = (1 - 2 * delta - tau) * p * bM[i + 1][j + 1]
                                + delta * (qX * bX[i + 1][j] + qY * bY[i][j + 1]);
                        bX[i][j] = (1 - epsilon - tau - gamma) * p * bM[i + 1][j + 1] + epsilon * qX * bX[i + 1][j]
                                + gamma * qY * bY[i][j + 1];
                        bY[i][j] = (1 - epsilon - tau - gamma) * p * bM[i + 1][j + 1] + epsilon * qY * bY[i][j + 1]
                                + gamma * qX * bX[i + 1][j];
                    }
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                res[i][j] = fM[i + 1][j + 1] * bM[i + 1][j + 1] / fE;
            }
        }
        return res;
    }
}
