/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bioinformaticsalgorithms;

/**
 *
 * @author pavelgulaev
 */
public class ForwardBackwardAlgo {

    static double[][] T = new double[][]{{0.8, 0.2}, {0.3, 0.7}};
    static double[][] transT = transposeMatrix(T);
    static double[][] o1 = new double[][]{{0.1, 0}, {0, 0.9}};
    static double[][] o2 = new double[][]{{0.9, 0}, {0, 0.1}};

    static double[][] fStart = new double[][]{{0.5, 0.5}};
    static double[][] bEnd = new double[][]{{1}, {1}};

    static Matrix[] f;
    static Matrix[] b;
    static double[][] res;

    public static void run(String seq) {
        f = new Matrix[seq.length() + 1];
        for (int i = 0; i < f.length; i++) {
            f[i] = new Matrix();
        }
        f[0].matrix = transposeMatrix(fStart);
        for (int i = 1; i <= seq.length(); i++) {
            Matrix Ot = new Matrix();
            if (seq.charAt(i - 1) == '0') {
                Ot.matrix = o1;
            } else {
                Ot.matrix = o2;
            }
            f[i].matrix = multiplyMatrix(multiplyMatrix(Ot.matrix, transT), f[i - 1].matrix);
            /*double sum = f[i].matrix[0][0] + f[i].matrix[1][0];
            f[i].matrix[0][0] /= sum;
            f[i].matrix[1][0] /= sum;*/
        }

        b = new Matrix[seq.length() + 1];
        for (int i = 0; i < b.length; i++) {
            b[i] = new Matrix();
        }
        b[seq.length()].matrix = bEnd;
        for (int i = seq.length() - 1; i >= 0; i--) {
            Matrix Ot = new Matrix();
            if (seq.charAt(i) == '0') {
                Ot.matrix = o1;
            } else {
                Ot.matrix = o2;
            }
            b[i].matrix = multiplyMatrix(multiplyMatrix(T, Ot.matrix), b[i + 1].matrix);
            /*double sum = b[i].matrix[0][0] + b[i].matrix[1][0];
            b[i].matrix[0][0] /= sum;
            b[i].matrix[1][0] /= sum;*/
        }
        res = new double[2][seq.length() + 1];
        double sum = f[f.length-1].matrix[0][0] + f[f.length-1].matrix[1][0];
        for (int i = 1; i <= seq.length(); i++) {
            res[0][i] = f[i].matrix[0][0] * b[i].matrix[0][0];
            res[1][i] = f[i].matrix[1][0] * b[i].matrix[1][0];
            //double sum = res[0][i] + res[1][i];
            res[0][i] /= sum;
            res[1][i] /= sum;
        }
        System.out.print("+ ");
        for (int i = 1; i < res[0].length; i++) {
            System.out.printf("%.2f ", res[0][i]);
        }
        System.out.println("");
        System.out.print("- ");
        for (int i = 1; i < res[0].length; i++) {
            System.out.printf("%.2f ", res[1][i]);
        }
    }

    static class Matrix {

        double[][] matrix;

    }

    public static double[][] transposeMatrix(double[][] m) {
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                temp[j][i] = m[i][j];
            }
        }
        return temp;
    }

    public static double[][] multiplyMatrix(double[][] A, double[][] B) {

        int aRows = A.length;
        int aColumns = A[0].length;
        int bRows = B.length;
        int bColumns = B[0].length;

        if (aColumns != bRows) {
            throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
        }

        double[][] C = new double[aRows][bColumns];

        for (int i = 0; i < aRows; i++) { // aRow
            for (int j = 0; j < bColumns; j++) { // bColumn
                for (int k = 0; k < aColumns; k++) { // aColumn
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }
}
