/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bioinformaticsalgorithms;

import java.util.LinkedList;

/**
 *
 * @author pavelgulaev
 */
public class AffineAlignment {

    private String first;
    private String second;
    private Pair[][] matrix;
    private int[][] ix;
    private int[][] iy;
    private int[][] prevEl;
    public int openGapWeight = -2;
    public int gapWeight = 1;
    public int mismatchWeight = 1;
    private int Weight = 1;
    private int matchWeight = 1;

    public AffineAlignment(String filename) {
        String[] strings = FastaReader.readTwoStrings(filename);
        first = strings[0];
        second = strings[1];
    }

    public AffineAlignment(String first, String second) {
        this.first = first;
        this.second = second;
    }

        static class Pair {
            int value;
            boolean opened = false;
        }
    
    public Pair[][] createAlignmentMatrix() {
        matrix = new Pair[first.length() + 1][second.length() + 1];
        ix = new int[first.length() + 1][second.length() + 1];
        iy = new int[first.length() + 1][second.length() + 1];
        int[][] backm = new int[first.length() + 1][second.length() + 1];
        int[][] backix = new int[first.length() + 1][second.length() + 1];
        int[][] backiy = new int[first.length() + 1][second.length() + 1];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] = new Pair();
            }
        }
        for (int i = 1; i < matrix.length; i++) {
            ix[i][0] = -openGapWeight - gapWeight * (i - 1);
            matrix[i][0].value = -openGapWeight - gapWeight * (i - 1);
            iy[i][0] = -10000000;
        }
        for (int i = 1; i < matrix[0].length; i++) {
            ix[0][i] = -10000000;
            matrix[0][i].value = -openGapWeight - gapWeight * (i - 1);
            iy[0][i] = -openGapWeight - gapWeight * (i - 1);
        }
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
                if (first.charAt(i - 1) != second.charAt(j - 1)) {
                    Weight = -mismatchWeight;
                } else {
                    Weight = matchWeight;
                }
//                Weight = Blosum.getDistance(first.charAt(i-1), second.charAt(j-1));
                int tempix = (matrix[i - 1][j].opened)?-100000000:matrix[i - 1][j].value;
                Integer[] tempx = {ix[i - 1][j] - gapWeight,
                    tempix - (openGapWeight)};
                ix[i][j] = max(tempx);
                backix[i][j] = maxIndex(tempx);
                int tempiy = (matrix[i][j - 1].opened)?-100000000:matrix[i][j - 1].value;
                Integer[] tempy = {iy[i][j - 1] - gapWeight,
                    tempiy - (openGapWeight)};
                iy[i][j] = max(tempy);
                backiy[i][j] = maxIndex(tempy);
                Integer[] tempcentr = {ix[i][j],
                    matrix[i - 1][j - 1].value + Weight,
                    iy[i][j]};
                matrix[i][j].value = max(tempcentr);
                backm[i][j] = maxIndex(tempcentr);
                if (backm[i][j] != 1) {
                    matrix[i][j].opened = true;
                }
            }
        }
        int i = matrix.length - 1;
        int j = matrix[0].length - 1;
        StringBuilder firstAl = new StringBuilder(first);
        StringBuilder secAl = new StringBuilder(second);
        System.out.println(Math.max(matrix[i][j].value, Math.max(ix[i][j], iy[i][j])));
        int backIndex = maxIndex(ix[i][j], matrix[i][j].value, iy[i][j]);
        while (i != 0 && j != 0) {
            switch (backIndex) {
                case 0:
                    if (backix[i][j] == 1) {
                        backIndex = 1;
                    }
                    i--;
                    secAl.replace(j, j, "-");
                    break;
                case 1:
                    if (backIndex == 1) {
                        if (backm[i][j] == 0) {
                            backIndex = 0;
                        } else {
                            if (backm[i][j] == 2) {
                                backIndex = 2;
                            } else {
                                i--;
                                j--;
                            }
                        }
                    }
                    break;
                case 2:
                    if (backiy[i][j] == 1) {
                        backIndex = 1;
                    }
                    j--;
                    firstAl.replace(i, i, "-");
                    break;
            }

        }
        for (int k = 0; k < i; k++) {
            secAl.replace(0, 0, "-");
        }
        for (int k = 0; k < j; k++) {
            firstAl.replace(0, 0, "-");
        }
        System.out.println(firstAl);
        System.out.println(secAl);
        return matrix;
    }

    // номер в переданном массиве
    private int maxIndex(Integer... a) {
        int max = Integer.MIN_VALUE;
        int maxIndex = -1;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    private int max(Integer... a) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }

}
