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
public class LocalAlignment {

    private String first;
    private String second;
    private int[][] matrix;
    private int[][] prevEl;
    public int gapWeight = 1;//51;
    public final int mismatchWeight = 2;//100;
    private int Weight = 1;//100;
    private final int matchWeight = 1;//100;

    public LocalAlignment(String filename) {
        String[] strings = FastaReader.readTwoStrings(filename);
        first = strings[0];
        second = strings[1];
    }

    public LocalAlignment(String first, String second) {
        this.first = first;
        this.second = second;
    }

    int max = 0;
    int maxI = 0;
    int maxJ = 0;

    public int[][] createAlignmentMatrix() {
        matrix = new int[first.length() + 1][second.length() + 1];
        prevEl = new int[first.length() + 1][second.length() + 1];
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][0] = 0;
        }
        for (int i = 0; i < matrix[0].length; i++) {
            matrix[0][i] = 0;
        }
        max = 0;
        maxI = 0;
        maxJ = 0;
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
                if (first.charAt(i - 1) != second.charAt(j - 1)) {
                    Weight = -mismatchWeight;
                } else {
                    Weight = matchWeight;
                }
                matrix[i][j] = Math.max(Math.max(
                        matrix[i - 1][j] - gapWeight,
                        matrix[i][j - 1] - gapWeight),
                        Math.max(matrix[i - 1][j - 1] + Weight, 0));
                if (matrix[i][j] >= max) {
                    max = matrix[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }
        return matrix;
    }

    public void printAlignment() {
        matrix = createAlignmentMatrix();
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(matrix[i][j]+" ");
            }
            System.out.println("");
        }
        String firstAl = "";
        String secondAl = "";
        int i = maxI;
        int j = maxJ;
        while (matrix[i][j] != 0) {
            if (first.charAt(i - 1) != second.charAt(j - 1)) {
                Weight = mismatchWeight;
            } else {
                Weight = -matchWeight;
            }

            if (matrix[i][j] + Weight == matrix[i - 1][j - 1]) {
                if (first.charAt(i - 1) == second.charAt(j - 1)) {
                    firstAl = first.toUpperCase().charAt(i - 1) + firstAl;
                    secondAl = second.toUpperCase().charAt(j - 1) + secondAl;
                } else {
                    firstAl = first.charAt(i - 1) + firstAl;
                    secondAl = second.charAt(j - 1) + secondAl;
                }
                i--;
                j--;
            } else {
                if (matrix[i][j] + gapWeight == matrix[i][j - 1]) {
                    firstAl = "-" + firstAl;
                    secondAl = second.charAt(j - 1) + secondAl;
                    j--;
                } else {
                    firstAl = first.charAt(i - 1) + firstAl;
                    secondAl = "-" + secondAl;
                    i--;
                }
            }
        }
        String firstStart = first.substring(0,i);
        String secondStart = second.substring(0,j);
        for (int k = 0; k < j - i; k++) {
            firstStart=" "+firstStart;
        }
        for (int k = 0; k < i - j; k++) {
            secondStart=" "+secondStart;
        }
        System.out.println(firstStart.toLowerCase()+firstAl+first.substring(maxI).toLowerCase());
        System.out.println(secondStart.toLowerCase()+secondAl+second.substring(maxJ).toLowerCase());
    }
}
