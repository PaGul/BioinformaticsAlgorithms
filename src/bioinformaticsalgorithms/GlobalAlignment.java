/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package bioinformaticsalgorithms;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author pavelgulaev
 */
public class GlobalAlignment {
    private String first;
    private String second;
    private int[][] matrix;
    private int[][] prevEl;
    public int gapWeight =51;
    public int mismatchWeight = 100;
    private int Weight = 100;
    private int matchWeight = 100;
    public GlobalAlignment(String filename) {
        String[] strings = FastaReader.readTwoStrings(filename);
        first = strings[0];
        second = strings[1];
    }

    public GlobalAlignment(String first, String second) {
        this.first = first;
        this.second = second;
    }
    
    public int[][] createAlignmentMatrix() {
        matrix = new int[first.length()+1][second.length()+1];
        prevEl = new int[first.length()+1][second.length()+1];
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][0] = -i;
        }
        for (int i = 0; i < matrix[0].length; i++) {
            matrix[0][i] = -i;
        }
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
                if (first.charAt(i-1)!=second.charAt(j-1)) {
                    Weight = -mismatchWeight;
                } else {
                    Weight = matchWeight;
                }
                matrix[i][j] = Math.max(Math.max(
                        matrix[i-1][j]-gapWeight, 
                        matrix[i][j-1]-gapWeight),
                        matrix[i-1][j-1]+Weight);
            }
        }
        return matrix;
    }
    
    public void printAlignment() {
        matrix = createAlignmentMatrix();
//        for (int i = 0; i < matrix.length; i++) {
//            for (int j = 0; j < matrix[0].length; j++) {
//                System.out.print(matrix[i][j]+" ");
//            }
//            System.out.println("");
//        }
        String firstAl = "";
        String secondAl = "";
        int i = matrix.length-1;
        int j = matrix[0].length-1;
        while ((i!=0) && (j!=0)) {
            if (first.charAt(i-1)!=second.charAt(j-1)) {
                    Weight = -mismatchWeight;
                } else {
                    Weight = matchWeight;
                }
            
            if (matrix[i][j]==matrix[i-1][j-1]+Weight) {
                if (first.charAt(i-1)==second.charAt(j-1)) {
                    firstAl=first.toUpperCase().charAt(i-1)+firstAl;
                    secondAl=second.toUpperCase().charAt(j-1)+secondAl;
                } else {
                    firstAl=first.toLowerCase().charAt(i-1)+firstAl;
                    secondAl=second.toLowerCase().charAt(j-1)+secondAl;
                }
                i--;
                j--;
            } else {
                if (matrix[i][j]+gapWeight==matrix[i][j-1]) {
                    firstAl="-"+firstAl;
                    secondAl=second.charAt(j-1)+secondAl;
                    j--;
                } else {
                    firstAl=first.charAt(i-1)+firstAl;
                    secondAl="-"+secondAl;
                    i--;
                }
            }
        }
        if (i>0) {
            for (int k = i; k > 0; k--) {
                firstAl=first.charAt(k-1)+firstAl;
                secondAl="-"+secondAl;
            }
        }
        if (j>0) {
            for (int k = j; k > 0; k--) {
                firstAl="-"+firstAl;
                secondAl=second.charAt(k-1)+secondAl;
            }
        }
        System.out.println(firstAl);
        System.out.println(secondAl);
    }
}