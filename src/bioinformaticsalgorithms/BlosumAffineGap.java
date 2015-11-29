/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package bioinformaticsalgorithms;

import java.util.Scanner;

/**
 *
 * @author pavelgulaev
 */
public class BlosumAffineGap {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        sc.next();
        String first = sc.next();
        sc.next();
        String second = sc.next();
        AffineAlignmentB aa = new AffineAlignmentB(first, second);
        aa.printAlignment();
    }
    
}

class AffineAlignmentB {

    private String first;
    private String second;
    private int[][] matrix;
    private int[][] ix;
    private int[][] iy;
    private int[][] prevEl;
    public int openGapWeight = 11;
    public int gapWeight = 1;
//    public int mismatchWeight = 1;
    private int Weight;
//    private int matchWeight = 1;

//    public AffineAlignment(String filename) {
//        String[] strings = FastaReader.readTwoStrings(filename);
//        first = strings[0];
//        second = strings[1];
//    }

    public AffineAlignmentB(String first, String second) {
        this.first = first;
        this.second = second;
    }
    
    public AffineAlignmentB(String filename) {
        String[] strings = FastaReader.readTwoStrings(filename);
        first = strings[0];
        second = strings[1];
    }
    int max = Integer.MIN_VALUE;
    public int[][] createAlignmentMatrix() {
        matrix = new int[first.length() + 1][second.length() + 1];
        ix = new int[first.length() + 1][second.length() + 1];
        iy = new int[first.length() + 1][second.length() + 1];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
//                matrix[i][j] = -100000;
//                ix[i][j] = -100000;
//                iy[i][j] = -100000;
            }
        }
        for (int i = 1; i < matrix.length; i++) {
            matrix[i][0] = -openGapWeight-i*gapWeight;
            ix[i][0] = -openGapWeight - gapWeight * (i);
        }
        for (int i = 1; i < matrix[0].length; i++) {
            matrix[0][i] = -openGapWeight-i*gapWeight;
            iy[0][i] = -openGapWeight - gapWeight * (i);
        }
    
        matrix[0][0] = 0;
//        ix[0][0] = 0;//-openGapWeight;
//        iy[0][0] = 0;//-openGapWeight;
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
//                if (first.charAt(i - 1) != second.charAt(j - 1)) {
//                    Weight = -mismatchWeight;
//                } else {
//                    Weight = matchWeight;
//                }
                Weight = Blosum.getDistance(first.charAt(i - 1), second.charAt(j - 1));
                ix[i][j] = Math.max(
                        matrix[i - 1][j] - (openGapWeight + gapWeight),
//                        Math.max(iy[i - 1][j] - (openGapWeight + gapWeight), 
                        ix[i - 1][j] - gapWeight);
                iy[i][j] = Math.max(
                        matrix[i][j - 1] - (openGapWeight + gapWeight), 
//                        Math.max(ix[i][j - 1] - (openGapWeight + gapWeight), 
                        iy[i][j - 1] - gapWeight);
                matrix[i][j] = Math.max(Math.max(
                        matrix[i - 1][j - 1] + Weight,
                        ix[i][j]),
                        iy[i][j]);
            }
        }
        System.out.println(matrix[matrix.length-1][matrix[0].length-1]);
        return matrix;
    }

    public void printAlignment() {
        matrix = createAlignmentMatrix();
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println("");
        }
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(ix[i][j] + " ");
            }
            System.out.println("");
        }
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(iy[i][j] + " ");
            }
            System.out.println("");
        }
        String firstAl = "";
        String secondAl = "";
        int i = matrix.length - 1;
        int j = matrix[0].length - 1;
        while ((i != 0) && (j != 0)) {
            int max = Math.max(matrix[i][j], Math.max(ix[i][j], iy[i][j]));
            if (max == matrix[i][j]) {
                firstAl = first.charAt(i - 1) + firstAl;
                secondAl = second.charAt(j - 1) + secondAl;
                i--;
                j--;
            }
            else {
                if (max == ix[i][j]) {
                    firstAl = first.charAt(i - 1) + firstAl;
                    secondAl = "-" + secondAl;
                    i--;
                    
                } else {
                    firstAl = "-" + firstAl;
                    secondAl = second.charAt(j - 1) + secondAl;
                    j--;
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
//        System.out.println(max);
        System.out.println(firstAl);
        System.out.println(secondAl);
    }
}


class Blosum{

    /*
     * Array representation of Blosum-62 matrix 
     * Refer to above matrix for corrseponding amino acids
     * i.e. score(A, R) corresponds to  matrix[0][1]=matrix[1][0]=-1
    */  
    private static final int[][] matrix = {
	{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
	{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
	{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
	{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
	{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
	{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
	{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
	{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
	{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
	{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
	{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
	{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
	{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
	{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
	{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
	{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
	{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
	{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
	{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
	{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};

	// quick and dirty equivalent of typesafe enum pattern, can also use HashMap
    // or even better, EnumMap in Java 5. 
    // This code is for Java 1.4.2, so we will stick to the simple implementation
    private static int getIndex(char a) {
    	// check for upper and lowercase characters
    	switch ((String.valueOf(a)).toUpperCase().charAt(0)) {
	    	case 'A': return 0;
	    	case 'R': return 1;
	    	case 'N': return 2;
	    	case 'D': return 3;
	    	case 'C': return 4;
	    	case 'Q': return 5;
	    	case 'E': return 6;
	    	case 'G': return 7;
	    	case 'H': return 8;
	    	case 'I': return 9;
	    	case 'L': return 10;
	    	case 'K': return 11;
	    	case 'M': return 12;
	    	case 'F': return 13;
	    	case 'P': return 14;
	    	case 'S': return 15;
	    	case 'T': return 16;
	    	case 'W': return 17;
	    	case 'Y': return 18;
	    	case 'V': return 19;
	    	default: return -1;
    	}
    }
    
    private static int getDistance(int i, int j) {
    	return matrix[i][j];
    }

    public static int getDistance(char a1, char a2) {
    	// toUpper
    	return getDistance(getIndex(a1), getIndex(a2));  	
    }
}

