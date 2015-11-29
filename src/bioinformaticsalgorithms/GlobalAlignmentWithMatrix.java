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
public class GlobalAlignmentWithMatrix {
    private String first;
    private String second;
    private int[][] matrix;
    private int[][] wMatrix = {{100,-100,-100,-100,-49},
                                {0,100,-100,-100,-49},
                                {0,0,100,-97,-49},
                                {0,0,0,100,-49}};
    public GlobalAlignmentWithMatrix(String filename) {
        String[] strings = FastaReader.readTwoStrings(filename);
        first = strings[0];
        second = strings[1];
    }

    
    public GlobalAlignmentWithMatrix(String first, String second) {
        this.first = first;
        this.second = second;
    }

    public void setMatrix(int[][] matrix) {
        this.wMatrix = matrix;
    }
    
    private int aminoToNumber(char amino) {
        switch(amino) {
            case 'A':
                return 0;
            case 'T':
                return 1;
            case 'G':
                return 2;
            case 'C':
                return 3;
            case '-':
                return 4;
        }
        return -1;
    }
    
    public int[][] createAlignmentMatrix() {
        matrix = new int[first.length()+1][second.length()+1];
        matrix[0][0] = 0;
        for (int i = 1; i < matrix.length; i++) {
            matrix[i][0] = wMatrix[aminoToNumber(first.charAt(i-1))][aminoToNumber('-')];
        }
        for (int i = 1; i < matrix[0].length; i++) {
            matrix[0][i] = wMatrix[aminoToNumber(second.charAt(i-1))][aminoToNumber('-')];
        }
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
                int fMatrCoord = aminoToNumber(first.charAt(i-1));
                int sMatrCoord = aminoToNumber(second.charAt(j-1));
                int matchWeight = (fMatrCoord<sMatrCoord)?wMatrix[fMatrCoord][sMatrCoord]:
                        wMatrix[sMatrCoord][fMatrCoord];
                int gapCoord = aminoToNumber('-');
                int insWeight= (fMatrCoord<gapCoord)?wMatrix[fMatrCoord][gapCoord]:
                        wMatrix[gapCoord][fMatrCoord];
                int delweight = (sMatrCoord<gapCoord)?wMatrix[sMatrCoord][gapCoord]:
                        wMatrix[gapCoord][sMatrCoord];
                matrix[i][j] = Math.max(Math.max(
                        matrix[i-1][j]+insWeight, 
                        matrix[i][j-1]+delweight),
                        matrix[i-1][j-1]+matchWeight);
            }
        }
        return matrix;
    }
    
     public void printAlignment() {
        matrix = createAlignmentMatrix();
        String firstAl = "";
        String secondAl = "";
        int i = matrix.length-1;
        int j = matrix[0].length-1;
        while ((i!=0) && (j!=0)) {
            int fMatrCoord = aminoToNumber(first.charAt(i-1));
                int sMatrCoord = aminoToNumber(second.charAt(j-1));
                int matchWeight = (fMatrCoord<sMatrCoord)?wMatrix[fMatrCoord][sMatrCoord]:
                        wMatrix[sMatrCoord][fMatrCoord];
                int gapCoord = aminoToNumber('-');
                int insWeight= (fMatrCoord<gapCoord)?wMatrix[fMatrCoord][gapCoord]:
                        wMatrix[gapCoord][fMatrCoord];
                int delweight = (sMatrCoord<gapCoord)?wMatrix[sMatrCoord][gapCoord]:
                        wMatrix[gapCoord][sMatrCoord];
            if (matrix[i][j]-matchWeight==matrix[i-1][j-1]) {
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
                if (matrix[i][j]-delweight==matrix[i][j-1]) {
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
