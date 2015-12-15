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
public class RNAFolding {

    private static String code;
    static int[][] matrix;
    public static void run(String f) {
        code = f;
        matrix = new int[f.length()][f.length()];
        
        for (int i = 1, j = 0; i < f.length() && j < f.length(); i++, j++) {
            matrix[i][j] = 0;
        }
        for (int i = 0, j = 0; i < f.length() && j < f.length(); i++, j++) {
            matrix[i][j] = 0;
        }
        for (int iter = 3; iter < f.length(); iter++) {
            for (int i = 0, j = iter; i < f.length() && j < f.length(); i++, j++) {
                int firstCand = (basePair(i, j)) ? (matrix[i + 1][j - 1] + 1) : (matrix[i + 1][j - 1]);
                int max = Integer.MIN_VALUE;
                for (int k = i + 1; k < j; k++) {
                    if (matrix[i][k] + matrix[k + 1][j] > max) {
                        max = matrix[i][k] + matrix[k + 1][j];
                    }
                }
                matrix[i][j] = Math.max(Math.max(firstCand, max), Math.max(matrix[i + 1][j], matrix[i][j - 1]));
            }
        }
//        for (int i = 0; i < f.length(); i++) {
//            for (int j = 0; j < f.length(); j++) {
//                System.out.print(matrix[i][j] + " ");
//            }
//            System.out.println("");
//        }

        class Coord {

            int i;
            int j;
            int value;

            public Coord(int i, int j) {
                this.i = i;
                this.j = j;
            }

            public Coord(int i, int j, int value) {
                this.i = i;
                this.j = j;
                this.value = value;
            }

            @Override
            public String toString() {
                return ""+i+" "+j+": "+ f.charAt(i)+"=="+f.charAt(j);
            }
            
        }

        LinkedList<Coord> stack = new LinkedList<>();
        LinkedList<Coord> res = new LinkedList<>();
        stack.push(new Coord(0, f.length() - 1));
        while (!stack.isEmpty()) {
            Coord temp = stack.pop();
            int i = temp.i;
            int j = temp.j;
            if (i >= j) {
                continue;
            } else {
                if (matrix[i + 1][j] == matrix[i][j]) {
                    stack.push(new Coord(i + 1, j));
                } else {
                    if (matrix[i][j - 1] == matrix[i][j]) {
                        stack.push(new Coord(i, j - 1));
                    } else {
                        int basePair = (basePair(i, j)) ? 1 : 0;
                        if (matrix[i + 1][j - 1] + basePair == matrix[i][j]) {
                            res.add(new Coord(i, j, matrix[i][j]));
                            stack.push(new Coord(i + 1, j - 1));
                        } else {
                            for (int k = i + 1; k < j; k++) {
                                if (matrix[i][k] + matrix[k + 1][j] == matrix[i][j]) {
                                    stack.push(new Coord(k+1,j));
                                    stack.push(new Coord(i,k));
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (Coord coord : res) {
            System.out.println(coord);
        }
    }
    
    private void traceback(int i, int j) {
        if (j<=i) {
            return;
        }
        else {
            if (matrix[i][j]==matrix[i][j-1]) {
                traceback(i, j-1);
                return;
            } else {
                
            }
        }
    }

    private static boolean basePair(int i, int j) {
        if ((code.charAt(i) == 'U' && code.charAt(j) == 'A') || (code.charAt(i) == 'A' && code.charAt(j) == 'U')) {
            return true;
        }
        if ((code.charAt(i) == 'G' && code.charAt(j) == 'C') || (code.charAt(i) == 'C' && code.charAt(j) == 'G')) {
            return true;
        }
        return false;
    }

}
