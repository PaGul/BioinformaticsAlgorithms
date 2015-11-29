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
public class ViterbiAlgo {

    static double piH = Math.log(0.5) / Math.log(2);
    static double piL = Math.log(0.5) / Math.log(2);

    static class Emmisions {

        double zero;
        double one;

        public Emmisions(double zero, double one) {
            this.zero = zero;
            this.one = one;
        }

        double getValue(char c) {
            switch (c) {
                case '0':
                    return zero;
                case '1':
                    return one;
                default:
                    throw new AssertionError();
            }
        }
    }
    static Emmisions pL;
    static Emmisions pH;
    static double pLL = Math.log(0.5) / Math.log(2);
    static double pHH = Math.log(0.5) / Math.log(2);
    static double pHL = Math.log(0.5) / Math.log(2);
    static double pLH = Math.log(0.5) / Math.log(2);

    static {
        pH = new Emmisions(Math.log(0.5) / Math.log(2), Math.log(0.5) / Math.log(2));
        pL = new Emmisions(Math.log(0.9) / Math.log(2), Math.log(0.1) / Math.log(2));
    }

    static String logHiddenSeq(String seq) {
        double[][] matrix = new double[2][seq.length()];
        boolean[][] whereAreYouFrom = new boolean[2][seq.length()];
        matrix[0][0] = piH + pH.getValue(seq.charAt(0));
        matrix[1][0] = piL + pL.getValue(seq.charAt(0));
        for (int i = 1; i < seq.length(); i++) {
            matrix[0][i] = pH.getValue(seq.charAt(i)) + Math.max(matrix[0][i - 1] + pHH, matrix[1][i - 1] + pLH);
            whereAreYouFrom[0][i] = (matrix[0][i - 1] + pHH > matrix[1][i - 1] + pLH);
            matrix[1][i] = pL.getValue(seq.charAt(i)) + Math.max(matrix[0][i - 1] + pHL, matrix[1][i - 1] + pLL);
            whereAreYouFrom[1][i] = (matrix[0][i - 1] + pHL > matrix[1][i - 1] + pLL);
        }
        StringBuilder res = new StringBuilder();
        int iter = seq.length() - 1;
        int currIndex = (matrix[0][seq.length() - 1] > matrix[1][seq.length() - 1]) ? 0 : 1;
        res.append((currIndex == 0) ? "+" : "-");
        while (iter >= 1) {
            if (whereAreYouFrom[currIndex][iter]) {
                res.append("+");
                currIndex = 0;
            } else {
                res.append("-");
                currIndex = 1;
            }
            iter--;
        }
        return res.reverse().toString();
    }
}
