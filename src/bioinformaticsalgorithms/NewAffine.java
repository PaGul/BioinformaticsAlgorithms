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
public class NewAffine {

    private String f;
    private String s;

    public NewAffine(String filename) {
        String[] strings = FastaReader.readTwoStrings(filename);
        f = strings[0];
        s = strings[1];
    }

    public NewAffine(String first, String second) {
        this.f = first;
        this.s = second;
    }
    int sigma = 11;
    int epsilon = 1;

    public void affine() {
        int[][][] S = new int[3][f.length() + 1][s.length() + 1];
        int[][][] backtrack = new int[3][f.length() + 1][s.length() + 1];
        for (int i = 1; i < f.length() + 1; i++) {
            S[0][i][0] = -sigma - (i - 1) * epsilon;
            S[1][i][0] = -sigma - (i - 1) * epsilon;
            S[2][i][0] = -10 * sigma; //-inf
        }
        for (int j = 1; j < s.length() + 1; j++) {
            S[2][0][j] = -sigma - (j - 1) * epsilon;
            S[1][0][j] = -sigma - (j - 1) * epsilon;
            S[0][0][j] = -10 * sigma;
        }
        for (int i = 0; i < f.length() + 1; i++) {
            for (int j = 0; j < s.length() + 1; j++) {
                S[0][i][j] = Math.max(S[0][i - 1][j] - epsilon, S[1][i - 1][j] - sigma);
                backtrack[0][i][j] = (S[0][i][j] == S[0][i - 1][j] - epsilon) ? 0 : 1;

                S[2][i][j] = Math.max(S[2][i][j - 1] - epsilon, S[1][i][j - 1] - sigma);
                backtrack[2][i][j] = (S[2][i][j] == S[2][i][j - 1] - epsilon) ? 0 : 1;
                int Weight = Blosum.getDistance(f.charAt(i - 1), s.charAt(j - 1));
                S[1][i][j] = Math.max(S[0][i][j], Math.max(S[1][i - 1][j - 1] + Weight, S[2][i][j]));
                if (S[1][i][j] == S[0][i][j]) {
                    backtrack[1][i][j] = 0;
                } else {
                    if (S[1][i][j] == S[1][i - 1][j - 1] + Weight) {
                        backtrack[1][i][j] = 1;
                    } else {
                        backtrack[1][i][j] = 2;
                    }
                }
            }
        }
        int i =  f.length();
        int j = s.length();
        String v_aligned = f;
        String w_aligned = s;
        int max_score = Math.max(S[0][i][j], Math.max(S[1][i][j], S[2][i][j]));
        int backtrack_matrix = 0;
        if (max_score == S[0][i][j]) {
            backtrack_matrix = 0;
        } else {
            if (max_score == S[1][i][j]) {
                backtrack_matrix = 1;
            } else {
                backtrack_matrix = 2;
            }
        }


    }
}
