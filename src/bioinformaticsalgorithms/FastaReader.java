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
public class FastaReader {
    public static String[] readTwoStrings(String filename) {
        String first = "";
        String second = "";
        try {
            BufferedReader in = new BufferedReader(
                new FileReader(filename));
        int temp;
        String code = "";
        in.readLine();
        List<String> input = new LinkedList<>();
        while ((temp = in.read()) != -1) {
            char letter = (char) temp;
            if (letter == '>') {
                if (!code.equals("")) {
                    input.add(code);
                }
                code = "";
                in.readLine();
            } else {
                code += letter;
                code += in.readLine();
            }
        }
        input.add(code);
        first = input.get(0);
        second = input.get(1);
        } catch (IOException e) {
            System.out.println("Problems with fasta file");
        }
        return new String[]{first, second};
    }
}
