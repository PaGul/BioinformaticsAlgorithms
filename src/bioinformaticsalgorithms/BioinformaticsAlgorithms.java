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
public class BioinformaticsAlgorithms {

    public static void main(String[] args) {
        WUPGMA w = new WUPGMA();
        System.out.println("WPGMA");
        w.runWPGMA();
        System.out.println("UPGMA");
        w.runUPGMA();
    }
    
}
