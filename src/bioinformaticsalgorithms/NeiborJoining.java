/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bioinformaticsalgorithms;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Set;

/**
 *
 * @author pavelgulaev
 */
public class NeiborJoining {

    public void run() {
//        LinkedList<Val> matrix = new LinkedList<>();

        Set<Node> letters = new HashSet<Node>(Arrays.asList(new Node[]{new Node("A"), new Node("B"), new Node("C"), new Node("D"), new Node("E"), new Node("F")}));
        HashMap<Node, Double> heightsOfNodes = new HashMap<>();
        LinkedList<Node> tree = new LinkedList<>();
        LinkedList<Node> currNodes = new LinkedList<>();
        for (Node letter : letters) {
            heightsOfNodes.put(letter, 0.0);
            currNodes.add(letter);
        }

        HashMap<Node, HashMap<Node, Double>> matrix = new HashMap<>();
//        HashMap<Node, Double> temp = new HashMap<>();
//        temp.put(new Node("B"), 5.0);
//        temp.put(new Node("C"), 5.0);
//        temp.put(new Node("D"), 5.0);
//        matrix.put(new Node("A"), temp);
//
//        temp = new HashMap<>();
//        temp.put(new Node("A"), 5.0);
//        temp.put(new Node("C"), 2.0);
//        temp.put(new Node("D"), 2.0);
//        matrix.put(new Node("B"), temp);
//
//        temp = new HashMap<>();
//        temp.put(new Node("A"), 5.0);
//        temp.put(new Node("B"), 2.0);
//        temp.put(new Node("D"), 1.0);
//        matrix.put(new Node("C"), temp);
//
//        temp = new HashMap<>();
//        temp.put(new Node("A"), 5.0);
//        temp.put(new Node("B"), 2.0);
//        temp.put(new Node("C"), 1.0);
//        matrix.put(new Node("D"), temp);
//
////        temp = new HashMap<>();
////        temp.put(new Node("A"), 8.0);
////        temp.put(new Node("B"), 9.0);
////        temp.put(new Node("C"), 7.0);
////        temp.put(new Node("D"), 3.0);
////        matrix.put(new Node("E"), temp);

        HashMap<Node, Double> temp = new HashMap<>();
        temp.put(new Node("B"), 5.0);
        temp.put(new Node("C"), 4.0);
        temp.put(new Node("D"), 7.0);
        temp.put(new Node("E"), 6.0);
        temp.put(new Node("F"), 8.0);
        matrix.put(new Node("A"), temp);

        temp = new HashMap<>();
        temp.put(new Node("A"), 5.0);
        temp.put(new Node("C"), 7.0);
        temp.put(new Node("D"), 10.0);
        temp.put(new Node("E"), 9.0);
        temp.put(new Node("F"), 11.0);
        matrix.put(new Node("B"), temp);

        temp = new HashMap<>();
        temp.put(new Node("A"), 4.0);
        temp.put(new Node("B"), 7.0);
        temp.put(new Node("D"), 7.0);
        temp.put(new Node("E"), 6.0);
        temp.put(new Node("F"), 8.0);
        matrix.put(new Node("C"), temp);

        temp = new HashMap<>();
        temp.put(new Node("A"), 7.0);
        temp.put(new Node("B"), 10.0);
        temp.put(new Node("C"), 7.0);
        temp.put(new Node("E"), 5.0);
        temp.put(new Node("F"), 9.0);
        matrix.put(new Node("D"), temp);

        temp = new HashMap<>();
        temp.put(new Node("A"), 6.0);
        temp.put(new Node("B"), 9.0);
        temp.put(new Node("C"), 6.0);
        temp.put(new Node("D"), 5.0);
        temp.put(new Node("F"), 8.0);
        matrix.put(new Node("E"), temp);

        temp = new HashMap<>();
        temp.put(new Node("A"), 8.0);
        temp.put(new Node("B"), 11.0);
        temp.put(new Node("C"), 8.0);
        temp.put(new Node("D"), 9.0);
        temp.put(new Node("E"), 8.0);
        matrix.put(new Node("F"), temp);
        Node newLetter = null;
        while (letters.size() > 2) {
            double min = Integer.MAX_VALUE;
            Value minVal = null;
            for (Node fLetter : letters) {
                for (Node sLetter : letters) {
                    if (fLetter.equals(sLetter)) {
                        continue;
                    }
                    double Q = (letters.size() - 2) * matrix.get(fLetter).get(sLetter);
                    for (Double distance : matrix.get(fLetter).values()) {
                        Q -= distance;
                    }
                    for (Double distance : matrix.get(sLetter).values()) {
                        Q -= distance;
                    }
                    if (Q < min) {
                        min = Q;
                        minVal = new Value(fLetter, sLetter, min);
                    }
                }
            }

            double fuDistance = 0.5 * matrix.get(minVal.first).get(minVal.second);
            double secondCooef = 0;
            for (Double distance : matrix.get(minVal.first).values()) {
                secondCooef += distance;
            }
            for (Double distance : matrix.get(minVal.second).values()) {
                secondCooef -= distance;
            }
            fuDistance -= (secondCooef * 0.5 / (letters.size() - 2));
            double guDistance = matrix.get(minVal.first).get(minVal.second) - fuDistance;
            minVal.first.height = guDistance;
            minVal.second.height = fuDistance;
            tree.add(minVal.first);
            tree.add(minVal.second);
            temp = new HashMap<>();
            newLetter = new Node(minVal.first, minVal.second);
            matrix.put(new Node(minVal.first, minVal.second), temp);

            for (Node node : letters) {
                if (node.equals(minVal.first) || node.equals(minVal.second)) {
                    continue;
                }
                double dist = 0.5 * (matrix.get(minVal.first).get(node) + matrix.get(minVal.second).get(node)
                        - matrix.get(minVal.first).get(minVal.second));
                temp = matrix.get(node);
                temp.put(newLetter, dist);
                matrix.put(node, temp);
                temp = matrix.get(newLetter);
                temp.put(node, dist);
                matrix.put(newLetter, temp);
            }
            matrix.remove(minVal.first);
            matrix.remove(minVal.second);
            letters.remove(minVal.first);
            letters.remove(minVal.second);
            letters.add(newLetter);
            for (HashMap<Node, Double> hashMap : matrix.values()) {
                hashMap.remove(minVal.first);
                hashMap.remove(minVal.second);
            }
        }
        Entry<Node, HashMap<Node, Double>> lastVal = matrix.entrySet().iterator().next();
        Entry<Node, Double> lastVal2 = lastVal.getValue().entrySet().iterator().next();
        lastVal2.getKey().height = lastVal2.getValue() / 2;
        lastVal.getKey().height = lastVal2.getValue() / 2;
        Node root = new Node(lastVal.getKey(), lastVal2.getKey(), 0);
        tree.add(root);
        System.out.println(inOrderNewick(root));
    }

    public String inOrderNewick(Node root) {
        if (root.hasChild) {
            String output = "";
            output += "(";
            output += inOrderNewick(root.child1);
            output += ",";
            output += inOrderNewick(root.child2);
            output += ")";
            if (root.height != 0) {
                output += ":" + root.height;
            }
            return output;
        } else {
            return root.getSeq();
        }
    }
}

class Val {

    String first;
    String second;
    double dist;

    public Val(String first, String second, double dist) {
        this.first = first;
        this.second = second;
        this.dist = dist;
    }

}
