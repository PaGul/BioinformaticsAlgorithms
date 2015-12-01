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
import java.util.PriorityQueue;
import java.util.Set;

/**
 *
 * @author pavelgulaev
 */
public class WUPGMA {

    public void runUPGMA() {
        PriorityQueue<Value> matrix = new PriorityQueue<Value>(new Comparator<Value>() {

            @Override
            public int compare(Value o1, Value o2) {
                if (o1.dist < o2.dist) {
                    return -1;
                } else {
                    if (o1.dist == o2.dist) {
                        return 0;
                    } else {
                        return 1;
                    }
                }
            }

        });
        Set<String> letters = new HashSet<String>(Arrays.asList(new String[]{"A", "B", "C", "D", "E", "F"}));
//        Set<String> letters = new HashSet<String>(Arrays.asList(new String[]{"K", "L", "M", "N"}));
        HashMap<String, Double> heightsOfNodes = new HashMap<>();
        LinkedList<Node> currNodes = new LinkedList<>();
        for (String letter : letters) {
            heightsOfNodes.put(letter, 0.0);
            currNodes.add(new Node(letter));
        }
//        matrix.add(new Value("K", "L", 16));
//        matrix.add(new Value("K", "M", 16));
//        matrix.add(new Value("K", "N", 10));
//        matrix.add(new Value("L", "M", 8));
//        matrix.add(new Value("L", "N", 8));
//        matrix.add(new Value("M", "N", 4));
        matrix.add(new Value("B", "A", 4));
        matrix.add(new Value("C", "A", 4));
        matrix.add(new Value("D", "A", 7));
        matrix.add(new Value("E", "A", 6));
        matrix.add(new Value("F", "A", 8));
        matrix.add(new Value("C", "B", 7));
        matrix.add(new Value("D", "B", 10));
        matrix.add(new Value("E", "B", 9));
        matrix.add(new Value("F", "B", 11));
        matrix.add(new Value("D", "C", 7));
        matrix.add(new Value("E", "C", 6));
        matrix.add(new Value("F", "C", 8));
        matrix.add(new Value("E", "D", 5));
        matrix.add(new Value("F", "D", 9));
        matrix.add(new Value("F", "E", 8));

        Node newNode = null;
        LinkedList<Node> tree = new LinkedList<>();
        while (!matrix.isEmpty()) {
            Value min = matrix.poll();
            newNode = new Node(min.first, min.second);

            String nameOfNode = min.first.letter + min.second.letter;
            heightsOfNodes.put(nameOfNode, min.dist / 2);
            min.first.height = min.dist / 2 - heightsOfNodes.get(min.first.letter);
            tree.add(min.first);
            min.second.height = min.dist / 2 - heightsOfNodes.get(min.second.letter);
            tree.add(min.second);
            for (Node c : currNodes) {
                if (c.letter.equals(min.first.letter) || c.letter.equals(min.second.letter)) {
                    continue;
                }
                double cAndNewNode = countDistUPGMA(matrix, min.first.letter, min.second.letter, c.letter);
                matrix.add(new Value(c, newNode, cAndNewNode));
            }
            Iterator<Node> it = currNodes.iterator();
            while (it.hasNext()) {
                Node node = it.next();
                if ((node.letter.equals(min.first.letter)) || (node.letter.equals(min.second.letter))) {
                    it.remove();
                }
            }
            Iterator<Value> it2 = matrix.iterator();
            while (it2.hasNext()) {
                Value value = it2.next();
                if ((value.first.letter.equals(min.first.letter))
                        || (value.first.letter.equals(min.second.letter))
                        || (value.second.letter.equals(min.first.letter))
                        || (value.second.letter.equals(min.second.letter))) {
                    it2.remove();
                }
            }
            currNodes.add(newNode);
        }
        tree.add(newNode);
        System.out.println(inOrderNewick(newNode));
    }

    public void runWPGMA() {
        PriorityQueue<Value> matrix = new PriorityQueue<Value>(new Comparator<Value>() {

            @Override
            public int compare(Value o1, Value o2) {
                if (o1.dist < o2.dist) {
                    return -1;
                } else {
                    if (o1.dist == o2.dist) {
                        return 0;
                    } else {
                        return 1;
                    }
                }
            }

        });
        Set<String> letters = new HashSet<String>(Arrays.asList(new String[]{"A", "B", "C", "D", "E", "F"}));
        HashMap<String, Double> heightsOfNodes = new HashMap<>();
        LinkedList<Node> currNodes = new LinkedList<>();
        for (String letter : letters) {
            heightsOfNodes.put(letter, 0.0);
            currNodes.add(new Node(letter));
        }

        matrix.add(new Value("B", "A", 5));
        matrix.add(new Value("C", "A", 4));
        matrix.add(new Value("D", "A", 7));
        matrix.add(new Value("E", "A", 6));
        matrix.add(new Value("F", "A", 8));
        matrix.add(new Value("C", "B", 7));
        matrix.add(new Value("D", "B", 10));
        matrix.add(new Value("E", "B", 9));
        matrix.add(new Value("F", "B", 11));
        matrix.add(new Value("D", "C", 7));
        matrix.add(new Value("E", "C", 6));
        matrix.add(new Value("F", "C", 8));
        matrix.add(new Value("E", "D", 5));
        matrix.add(new Value("F", "D", 9));
        matrix.add(new Value("F", "E", 8));

        Node newNode = null;
        LinkedList<Node> tree = new LinkedList<>();
        while (!matrix.isEmpty()) {
            Value min = matrix.poll();
            newNode = new Node(min.first, min.second);

            String nameOfNode = min.first.letter + min.second.letter;
            heightsOfNodes.put(nameOfNode, min.dist / 2);
            min.first.height = min.dist / 2 - heightsOfNodes.get(min.first.letter);
            tree.add(min.first);
            min.second.height = min.dist / 2 - heightsOfNodes.get(min.second.letter);
            tree.add(min.second);
            for (Node c : currNodes) {
                if (c.letter.equals(min.first.letter) || c.letter.equals(min.second.letter)) {
                    continue;
                }
                double cAndNewNode = countDistWPGMA(matrix, min.first.letter, min.second.letter, c.letter);
                matrix.add(new Value(c, newNode, cAndNewNode));
            }
            Iterator<Node> it = currNodes.iterator();
            while (it.hasNext()) {
                Node node = it.next();
                if ((node.letter.equals(min.first.letter)) || (node.letter.equals(min.second.letter))) {
                    it.remove();
                }
            }
            Iterator<Value> it2 = matrix.iterator();
            while (it2.hasNext()) {
                Value value = it2.next();
                if ((value.first.letter.equals(min.first.letter))
                        || (value.first.letter.equals(min.second.letter))
                        || (value.second.letter.equals(min.first.letter))
                        || (value.second.letter.equals(min.second.letter))) {
                    it2.remove();
                }
            }
            currNodes.add(newNode);
        }
        tree.add(newNode);
        System.out.println(inOrderNewick(newNode));
    }

    public double countDistWPGMA(PriorityQueue<Value> matrix, String first, String second, String c) {
        double fComp = 0;
        double sComp = 0;
        for (Value value : matrix) {
            if ((value.first.letter.equals(first) && value.second.letter.equals(c))
                    || (value.second.letter.equals(first) && value.first.letter.equals(c))) {
                fComp = value.dist;
            }
            if ((value.first.letter.equals(second) && value.second.letter.equals(c))
                    || (value.second.letter.equals(second) && value.first.letter.equals(c))) {
                sComp = value.dist;
            }
        }
        return (fComp + sComp) / 2;
    }

    public double countDistUPGMA(PriorityQueue<Value> matrix, String first, String second, String c) {
        double fComp = 0;
        int fChildren = 0;
        double sComp = 0;
        int sChildren = 0;
        for (Value value : matrix) {
            if ((value.first.letter.equals(first) && value.second.letter.equals(c))) {
                fComp = value.dist;
                fChildren = value.first.children;
            }
            if (value.second.letter.equals(first) && value.first.letter.equals(c)) {
                fComp = value.dist;
                fChildren = value.second.children;
            }
            if (value.first.letter.equals(second) && value.second.letter.equals(c)) {
                sComp = value.dist;
                sChildren = value.first.children;
            }
            if (value.second.letter.equals(second) && value.first.letter.equals(c)) {
                sComp = value.dist;
                sChildren = value.second.children;
            }
        }
        return (fComp*fChildren + sComp*sChildren) / (fChildren+sChildren);
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

class Value {

//    public LinkedList<String> first = new LinkedList<String>();
//    public LinkedList<String> second = new LinkedList<String>();
    Node first;
    Node second;
    public double dist;

    public Value(String first, String second, double dist) {
        this.first = new Node(first);
        this.second = new Node(second);
        this.dist = dist;
    }

    public Value(Node first, Node second, double dist) {
        this.first = first;
        this.second = second;
        this.dist = dist;
    }

}

class Node {

    public boolean hasChild;
    public Node child1;
    public Node child2;
    public double height;
    public String letter;
    public int children = 1;

    public Node(String letter) {
        this(letter, 0);
    }

    public Node(String letter, double height) {
        this.hasChild = false;
        this.letter = letter;
        this.height = height;
    }

    public Node(Node child1, Node child2, double height) {
        this.hasChild = true;
        this.child1 = child1;
        this.child2 = child2;
        this.height = height;
    }

    public Node(Node child1, Node child2) {
        this.hasChild = true;
        this.child1 = child1;
        this.child2 = child2;
        this.children = child1.children + child2.children;
        this.letter = child1.letter + child2.letter;
//        this.height = height;
    }

    public String getSeq() {
        return letter + ":" + height;
    }
}
