package com.company;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

/**
 * Created by Sequoia on 3/10/2017.
 */
public class FileOutput {

    public static void print(double number) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new File("some.txt"));
        pw.print((double) number + " ");
    }

    public static void println(double number) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new File("some.txt"));
        pw.println((double) number);
    }

    public static void println() throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new File("some.txt"));
        pw.println();
    }
}
