package com.company;

import jacobi.Coefficients;
import jacobi.FinalCoordinates;
import jacobi.InitialCoordinates;
import jacobi.Sphere;

/**
 * Created by Sequoia on 3/11/2017.
 */
public class Main {

    public static void main(String[] args) {
        Coefficients coef = new Coefficients(3.11, 0.66,11.4);
        InitialCoordinates start = new InitialCoordinates(0.5, 0, 0);
        FinalCoordinates end = new FinalCoordinates(2, 2*Math.PI, 30);

        Sphere sphere = new Sphere(coef, start, end);
        sphere.setIterations(10,10,10);
        sphere.setAmbientTemperature();
        sphere.init();
        sphere.run();
    }
}
