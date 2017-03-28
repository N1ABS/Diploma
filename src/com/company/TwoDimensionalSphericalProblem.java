package com.company;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class TwoDimensionalSphericalProblem {

    PrintWriter pw;
    PrintWriter d;

    /**
     * R - count of iterations by RADIUS (R)
     * F - count of iterations by ANGLE (F)
     * T - count of iterations by TIME (T)
     */
    private int R, F, T;

    private final double EPSILON = 0.1;

    double R_initial = 0.5, F_initial = 0, T_initial = 0;
    double R_endpoint = 1.0, F_endpoint = 2*Math.PI, T_endpoint = 1.0;

    double dR, dF, dT;

    double Cp;

    double lambda;
    double alpha;

    double[][][] Y;
    double[][] TemperatureOfRoom;
    double[][] TemperatureOfAir;
    double[][] Rs;

    ArrayList<double[][][]> Ys = new ArrayList<>();

    public TwoDimensionalSphericalProblem() throws FileNotFoundException {
        pw = new PrintWriter(new File("some.txt"));

        R = 6;
        F = 6;
        T = 10000;

        Cp = 3.11;
        lambda = 0.66;
        alpha = 11.4;

        Y = new double[R + 1][F + 1][T + 1];
        TemperatureOfRoom = new double[R + 1][T + 1];
        TemperatureOfAir = new double[R + 1][T + 1];
        Rs = new double[R + 1][T + 1];

        dR = (R_endpoint - R_initial)/ R;
        dF = (F_endpoint - F_initial)/ F;
        dT = (T_endpoint - T_initial)/ T;
    }

    public void run(){



        for (int t = 0; t <= T; t++) {
            for (int f = 0; f < R; f++) {
                TemperatureOfAir[f][t] = -15;
                TemperatureOfRoom[f][t] = 20;
            }
        }

        for (int t = 0; t < T; t++) {

            for (int r = 0; r <= R; r++) {
                Y[r][0][0] = 15;
                Y[r][1][0] = 15;
                Y[r][F][0] = 15;
                Y[r][F - 1][0] = 15;
            }

            for (int f = 0; f <= F; f++) {
                Y[0][f][t+1] = (lambda*Math.pow(R_initial, 2)*Y[1][f][t+1] + alpha*dR*TemperatureOfRoom[f][t+1])
                        /(lambda*Math.pow(R_initial, 2) + alpha*dR);
                Y[R][f][t+1] = (lambda*Math.pow(R_endpoint, 2)*Y[R-1][f][t+1] + alpha*dR*TemperatureOfAir[f][t+1])
                        /(lambda*Math.pow(R_endpoint, 2) + alpha*dR);
            }


            for (int r = 1; r < R; r++) {
                double Pi = (lambda * dT) / (Cp * Math.pow(R_initial+r*dR, 2) * dR * dF);
                System.out.println("Kourant on ["+t+"]["+r+"] = " + Pi);

                for (int f = 1; f < F; f++) {
                    Y[r][f][t + 1] = (Pi *
                            ((Math.pow(R_initial+r*dR+dR/2, 2)*dF/dR)*Y[r+1][f][t+1] -
                                    (Math.pow(R_initial+r*dR-dR/2, 2)*dF/dR)*Y[r-1][f][t+1])+
                            (((1+1/Math.tan(F_initial+f*dF+dF/2))*dR/dF)*Y[r][f+1][t+1] -
                                    ((1+1/Math.tan(F_initial+f*dF-dF/2))*dR/dF)*Y[r][f-1][t+1])
                    )/(
                            1 +
                            Pi*(Math.pow(R_initial+r*dR+dR/2, 2)-Math.pow(R_initial+r*dR-dR/2, 2))*dF/dR +
                            Pi*(1/Math.tan(F_initial+f*dF+dF/2) - 1/Math.tan(F_initial+f*dF-dF/2))*dR/dF
                    ) + Y[r][f][t];
                }
            }


        }

        for (int r = 0; r <= R; r++) {
            for (int f = 0; f <= F; f++) {
                if (f == 9) System.out.print("\t");
                System.out.print(Y[r][f][1] + "\t");
            }
            System.out.println();
        }

    }

    public void alternateRun() throws FileNotFoundException {

        for (int t = 0; t <= T; t++) {
            for (int f = 0; f < R; f++) {
                TemperatureOfAir[f][t] = -15;
                TemperatureOfRoom[f][t] = 20;
            }
        }

        double Arf[][] = new double[R+1][F+1];
        double Crf[][] = new double[R+1][F+1];
        double Drf[][] = new double[R+1][F+1];
        double Erf[][] = new double[R+1][F+1];
        double Krf[][] = new double[R+1][F+1];


        System.out.println(Math.sin(Math.toRadians(45)) + "  " + Math.sin(Math.PI/4));
        System.out.println(Math.sin(Math.toRadians(270)) + "  " + Math.sin(3*Math.PI/2));
        System.out.println(Math.sin(Math.toRadians(30)) + "  " + Math.sin(Math.PI/6));
        System.out.println(Math.sin(Math.toRadians(90))  + "  " + Math.sin(Math.PI/2));
        System.out.println(Math.sin(Math.toRadians(0))   + "  " + Math.sin(0*Math.PI/180));
        System.out.println();
        System.out.println();

        for (int r = 0; r <= R; r++) {
            for (int f = 0; f <= F; f++) {
                Arf[r][f] = Cp*Math.pow(R_initial+r*dR, 2)*Math.sin(Math.toRadians(F_initial+f*dF));
                Crf[r][f] = Math.sin(Math.toRadians(F_initial+f*dF))*Math.pow(R_initial+r*dR+dR/2, 2)*lambda*dT/Math.pow(dR,2);
                Drf[r][f] = Math.sin(Math.toRadians(F_initial+f*dF))*Math.pow(R_initial+r*dR-dR/2, 2)*lambda*dT/Math.pow(dR,2);
                Erf[r][f] = Math.sin(Math.toRadians(F_initial+f*dF+dF/2))*dT/Math.pow(dF,2);
                Krf[r][f] = Math.sin(Math.toRadians(F_initial+f*dF-dF/2))*dT/Math.pow(dF,2);

                //System.out.println(Math.sin((F_initial+f*dF)));
                //System.out.println(F_initial);
                //System.out.println(F_endpoint);
                System.out.print(Arf[r][f]+" ");
            }
            System.out.println();
        }
        System.out.println();


        int fin_t = 0;
        for (int t = 0; ; t++) {

            double MAX_R = Rs[0][0];


            for (int f = 0; f <= F; f++) {
                Y[0][f][t+1] = (lambda*Math.pow(R_initial, 2)*Y[1][f][t+1] + alpha*dR*TemperatureOfRoom[f][t+1])
                        /(lambda*Math.pow(R_initial, 2) + alpha*dR);
                Y[R][f][t+1] = (lambda*Math.pow(R_endpoint, 2)*Y[R-1][f][t+1] + alpha*dR*TemperatureOfAir[f][t+1])
                        /(lambda*Math.pow(R_endpoint, 2) + alpha*dR);
            }

            for (int r = 1; r < R; r++) {
                for (int f = 1; f < F; f++) {
                    Y[r][f][t+1] = (Crf[r][f]*Y[r+1][f][t+1] + Drf[r][f]*Y[r-1][f][t+1] +
                                    Erf[r][f]*Y[r][f+1][t+1] + Krf[r][f]*Y[r][f-1][t+1] + Arf[r][f]*Y[r][f][t])
                                   /
                                   (Arf[r][f] + Crf[r][f] + Drf[r][f] + Erf[r][f] + Krf[r][f]);

                    Rs[r][f] = (Math.abs(Y[r][f][t+1] - Y[r][f][t]));

                    if (MAX_R < Rs[r][f]) MAX_R = Rs[r][f];
                }
            }

            /*System.out.println("MAX_R = " + MAX_R);
            pw.print("MAX_R["+t+"] = ");
            pw.println((double)MAX_R);
            pw.flush();
            System.out.println("t = " + t);
            */
            if (MAX_R < EPSILON) {
                fin_t = t;
                break;
            }

        }


        pw.close();

        for (int t = fin_t-2; t < fin_t-1; t++) {
            System.out.println("Time " + t);
            for (int r = 0; r <= R; r++) {
                for (int f = 0; f <= F; f++) {
                    if (f == 9) System.out.print("\t");
                    System.out.print(Y[r][f][t] + "\t");
                }
                System.out.println();
            }

            System.out.println();
        }

    }


    public static void main(String[] args) throws FileNotFoundException {
        TwoDimensionalSphericalProblem pr = new TwoDimensionalSphericalProblem();
        pr.alternateRun();
    }

}
