package jacobi;

import com.sun.org.apache.xml.internal.security.Init;

/**
 * Created by Sequoia on 3/11/2017.
 */
public class Sphere {
    private int R, F, T;
    double dR, dF, dT;

    Coefficients coef;
    InitialCoordinates start;
    FinalCoordinates end;

    double[/*R*/][/*F*/][/*T*/] Yprev;
    double[/*R*/][/*F*/][/*T*/] Ynext;
    double[/*F*/][/*T*/] TemperatureOfRoom;
    double[/*F*/][/*T*/] TemperatureOfAir;

    double[/*F*/][/*T*/] Rs;

    double A[][] = new double[R+1][F+1];
    double C[][] = new double[R+1][F+1];
    double D[][] = new double[R+1][F+1];
    double E[][] = new double[R+1][F+1];
    double K[][] = new double[R+1][F+1];

    final double EPSILON = 0.1;
    int fin_t = 0;

    public Sphere(Coefficients coef, InitialCoordinates start, FinalCoordinates end){
        this.coef = coef;
        this.start = start;
        this.end = end;

        Yprev = new double [R + 1][F + 1][T + 1];
        Ynext = new double [R + 1][F + 1][T + 1];
        TemperatureOfRoom = new double [R + 1][T + 1];
        TemperatureOfAir = new double [R + 1][T + 1];
        Rs = new double [F + 1][T + 1];
    }

    public void setIterations(int R, int F, int T){
        this.R = R;
        this.F = F;
        this.T = T;
    }

    public void setAmbientTemperature() {
        for (int t = 0; t <= T; t++) {
            for (int f = 0; f <= F; f++) {
                TemperatureOfRoom[f][t] = 18;
                TemperatureOfAir[f][t] = -15;
            }
        }
    }

    private void setNumericalCoefficients() {
        for (int r = 0; r <= R; r++) {
            for (int f = 0; f <= F; f++) {
                A[r][f] = coef.Cp*Math.pow(start.r+r*dR, 2)*Math.sin(Math.toRadians(start.f+f*dF));
                C[r][f] = Math.sin(Math.toRadians(start.f+f*dF))*Math.pow(start.r+r*dR+dR/2, 2)*coef.lambda*dT/Math.pow(dR,2);
                D[r][f] = Math.sin(Math.toRadians(start.f+f*dF))*Math.pow(start.r+r*dR-dR/2, 2)*coef.lambda*dT/Math.pow(dR,2);
                E[r][f] = Math.sin(Math.toRadians(start.f+f*dF+dF/2))*dT/Math.pow(dF,2);
                K[r][f] = Math.sin(Math.toRadians(start.f+f*dF-dF/2))*dT/Math.pow(dF,2);
            }
        }
    }

    public void init() {
        dR = (end.r - start.r) / R;
        dF = (end.f - start.f) / F;
        dR = (end.t - start.t) / T;
    }

    public void run() {
        setNumericalCoefficients();

        for (int t = 0; t < T; t++) {

            while (true) {


                for (int f = 0; f <= F; f++) {
                    Ynext[0][f][t + 1] = (coef.lambda * Math.pow(start.r, 2) * Yprev[1][f][t + 1] + coef.alpha * dR * TemperatureOfRoom[f][t + 1])
                            / (coef.lambda * Math.pow(start.r, 2) + coef.alpha * dR);
                    Ynext[R][f][t + 1] = (coef.lambda * Math.pow(end.r, 2) * Yprev[R - 1][f][t + 1] + coef.alpha * dR * TemperatureOfAir[f][t + 1])
                            / (coef.lambda * Math.pow(end.r, 2) + coef.alpha * dR);
                }

                double MAX_R = 0;

                for (int r = 1; r < R; r++) {
                    for (int f = 1; f < F; f++) {
                        Ynext[r][f][t + 1] = (C[r][f] * Yprev[r + 1][f][t + 1] + D[r][f] * Yprev[r - 1][f][t + 1] +
                                E[r][f] * Yprev[r][f + 1][t + 1] + K[r][f] * Yprev[r][f - 1][t + 1] + A[r][f] * Yprev[r][f][t])
                                /
                                (A[r][f] + C[r][f] + D[r][f] + E[r][f] + K[r][f]);

                        Rs[r][f] = (Math.abs(Ynext[r][f][t + 1] - Yprev[r][f][t]));

                        MAX_R = (MAX_R < Rs[r][f]) ? Rs[r][f] : MAX_R;
                    }
                }

                if (MAX_R < EPSILON) {
                    fin_t = t;
                    break;
                }

                Yprev = Ynext;

            }

        }
    }
}
