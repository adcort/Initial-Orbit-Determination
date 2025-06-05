#include "doubler.h"
#include "Const.h"
#include "Matrix.h"
#include <cmath>
#include <iostream>

/**
 * @file doubler.cpp
 * @brief Calculates various parameters related to orbits and geometry for a double-range system.
 *
 * This function computes multiple parameters such as position vectors, norms, angles, and more,
 * based on various input parameters and provided matrices.
 *
 * @param cc1 Correction constant 1.
 * @param cc2 Correction constant 2.
 * @param magrsite1 Magnitude of site 1.
 * @param magrsite2 Magnitude of site 2.
 * @param magr1in Magnitude of range 1.
 * @param magr2in Magnitude of range 2.
 * @param los1 Line of sight state vector 1.
 * @param los2 Line of sight state vector 2.
 * @param los3 Line of sight state vector 3.
 * @param rsite1 Position matrix of site 1.
 * @param rsite2 Position matrix of site 2.
 * @param rsite3 Position matrix of site 3.
 * @param t1 Time 1.
 * @param t3 Time 3.
 * @param direct Character indicating the direction.
 * @param r2 Resulting position matrix 2 (output).
 * @param r3 Resulting position matrix 3 (output).
 * @param f1 Resulting value 1.
 * @param f2 Resulting value 2.
 * @param q1 Resulting value 3.
 * @param magr1 Magnitude of range 1 (output).
 * @param magr2 Magnitude of range 2 (output).
 * @param a Calculated orbital parameter (output).
 * @param deltae32 Delta Epsilon 32 (output).
 */

void doubler(double cc1, double cc2, double& magrsite1, double& magrsite2, double magr1in, double magr2in,
             Matrix& los1, Matrix& los2, Matrix& los3, Matrix& rsite1, Matrix& rsite2, Matrix& rsite3,
             double t1, double t3, char direct, Matrix& r2, Matrix& r3, double& f1, double& f2, double& q1,
             double& magr1, double& magr2, double& a, double& deltae32) {

    double rho1 = (-cc1 + sqrt(cc1*cc1 - 4*(magrsite1*magrsite1 - magr1in*magr1in))) / 2.0;
    double rho2 = (-cc2 + sqrt(cc2*cc2 - 4*(magrsite2*magrsite2 - magr2in*magr2in))) / 2.0;

    Matrix r1 = los1 * rho1 + rsite1;
    r2 = los2 * rho2 + rsite2;

    magr1 = r1.normaVector();  // Norma del vector r1
    magr2 = r2.normaVector();  // Norma del vector r2

    Matrix w(1,3);

    if (direct == 'y') {
        w = (r1.productoVectorial(r2)) / (magr1 * magr2);
    } else {
        w = r1.productoVectorial(r2) / (-magr1 * magr2);
    }

    double rho3 = -rsite3.productoEscalar(w) / los3.productoEscalar(w);
    r3 = los3 * rho3 + rsite3;
    double magr3 = r3.normaVector();  // Norma del vector r3

    double cosdv21 = r2.productoEscalar(r1) / (magr2 * magr1);
    double sindv21 = (r2.productoVectorial(r1)).normaVector() / (magr2 * magr1);
    double dv21 = atan2(sindv21, cosdv21);

    double cosdv31 = r3.productoEscalar(r1) / (magr3 * magr1);
    double sindv31 = sqrt(1.0 - cosdv31 * cosdv31);
    double dv31 = atan2(sindv31, cosdv31);

    double cosdv32 = r3.productoEscalar(r2) / (magr3 * magr2);
    double sindv32 = ((r3.productoVectorial(r2)).normaVector()) / (magr3 * magr2);
    double dv32 = atan2(sindv32, cosdv32);

    double p, e, s, c;

    if (dv31 > M_PI) {
        double c1 = (magr2 * sindv32) / (magr1 * sindv31);
        double c3 = (magr2 * sindv21) / (magr3 * sindv31);
        p = (c1 * magr1 + c3 * magr3 - magr2) / (c1 + c3 - 1);
    } else {
        double c1 = (magr1 * sindv31) / (magr2 * sindv32);
        double c3 = (magr1 * sindv21) / (magr3 * sindv32);
        p = (c3 * magr3 - c1 * magr2 + magr1) / (-c1 + c3 + 1);
    }

    double ecosv1 = p / magr1 - 1;
    double ecosv2 = p / magr2 - 1;
    double ecosv3 = p / magr3 - 1;
    double esinv2;
    if (dv21 != M_PI) {
        esinv2 = (-cosdv21 * ecosv2 + ecosv1) / sindv21;
    } else {
        esinv2 = (cosdv32 * ecosv2 - ecosv3) / sindv31;
    }
    e = sqrt(ecosv2*ecosv2+esinv2*esinv2);
    a = p/(1-e*e);
    double deltam12=0.0; double deltam32 = 0.0; double n=0.0;
    if (e * e < 0.99) {
        n = sqrt(Const::GM_Earth / (a * a * a));
        s = magr2 / p * sqrt(1 - e * e) * esinv2;
        c = magr2 / p * (e * e + ecosv2);

        double sinde32 = magr3 / sqrt(a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s;
        double cosde32 = 1 - magr2 * magr3 / (a * p) * (1 - cosdv32);
        deltae32 = atan2(sinde32, cosde32);

        double sinde21 = magr1 / sqrt(a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s;
        double cosde21 = 1 - magr2 * magr1 / (a * p) * (1 - cosdv21);
        double deltae21 = atan2(sinde21, cosde21);

        deltam32 = deltae32 + 2 * s * sin(deltae32 / 2)*sin(deltae32 / 2) - c * sin(deltae32);
        deltam12 = -deltae21 + 2 * s * sin(deltae21 / 2)*sin(deltae21 / 2)  + c * sin(deltae21);

    } else {
        n = sqrt(Const::GM_Earth / -(a * a * a));

        s = magr2 / p * sqrt(e * e - 1) * esinv2;
        c = magr2 / p * (e * e + ecosv2);

        double sindh32 = magr3 / sqrt(-a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s;
        double sindh21 = magr1 / sqrt(-a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s;

        double deltah32 = log(sindh32 + sqrt(sindh32 * sindh32 + 1));
        double deltah21 = log(sindh21 + sqrt(sindh21 * sindh21 + 1));

        deltam32 = -deltah32 + 2 * s * sinh(deltah32 / 2)*sinh(deltah32 / 2) + c * sinh(deltah32);
        deltam12 = deltah21 + 2 * s * sinh(deltah21 / 2)*sinh(deltah21 / 2) - c * sinh(deltah21);

        deltae32 = deltah32;
    }

    f1 = t1 - deltam12 / n;
    f2 = t3 - deltam32 / n;

    q1 = sqrt(f1 * f1 + f2 * f2);

}