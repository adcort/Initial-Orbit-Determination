#include "hgibbs.h"
#include <cmath>
#include "Matrix.h"
#include "Const.h"

/**
 * @file hgibbs.cpp
 * @brief Implementation of the Herrick-Gibbs approximation for orbit determination.
 *
 * This function implements the Herrick-Gibbs approximation for orbit determination,
 * finding the mean velocity vector for three given position vectors.
 *
 * @param r1 First position vector in IJK coordinates (in meters).
 * @param r2 Second position vector in IJK coordinates (in meters).
 * @param r3 Third position vector in IJK coordinates (in meters).
 * @param Mjd1 Julian date of the first observation (days since 4713 BC).
 * @param Mjd2 Julian date of the second observation (days since 4713 BC).
 * @param Mjd3 Julian date of the third observation (days since 4713 BC).
 * @param v2 Velocity vector at position r2 (in m/s).
 * @param theta Angle between position vectors (in radians).
 * @param error Success indicator ('ok' if successful).
 */

void hgibbs(const std::vector<double>& r1, const std::vector<double>& r2, const std::vector<double>& r3,
            double Mjd1, double Mjd2, double Mjd3,
            std::vector<double>& v2, double& theta, double& theta1, double& copa, std::string& error){
    error = "ok";
    theta = 0.0;
    theta1 = 0.0;
    double magr1 = Matrix::normaVector(r1);
    double magr2 = Matrix::normaVector(r2);
    double magr3 = Matrix::normaVector(r3);

    v2.assign(3, 0.0);

    double tolangle = 0.01745329251994;
    double dt21 = (Mjd2 - Mjd1) * 86400.0;
    double dt31 = (Mjd3 - Mjd1) * 86400.0;
    double dt32 = (Mjd3 - Mjd2) * 86400.0;

    std::vector<double> p = Matrix::productoVectorial(r2,r3);
    std::vector<double> pn = Matrix::unit(p);
    std::vector<double> r1n = Matrix::unit(r1);
    copa = asin(Matrix::productoEscalar(pn, r1n));

    if (fabs(Matrix::productoEscalar(r1n,pn)) > 0.017452406) {
        error = "not coplanar";
    }

    theta = Matrix::angl(r1, r2);
    theta1 = Matrix::angl(r2, r3);

    if (theta > tolangle || theta1 > tolangle){
        error = "angl > 1ø";
    }

    double term1 = -dt32 * (1.0 / (dt21 * dt31) + Const::GM_Earth / (12.0 * magr1 * magr1 * magr1));
    double term2 = (dt32 - dt21) * (1.0 / (dt21 * dt32) + Const::GM_Earth / (12.0 * magr2 * magr2 * magr2));
    double term3 = dt21 * (1.0 / (dt32 * dt31) + Const::GM_Earth / (12.0 * magr3 * magr3 * magr3));

    for (int i = 0; i < 3; i++)
        v2[i] = term1 * r1[i] + term2 * r2[i] + term3 * r3[i];
}