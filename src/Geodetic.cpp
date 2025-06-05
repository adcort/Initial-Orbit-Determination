#include "Geodetic.h"
#include "Const.h"
#include "Matrix.h"
#include <limits>
#include <cmath>
#include <stdexcept>

/**
 * @file Geodetic.cpp
 * @brief Converts a given position vector to geodetic coordinates.
 *
 * @details This function calculates geodetic coordinates (longitude, latitude, and altitude)
 * from a given position vector in Cartesian coordinates.
 *
 * @param r Vector of Cartesian coordinates [x, y, z] in meters.
 * @param lon Geodetic longitude in radians.
 * @param lat Geodetic latitude in radians.
 * @param h Height above the ellipsoid in meters.
 *
 * @return void
 *
 * @note Last modified: 2015/08/12 by M. Mahooti
 */

void Geodetic(double& lon, double& lat, double& h, const std::vector<double>& r) {
    //Variables locales
    const double eps = std::numeric_limits<double>::epsilon();
    const double R_equ = Const::R_Earth;
    const double f = Const::f_Earth;
    const double epsRequ = eps * R_equ; // Convergence criterion
    const double e2 = f * (2.0 - f); // Square of eccentricity

    double X = r[0]; // Cartesian coordinates
    double Y = r[1];
    double Z = r[2];
    double rho2 = X * X + Y * Y; // Square of distance from z-axis

    // Check validity of input data
    // Calculamos la norma dentro del if
    if (Matrix::normaVector(r) == 0.0) {
        lon = 0.0;
        lat = 0.0;
        h = -Const::R_Earth;
        throw std::runtime_error("invalid input in Geodetic constructor\n");
    }

    // Iteración
    double dZ = e2 * Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;
    while (true) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        SinPhi = ZdZ / Nh; // Sine of geodetic latitude
        N = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        dZ_new = N * e2 * SinPhi;
        if (abs(dZ - dZ_new) < epsRequ) {
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2(Y, X);
    lat = atan2(ZdZ, sqrt(rho2));
    h = Nh - N;
}