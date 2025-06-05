#include "elements.h"
#include <cmath>
#include "Const.h"
#include "Matrix.h"
/**
 * @file elements.cpp
 * @brief Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits.
 *
 * @details This function computes the osculating Keplerian elements (semilatus rectum, semimajor axis, eccentricity,
 * inclination, longitude of the ascending node, argument of pericenter, and mean anomaly) from the satellite state
 * vector for elliptic orbits.
 *
 * @param y State vector (x, y, z, vx, vy, vz) in meters and meters per second.
 *
 * @param p Semilatus rectum in meters.
 * @param a Semimajor axis in meters.
 * @param e Eccentricity.
 * @param i Inclination in radians.
 * @param Omega Longitude of the ascending node in radians.
 * @param omega Argument of pericenter in radians.
 * @param M Mean anomaly in radians.
 *
 * @note The function cannot be used with state vectors describing a circular or non-inclined orbit.
 *
 * @return void
 *
 * @note Last modified: 2015/08/12 by M. Mahooti
 */
void elements(const std::vector<double>& y, double& p, double& a, double& e, double& i, double& Omega, double& omega, double& M) {
    const double pi2 = 2 * M_PI;

    std::vector<double> r(y.begin(), y.begin() + 3);  // Position
    std::vector<double> v(y.begin() + 3, y.end());    // Velocity

    std::vector<double> h = Matrix::productoVectorial(r, v);              // Areal velocity
    double magh = Matrix::normaVector(h);
    p = magh * magh / Const::GM_Earth;
    double H = Matrix::normaVector(h);

    Omega = atan2(h[0], -h[1]);                       // Long. ascend. node
    Omega = fmod(Omega, pi2);
    i = atan2(sqrt(h[0] * h[0] + h[1] * h[1]), h[2]); // Inclination
    double u = atan2(r[2] * H, -r[0] * h[1] + r[1] * h[0]);  // Arg. of latitude

    double R = Matrix::normaVector(r);                               // Distance
    a = 1 / (2 / R - Matrix::productoEscalar(v, v) / Const::GM_Earth);   // Semi-major axis

    double eCosE = 1 - R / a;                         // e*cos(E)
    double eSinE = Matrix::productoEscalar(r, v) / sqrt(Const::GM_Earth * a); // e*sin(E)

    double e2 = eCosE * eCosE + eSinE * eSinE;
    e = sqrt(e2);                                     // Eccentricity
    double E = atan2(eSinE, eCosE);                   // Eccentric anomaly

    M = fmod(E - eSinE, pi2);                         // Mean anomaly
    double nu = atan2(sqrt(1.0 - e2) * eSinE, eCosE - e2); // True anomaly

    omega = fmod(u - nu, pi2);                        // Arg. of perihelion
}