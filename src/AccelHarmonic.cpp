#include "AccelHarmonic.h"
#include "Const.h"
#include "Legendre.h"
#include <cmath>

/**
 * @file AccelHarmonic.cpp
 * @brief Computes the acceleration due to the harmonic gravity field of the central body.
 *
 * This function calculates the acceleration experienced by a satellite due to the
 * gravitational field of the central body represented in terms of spherical harmonics.
 *
 * @param r     Satellite position vector in the inertial system.
 * @param E     Transformation matrix to body-fixed system.
 * @param n_max Maximum degree of spherical harmonics.
 * @param m_max Maximum order of spherical harmonics (must satisfy m_max <= n_max; m_max = 0 for zonals).
 * @param a     Output parameter: Acceleration vector (a = d^2r/dt^2).
 * @param GGM03S Gravity model coefficients matrix (should contain Cnm and Snm coefficients).
 *
 * @return void
 *
 * @note The function assumes that r and E are properly defined Matrix objects.
 * @note GGM03S should contain the coefficients Cnm and Snm required for the calculation.
 * @note Last modified: 2015/08/12 by M. Mahooti.
 */


void AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max, Matrix& a, Matrix& GGM03S) {

    const float r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    const float gm = 398600.4415e9;    // [m^3/s^2]; GGM03S

    // Body-fixed position
    Matrix r_bf = E * r;

    // Matrices para almacenar Cnm y Snm
    Matrix Cnm(n_max + 1, m_max + 1);
    Matrix Snm(n_max + 1, m_max + 1);
    // Llenar Cnm y Snm desde GGM03S
    for (int n = 0; n <= n_max; ++n) {
        for (int m = 0; m <= std::min(m_max, n); ++m) {
            Cnm.set(n, m, GGM03S.get(n,2));  // El tercer elemento de cada fila para Cnm
            Snm.set(n, m, GGM03S.get(n,3));  // El cuarto elemento de cada fila para Snm
        }
    }


    // Auxiliary quantities
    float d = r_bf.normaVector();                     // distance
    float latgc = std::asin(r_bf.get(0, 2) / d);
    float lon = std::atan2(r_bf.get(0, 1), r_bf.get(0, 0));


    Matrix pnm(n_max + 1, m_max + 1);
    Matrix dpnm(n_max + 1, m_max + 1);

    Legendre(n_max, m_max, latgc, pnm, dpnm);

    float dUdr = 0;
    float dUdlatgc = 0;
    float dUdlon = 0;
    float q1 = 0, q2 = 0, q3 = 0;

    // Calcular componentes de la aceleración
    for (int n = 0; n <= n_max; ++n) {
        float b1 = (-gm / (d * d)) * std::pow(r_ref / d, n) * (n + 1);
        float b2 = (gm / d) * std::pow(r_ref / d, n);
        float b3 = (gm / d) * std::pow(r_ref / d, n);

        for (int m = 0; m <= std::min(m_max, n); ++m) {
            q1 += pnm.get(n, m) * (Cnm.get(n, m) * std::cos(m * lon) + Snm.get(n, m) * std::sin(m * lon));
            q2 += dpnm.get(n, m) * (Cnm.get(n, m) * std::cos(m * lon) + Snm.get(n, m) * std::sin(m * lon));
            q3 += m * pnm.get(n, m) * (Snm.get(n, m) * std::cos(m * lon) - Cnm.get(n, m) * std::sin(m * lon));
        }

        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;
        q1 = q2 = q3 = 0;
    }

    // Body-fixed acceleration
    float r2xy = r_bf.get(0, 0) * r_bf.get(0, 0) + r_bf.get(0, 1) * r_bf.get(0, 1);
    float ax = (1 / d * dUdr - r_bf.get(0, 2) / (d * d * std::sqrt(r2xy)) * dUdlatgc) * r_bf.get(0, 0) - (1 / r2xy * dUdlon) * r_bf.get(0, 1);
    float ay = (1 / d * dUdr - r_bf.get(0, 2) / (d * d * std::sqrt(r2xy)) * dUdlatgc) * r_bf.get(0, 1) + (1 / r2xy * dUdlon) * r_bf.get(0, 0);
    float az = 1 / d * dUdr * r_bf.get(0, 2) + std::sqrt(r2xy) / (d * d) * dUdlatgc;
    double a_bfv[] = { ax, ay, az };
    Matrix a_bf(3,1, a_bfv, 3);

    // Inertial acceleration
    a = E.transponer() * a_bf;
}