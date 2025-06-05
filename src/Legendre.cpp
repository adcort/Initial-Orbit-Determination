#include "Legendre.h"
#include <cmath>

/**
 * @file Legendre.cpp
 * @brief Computes the coefficient matrices of Legendre polynomials and their derivatives.
 *
 * This function computes the coefficient matrices of Legendre polynomials `pnm`
 * and their derivatives `dpnm` up to a given order `n` and degree `m`, evaluated at angle `fi`.
 *
 * @param n Maximum order of the Legendre polynomials.
 * @param m Maximum degree of the Legendre polynomials.
 * @param fi Angle in radians at which the Legendre polynomials are evaluated.
 * @param pnm Reference to a matrix to store the coefficients of Legendre polynomials.
 * @param dpnm Reference to a matrix to store the coefficients of derivatives of Legendre polynomials.
 */

void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm) {


    // Asignar valores iniciales
    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2, 2) = std::sqrt(3) * std::cos(fi);
    dpnm(2, 2) = -std::sqrt(3) * std::sin(fi);

    // Diagonal coefficients
    for (int i = 2; i <= n; ++i) {
        pnm(i + 1, i + 1) = std::sqrt((2 * i + 1) / (2.0 * i)) * std::cos(fi) * pnm(i, i);
    }
    for (int i = 2; i <= n; ++i) {
        dpnm(i + 1, i + 1) = std::sqrt((2 * i + 1) / (2.0 * i)) * (std::cos(fi) * dpnm(i, i) - std::sin(fi) * pnm(i, i));
    }

    // Horizontal first step coefficients
    for (int i = 1; i <= n; ++i) {
        pnm(i + 1, i) = std::sqrt(2 * i + 1) * std::sin(fi) * pnm(i, i);
    }
    for (int i = 1; i <= n; ++i) {
        dpnm(i + 1, i) = std::sqrt(2 * i + 1) * (std::cos(fi) * pnm(i, i) + std::sin(fi) * dpnm(i, i));
    }

    // Horizontal second step coefficients
    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; ++i) {
            pnm(i + 1, j + 1) = std::sqrt((2 * i + 1) / ((i - j) * (i + j))) * (std::sqrt(2 * i - 1) * std::sin(fi) * pnm(i, j + 1) -
                                                                                std::sqrt(((i + j - 1) * (i - j - 1)) / (2.0 * i - 3)) * pnm(i - 1, j + 1));
        }
        j = j + 1;
        k = k + 1;
        if (j > m) break;
    }

    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; ++i) {
            dpnm(i + 1, j + 1) = std::sqrt((2 * i + 1) / ((i - j) * (i + j))) * (std::sqrt(2 * i - 1) * std::sin(fi) * dpnm(i, j + 1) +
                                                                                 std::sqrt(2 * i - 1) * std::cos(fi) * pnm(i, j + 1) - std::sqrt(((i + j - 1) * (i - j - 1)) / (2.0 * i - 3)) * dpnm(i - 1, j + 1));
        }
        j = j + 1;
        k = k + 1;
        if (j > m) break;
    }
}