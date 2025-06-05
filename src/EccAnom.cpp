#include "EccAnom.h"
#include "Const.h"
#include <cmath>
#include <limits>
#include <stdexcept>

/**
 * @file EccAnom.cpp
 * @brief Computes the eccentric anomaly for elliptic orbits.
 *
 * This function computes the eccentric anomaly for elliptic orbits based on the given mean anomaly
 * and eccentricity.
 *
 * @param M Mean anomaly in radians.
 * @param e Eccentricity of the orbit (must be in the range [0, 1]).
 *
 * @return Eccentric anomaly in radians.
 *
 * @note Last modified: 2015/08/12 by M. Mahooti.
 */

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;

    // Starting value
    M = fmod(M, Const::doublepi);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = Const::pi;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // Iteration
    //Aquí usamos en C++ el equivalente a eps de Matlab que es epsilon de la librería numeric_limits.
    //Devuelve el menor entero positivo que podemos sumar a 1.0 y que sea diferente a 1.0. Se usa para
    //calcular criterios de convergencia y ver que una diferencia es menor que un epsilon (arbitrario) dado.
    while (fabs(f) > 1e2 * std::numeric_limits<double>::epsilon()) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i = i + 1;
        if (i == maxit) {
            throw std::runtime_error("convergence problems in EccAnom");
        }
    }

    return E;
}

