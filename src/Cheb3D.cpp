#include "Cheb3D.h"
#include <stdexcept>

/**
 * @file Cheb3D.cpp
 * @brief Chebyshev approximation of 3-dimensional vectors.
 *
 * This function computes the Chebyshev approximation for 3-dimensional vectors.
 *
 * @param t  The time at which the approximation is evaluated.
 * @param N  Number of coefficients.
 * @param Ta Begin interval.
 * @param Tb End interval.
 * @param Cx Coefficients of the Chebyshev polynomial (x-coordinate).
 * @param Cy Coefficients of the Chebyshev polynomial (y-coordinate).
 * @param Cz Coefficients of the Chebyshev polynomial (z-coordinate).
 *
 * @return An array containing the approximated 3-dimensional vector.
 *
 * @throws std::out_of_range If the time t is out of the interval [Ta, Tb].
 *
 * @note Last modified: 2018/01/27 by M. Mahooti
 */

std::vector<double> Cheb3D(double t, int N, double Ta, double Tb, const std::vector<double>& Cx, const std::vector<double>& Cy, const std::vector<double>& Cz) {
    // Check validity
    if (t < Ta || t > Tb) {
        throw std::runtime_error("ERROR: Time out of range in Cheb3D::Value\n"); //Lanzamos excepción
    }

    // Clenshaw algorithm
    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    //Hacemos esto para sustituir la función zeros de Matlab
    //que rellena una fila de 3 entradas con ceros
    std::vector<double> f1 = {0.0, 0.0, 0.0};
    std::vector<double> f2 = {0.0, 0.0, 0.0};

    for (int i = N; i >= 2; --i) {
        std::vector<double> old_f1 = f1;
        //Aquí hacemos el producto para cada una de las componentes ya que solo son 3
        f1[0] = 2 * tau * f1[0] - f2[0] + Cx[i];
        f1[1] = 2 * tau * f1[1] - f2[1] + Cy[i];
        f1[2] = 2 * tau * f1[2] - f2[2] + Cz[i];
        f2 = old_f1;
    }

    std::vector<double> ChebApp(3);
    //Repetimos el método anterior para cada entrada
    ChebApp[0] = tau * f1[0] - f2[0] + Cx[0];
    ChebApp[1] = tau * f1[1] - f2[1] + Cy[0];
    ChebApp[2] = tau * f1[2] - f2[2] + Cz[0];

    return ChebApp;
}