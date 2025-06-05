#include "MeasUpdate.h"

/**
 * @file MeasUpdate.cpp
 * @brief Performs a measurement update using the Kalman filter.
 *
 * This function updates the state and covariance using the Extended Kalman Filter.
 * It calculates the Kalman gain, updates the estimated state, and updates the state
 * covariance based on the received measurement and predicted measurement model.
 *
 * @param x     Estimated state vector before and updated state vector after the update.
 * @param z     Measurement vector.
 * @param g     Predicted measurement vector from the model.
 * @param s     Vector containing the standard deviations of measurements.
 * @param G     Jacobian matrix of the measurement function.
 * @param P     State covariance matrix before and after the update.
 * @param n     Dimension of the state vector.
 * @param K     Kalman gain calculated during the update.
 *
 * @note This function assumes that input matrices and vectors are correctly sized
 *       and that matrix inversion and multiplication operations are valid for the matrices involved.
 */

void MeasUpdate(Matrix& x,  Matrix& z,  Matrix& g,  std::vector<double>& s,  Matrix& G, Matrix& P, int n, Matrix& K) {
    // Dimensión de la medición
    int m = z.getNumFilas();

    // Inicializar matriz de covarianza inversa de la medición
    Matrix Inv_W(m, m);
    for (int i = 0; i < m; ++i) {
        Inv_W.set(i, i, s[i] * s[i]);
    }

    // Calcular la ganancia de Kalman
    Matrix Gt = G.transponer();
    Matrix temp = Inv_W + G * P * Gt;
    Matrix tempInv = temp.inversa();
    K = P * Gt * tempInv;

    // Actualizar el estado
    x = x + K * (z - g);

    // Actualizar la covarianza del estado
    Matrix I = Matrix::identidad(n);
    P = (I - K * G) * P;
}