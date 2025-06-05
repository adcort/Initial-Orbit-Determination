#include <cstdio>
#include "TimeUpdate.h"

/**
 * @file TimeUpdate.cpp
 * @brief Performs a time update on the state covariance matrix.
 *
 * This function performs a time update on the state covariance matrix `P`
 * using the state transition matrix `Phi` and a constant `Qdt` representing
 * the process noise. If only `P` and `Phi` are provided, `Qdt` defaults to `0.0`.
 *
 * @param P Reference to the state covariance matrix.
 * @param Phi Reference to the state transition matrix.
 * @param Qdt Constant representing the process noise (default value is 0.0).
 * @return The updated state covariance matrix.
 */

Matrix TimeUpdate(Matrix &P, Matrix &Phi, double Qdt){
    //If nargin==2, Qdt tomará por defecto el valor 0.0.
    return Phi*P*Phi.transponer() + Qdt;
}