#include "PoleMatrix.h"
#include "R_y.h"
#include "R_x.h"
/**
 * @file PoleMatrix.cpp
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date.
 *
 * This function computes the pole matrix based on the pole coordinates (xp, yp).
 *
 * @param xp Pole coordinate x
 * @param yp Pole coordinate y
 * @return PoleMat Pole matrix
 *
 * @note Last modified: 2015/08/12 by M. Mahooti
 */

Matrix PoleMatrix(double xp, double yp){

    Matrix PoleMat = R_y(-xp) * R_x(-yp);
    return PoleMat;
}