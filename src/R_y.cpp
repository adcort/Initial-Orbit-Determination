#include "R_y.h"
#include <cmath>

/**
 * @file R_y.cpp
 * @brief Implementation of the function to generate a rotation matrix around the Y-axis.
 *
 * This function calculates the 3x3 rotation matrix that rotates a vector
 * in three-dimensional space around the Y-axis by a given angle.
 *
 * @param angle The rotation angle in radians.
 * @return The 3x3 rotation matrix.
 */

Matrix R_y(double angle){
    double C = cos(angle);
    double S = sin(angle);

    //rotmat = zeros(3,3);
    Matrix rotmat(3,3);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;

    return rotmat;
}
