#include "NutMatrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"
/**
 * @file NutMatrix.cpp
 * @brief Transformation from mean to true equator and equinox.
 *
 * This function computes the nutation matrix for a given Modified Julian Date (Terrestrial Time).
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return NutMat Nutation matrix
 *
 * @note Last modified: 2015/08/12 by M. Mahooti
 */

Matrix NutMatrix(double Mjd_TT){
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity (Mjd_TT);

    // Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles (Mjd_TT, dpsi, deps);

    // Transformation from mean to true equator and equinox
    Matrix NutMat = R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);
    return NutMat;
}