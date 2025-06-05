#include "GHAMatrix.h"
#include "gast.h"
#include "R_z.h"

/**
 * @file GHAMatrix.cpp
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 *
 * @details
 * This function performs the transformation from the true equator and equinox system to the Earth equator and Greenwich meridian system.
 *
 * @param[in] Mjd_UT1    Modified Julian Date UT1
 * @param[out] GHAmat    Greenwich Hour Angle matrix
 *
 * @note
 * Last modified: 2015/08/12   M. Mahooti
 */

Matrix GHAMatrix(double Mjd_UT1) {
    return R_z(gast(Mjd_UT1));
}