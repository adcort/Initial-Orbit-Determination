#include <cmath>
#include "gast.h"
#include "gmst.h"
#include "EqnEquinox.h"
#include "Const.h"

/**
 * @file gast.cpp
 * @brief Greenwich Apparent Sidereal Time (GAST)
 *
 * @details
 * Purpose:
 *   Calculates the Greenwich Apparent Sidereal Time (GAST) for a given
 *   Modified Julian Date (Mjd_UT1).
 *
 * Input:
 *   - Mjd_UT1: Modified Julian Date UT1
 *
 * Output:
 *   - gstime: GAST in radians
 *
 * Last modified: 2015/08/12 by M. Mahooti
 *
 */

double gast(double Mjd_UT1){
    return fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*Const::pi );
}