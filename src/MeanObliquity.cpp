#include "MeanObliquity.h"
#include "Const.h"
#include <cmath>

/**
 * @file MeanObliquity.cpp
 * @brief Computes the mean obliquity of the ecliptic
 *
 * @details
 * This function computes the mean obliquity of the ecliptic.
 *
 * @param[in] Mjd_TT    Modified Julian Date (Terrestrial Time)
 * @return double       Mean obliquity of the ecliptic [rad]
 *
 * @note
 * Last modified: 2015/08/12   M. Mahooti
 */

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - Const::MJD_J2000) / 36525.0;

    double MOblq = Const::Rad * (84381.448 / 3600.0 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);

    return MOblq;
}