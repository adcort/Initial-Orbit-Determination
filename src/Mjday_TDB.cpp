#include "Mjday_TDB.h"
#include <cmath>
/**
 * @file Mjday_TDB.cpp
 * @brief Computes the Modified Julian Date for barycentric dynamical time
 *
 * @details
 * This function computes the Modified Julian Date for barycentric dynamical time (TDB).
 *
 * @param[in] Mjd_TT    Modified Julian Date (Terrestrial Time)
 * @return double       Modified Julian Date (TDB)
 *
 * @note
 * Reference: Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill; New York; 3rd edition (2007).
 * Last modified: 2015/08/12   M. Mahooti
 */

double Mjday_TDB(double Mjd_TT) {
    // Compute Julian Centuries of TT
    double T_TT = (Mjd_TT - 51544.5) / 36525;

    // Compute Modified Julian Date of TDB
    double Mjd_TDB = Mjd_TT + (0.001658 * sin(628.3076 * T_TT + 6.2401) +
                               0.000022 * sin(575.3385 * T_TT + 4.2970) +
                               0.000014 * sin(1256.6152 * T_TT + 6.1969) +
                               0.000005 * sin(606.9777 * T_TT + 4.0212) +
                               0.000005 * sin(52.9691 * T_TT + 0.4444) +
                               0.000002 * sin(21.3299 * T_TT + 5.5431) +
                               0.000010 * sin(628.3076 * T_TT + 4.2490)) / 86400;

    return Mjd_TDB;
}