#include "EqnEquinox.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include <cmath>

/**
 * @file EqnEquinox.cpp
 * @brief Computation of the equation of the equinoxes
 *
 * @details
 * The equation of the equinoxes \f$d\psi\cos(\epsilon)\f$ is the right ascension of
 * the mean equinox referred to the true equator and equinox and is equal
 * to the difference between apparent and mean sidereal time.
 *
 * @param[in] Mjd_TT    Modified Julian Date (Terrestrial Time)
 * @return double       Equation of the equinoxes
 *
 * @note
 * Last modified:   2015/08/12   M. Mahooti
 */

double EqnEquinox(double Mjd_TT){

    double dpsi, deps;
    // Nutation in longitude and obliquity
    NutAngles(Mjd_TT, dpsi, deps);
    // Equation of the equinoxes
    return dpsi * cos(MeanObliquity(Mjd_TT));
}