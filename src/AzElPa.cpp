#include "AzElPa.h"
#include "Const.h"
#include <cmath>
/**
 * @file AzElPa.cpp
 * @brief Computes azimuth, elevation and partials from local tangent coordinates.
 *
 * Computes azimuth (A), elevation (E) and their partial derivatives (partials)
 * with respect to local tangent coordinates (East-North-Zenith frame).
 *
 * @param s Topocentric local tangent coordinates [East, North, Zenith].
 * @return std::tuple<double, double, std::array<double, 3>, std::array<double, 3>>
 *         A tuple containing:
 *         - A: Azimuth [rad]
 *         - E: Elevation [rad]
 *         - dAds: Partials of azimuth w.r.t. s
 *         - dEds: Partials of elevation w.r.t. s
 *
 * @note Assumes s is a 3-element vector with components in the East-North-Zenith frame.
 * @note Azimuth A ranges from 0 to 2*pi.
 * @note Last modified: 2015/08/12 by M. Mahooti
 */

void AzElPa(Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds){

    //pi2 está declarada en global.h
    double rho = sqrt(s.get(0,0) * s.get(0,0) + s.get(0,1) * s.get(0,1));

    //Angles
    Az = atan2(s.get(0,0), s.get(0,1));

    if(Az<0.0){
        Az = Az+Const::doublepi;
    }

    El = atan(s.get(0,2)/rho);

    //Calculamos dot ya que esta función no existe en C++
    double dot = s.get(0,0) * s.get(0,0) + s.get(0,1) * s.get(0,1) + s.get(0,2) * s.get(0,2);

    //Partials
    dAds.set(0,0, s.get(0,1) / (rho * rho));
    dAds.set(0,1, -s.get(0,0) / (rho * rho));
    dAds.set(0,2, 0.0);
    dEds.set(0,0, -s.get(0,0)*s.get(0,2) / rho / dot);
    dEds.set(0,1, -s.get(0,1)*s.get(0,2) / rho / dot);
    dEds.set(0,2, rho / dot );
}
