#include "Position.h"

/**
 * @file Position.cpp
 * @brief Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 * @author M. Mahooti
 * @date 2015/08/12
 */


void Position(double lon, double lat, double h, double result[3]){
    double R_equ = Const::R_Earth;
    double f = Const::f_Earth;

    double e2 = f * (2.0 - f);   // Square of eccentricity
    double CosLat = cos(lat);    // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    result[0] = (N + h) * CosLat * cos(lon);
    result[1] = (N + h) * CosLat * sin(lon);
    result[2] = ((1.0 - e2) * N + h) * SinLat;
}