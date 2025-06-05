#include "LTC.h"
#include "R_y.h"
#include "R_z.h"
/**
 * @file LTC.cpp
 * @brief Transformation from Greenwich meridian system to local tangent coordinates
 *
 * @details
 * This function performs the transformation from the Greenwich meridian system to local tangent coordinates.
 *
 * @param[in] lon    Geodetic East longitude [rad]
 * @param[in] lat    Geodetic latitude [rad]
 * @param[out] M     Rotation matrix from the Earth equator and Greenwich meridian
 *                   to the local tangent (East-North-Zenith) coordinate system
 *
 * @note
 * Last modified: 2015/08/12   M. Mahooti
 */
Matrix LTC(double lon, double lat){

    Matrix m = R_y(-lat) * R_z(lon);
    for (int j = 0; j < 3; ++j) {
        double Aux = m.get(0,j);
        m.set(0,j, m.get(1,j));
        m.set(1,j, m.get(2,j));
        m.set(2,j, Aux);
    }

    return m;
}