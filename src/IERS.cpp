#include "IERS.h"
#include "Const.h"
#include <cmath>

/**
 * @file IERS.cpp
 * @brief IERS: Management of IERS time and polar motion data
 *
 * This function retrieves IERS Earth rotation parameters based on the input
 * Modified Julian Date (Mjd_UTC) and the provided Earth Orientation Parameters (eop).
 *
 * @param eop         Matrix containing Earth Orientation Parameters (eop) data
 * @param Mjd_UTC     Modified Julian Date (UTC) for which parameters are requested
 * @param interp      Interpolation type ('n' for none, 'l' for linear)
 * @param x_pole      Output: Pole coordinate in X direction [rad]
 * @param y_pole      Output: Pole coordinate in Y direction [rad]
 * @param UT1_UTC     Output: UT1-UTC time difference [s]
 * @param LOD         Output: Length of day [s]
 * @param dpsi        Output: Nutation correction in longitude [rad]
 * @param deps        Output: Nutation correction in obliquity [rad]
 * @param dx_pole     Output: Rate of change of X pole coordinate [rad/day]
 * @param dy_pole     Output: Rate of change of Y pole coordinate [rad/day]
 * @param TAI_UTC     Output: TAI-UTC time difference [s]
 *
 * @note The function assumes the input eop matrix is structured such that each row
 *       contains the following parameters in sequence: MJD, x_pole, y_pole, UT1-UTC,
 *       LOD, dpsi, deps, dx_pole, dy_pole, TAI-UTC.
 *
 * @note If interp is 'l', linear interpolation is used between nearest dates in eop.
 *       If interp is 'n' or any other character, the nearest date's data from eop is used.
 *
 * @see const.h for unit conversions and constants used in the calculations.
 *
 * @return None. Outputs are passed by reference.
 *
 * @author M. Mahooti
 * @date 2018/02/01
 */

void IERS(Matrix& eop, double Mjd_UTC, char interp,
          double& x_pole, double& y_pole, double& UT1_UTC, double& LOD,
          double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {


    int i;
    double mfme, fixf;

    // Definir constantes si es necesario
    const double pi2 = Const::doublepi;

    // Si no se proporciona el tipo de interpolación, se asume 'n'
    if (interp != 'l' && interp != 'n') {
        interp = 'n';
    }

    if (interp == 'l') {
        // Interpolación lineal
        int mjd = static_cast<int>(floor(Mjd_UTC));
        for (i = 0; i < eop.getNumFilas(); i++) {
            if (mjd == static_cast<int>(eop.get(i, 3))) {
                break;
            }
        }

        Matrix preeop(1, 13);
        Matrix nexteop(1, 13);

        // Extraer datos de interpolación
        for (int j = 0; j < 13; j++) {
            preeop.set(0, j, eop.get(i, j));
            nexteop.set(0, j, eop.get(i + 1, j));
        }

        mfme = 1440.0 * (Mjd_UTC - floor(Mjd_UTC));
        fixf = mfme / 1440.0;

        // Interpolación de los parámetros de rotación de la Tierra
        x_pole = preeop.get(0, 4) + (nexteop.get(0, 4) - preeop.get(0, 4)) * fixf;
        y_pole = preeop.get(0, 5) + (nexteop.get(0, 5) - preeop.get(0, 5)) * fixf;
        UT1_UTC = preeop.get(0, 6) + (nexteop.get(0, 6) - preeop.get(0, 6)) * fixf;
        LOD = preeop.get(0, 7) + (nexteop.get(0, 7) - preeop.get(0, 7)) * fixf;
        dpsi = preeop.get(0, 8) + (nexteop.get(0, 8) - preeop.get(0, 8)) * fixf;
        deps = preeop.get(0, 9) + (nexteop.get(0, 9) - preeop.get(0, 9)) * fixf;
        dx_pole = preeop.get(0, 10) + (nexteop.get(0, 10) - preeop.get(0, 10)) * fixf;
        dy_pole = preeop.get(0, 11) + (nexteop.get(0, 11) - preeop.get(0, 11)) * fixf;
        TAI_UTC = preeop.get(0, 12);

        // Convertir a unidades adecuadas
        x_pole /= Const::Arcs;  // Coordenada del polo [rad]
        y_pole /= Const::Arcs;  // Coordenada del polo [rad]
        dpsi /= Const::Arcs;
        deps /= Const::Arcs;
        dx_pole /= Const::Arcs; // Coordenada del polo [rad]
        dy_pole /= Const::Arcs; // Coordenada del polo [rad]

    } else if (interp == 'n') {
        // Sin interpolación
        int mjd = static_cast<int>(floor(Mjd_UTC));
        for (i = 0; i < eop.getNumFilas(); i++) {
            if (mjd == static_cast<int>(eop.get(i, 3))) {
                break;
            }
        }

        // Extraer datos sin interpolación
        x_pole = eop.get(i, 4) / Const::Arcs;  // Coordenada del polo [rad]
        y_pole = eop.get(i, 5) / Const::Arcs;  // Coordenada del polo [rad]
        UT1_UTC = eop.get(i, 6);              // Diferencia de tiempo UT1-UTC [s]
        LOD = eop.get(i, 7);                  // Longitud del día [s]
        dpsi = eop.get(i, 8) / Const::Arcs;
        deps = eop.get(i, 9) / Const::Arcs;
        dx_pole = eop.get(i, 10) / Const::Arcs; // Coordenada del polo [rad]
        dy_pole = eop.get(i, 11) / Const::Arcs; // Coordenada del polo [rad]
        TAI_UTC = eop.get(i, 12);             // Diferencia de tiempo TAI-UTC [s]
    }
}