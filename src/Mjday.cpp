#include "Mjday.h"
#include "global.h"
#include <cmath>

/**
 * @file Mjday.cpp
 * @brief Converts a given date and time to the Modified Julian Date (MJD).
 *
 * @param year Year (e.g., 2024).
 * @param mon Month (1-12).
 * @param day Day (1-31).
 * @param hr Universal time hour (0-23).
 * @param min Universal time minute (0-59).
 * @param sec Universal time second (0-59).
 * @return double Modified Julian Date (MJD).
 */

double Mjday(int yr, int mon, int day, int hr, int min, double sec) {
    // Adaptamos así la función nargin de Matlab que cuenta el número de argumentos de entrada
    int numArgs = global::nargin(yr, mon, day, hr, min, sec);
    if (numArgs < 4) {
        hr = 0;
        min=0;
        sec = 0;
    }

    double jd = 367.0 * yr
                - std::floor((7 * (yr + std::floor((mon + 9) / 12.0))) * 0.25)
                + std::floor(275 * mon / 9.0)
                + day + 1721013.5
                + ((sec / 60.0 + min) / 60.0 + hr) / 24.0;

    double Mjd = jd - 2400000.5;

    return Mjd;
}
