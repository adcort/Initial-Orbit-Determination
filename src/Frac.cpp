#include "../include/Frac.h"
#include <cmath>

/**
 * @file Frac.cpp
 * @brief Fractional part of a number (y=x-[x])
 *
 */


double Frac(double x){
    return x-floor(x);
}