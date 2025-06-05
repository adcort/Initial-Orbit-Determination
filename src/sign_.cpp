#include "../include/sign_.h"
#include <cmath>
/**
 * @file sign_.cpp
 * @brief Returns the absolute value of a with the sign of b
 *
 * @details
 * This function returns the absolute value of parameter a but with the sign of parameter b.
 *
 * @param a Input value
 * @param b Sign reference
 * @return double Resulting value with the sign of b
 */
double sign_(double a, double b){
    if(b>=0.0) return abs(a);
    else return -abs(a);
}