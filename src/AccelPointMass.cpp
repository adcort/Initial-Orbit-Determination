#include "../include/AccelPointMass.h"
/**
 * @file AccelPointMass.cpp
 * @brief Computes the perturbational acceleration due to a point mass.
 *
 * @param r Satellite position vector.
 * @param s Point mass position vector.
 * @param GM Gravitational coefficient of point mass.
 * @return Acceleration (a = d^2r/dt^2).
 *
 * @note Last modified: 2018/01/27 by M. Mahooti.
 */
void AccelPointMass(Matrix& r, Matrix& s, double GM, Matrix& a){


    // Relative position vector of satellite w.r.t. point mass
    Matrix d(1,3);
    d = r-s;

    // Acceleration
    a = (d/(d.normaVector()*d.normaVector()*d.normaVector()) + s/(s.normaVector()*s.normaVector()*s.normaVector())) * (-GM); //Multiplicamos por el escalar al final para que no dé problemas.
}
