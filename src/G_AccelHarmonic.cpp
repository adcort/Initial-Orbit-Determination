#include <iostream>
#include "G_AccelHarmonic.h"
#include "AccelHarmonic.h"
/**
 * @file G_AccelHarmonic.cpp
 * @brief Computes the gradient of the Earth's harmonic gravity field.
 *
 * This function computes the gradient of the Earth's harmonic gravity field
 * based on the given satellite position vector and transformation matrix.
 * It uses the gravity model degree and order specified.
 *
 * @param r Satellite position vector in the true-of-date system.
 * @param U Transformation matrix to body-fixed system.
 * @param n Gravity model degree.
 * @param m Gravity model order.
 * @return Gradient (G=da/dr) in the true-of-date system.
 *
 * @note This function assumes that the transformation matrix U and
 *       satellite position vector r are correctly formatted.
 *
 * @see AccelHarmonic
 */
Matrix G_AccelHarmonic(Matrix &r, Matrix &U, int n_max, int m_max, Matrix & GGM03S){
    double d = 1.0; //Position increment [m]

    Matrix G(3,3);
    Matrix dr(3,1);
    Matrix da100(3,1);
    Matrix da200(3,1);
    Matrix da110(3,1);
    Matrix da210(3,1);
    Matrix da120(3,1);
    Matrix da220(3,1);
    Matrix da1(3,1);
    Matrix da2(3,1);
    Matrix da3(3,1);


    for(int j = 0; j < 3; ++j){
        dr.set(j,0,d);
    }
    Matrix rAux(3,1); rAux = r-dr/2;
    AccelHarmonic(rAux, U, n_max, m_max, da100, GGM03S);
    AccelHarmonic(rAux, U, n_max, m_max, da200, GGM03S);
    AccelHarmonic(rAux, U, n_max, m_max, da110, GGM03S);
    AccelHarmonic(rAux, U, n_max, m_max, da210, GGM03S);
    AccelHarmonic(rAux, U, n_max, m_max, da120, GGM03S);
    AccelHarmonic(rAux, U, n_max, m_max, da220, GGM03S);

    da1 = da100 - da200;
    da2 = da110 - da210;
    da3 = da120 - da220;

    for(int i = 0; i < 3; ++i){
        G.set(i,0, da1.get(i,0));
        G.set(i,1, da2.get(i,1));
        G.set(i,2, da3.get(i,2));
    }


    return G;
}