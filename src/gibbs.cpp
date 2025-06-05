
#include "gibbs.h"
#include "Matrix.h"
#include "Const.h"
#include <cmath>
/**
 * @file gibbs.cpp
 * \brief Performs Gibbs method of orbit determination.
 *
 * This function calculates the velocity at the midpoint of 3 given
 * position vectors using the Gibbs method.
 *
 * \param r1 ijk position vector #1 [m]
 * \param r2 ijk position vector #2 [m]
 * \param r3 ijk position vector #3 [m]
 * \param v2 Output ijk velocity vector for r2 [m/s]
 * \param theta Output angle between vectors [rad]
 * \param error Output flag indicating success ('ok', ...)
 */


void gibbs(std::vector<double>& r1,  std::vector<double>& r2,  std::vector<double>& r3,
           std::vector<double>& v2, double& theta, double& theta1, double& copa, std::string& error){

    double small= 0.00001;
    theta= 0.0;
    error = "          ok";
    theta1= 0.0;

    double magr1 = Matrix::normaVector(r1);
    double magr2 = Matrix::normaVector(r2);
    double magr3 = Matrix::normaVector(r3);

    std::vector<double> p = Matrix::productoVectorial(r2, r3);
    std::vector<double> q = Matrix::productoVectorial(r3, r1);
    std::vector<double> w = Matrix::productoVectorial(r1, r2);
    std::vector<double> pn = Matrix::unit(p);
    std::vector<double> r1n = Matrix::unit(r1);
    copa = asin(Matrix::productoEscalar(pn, r1n));

    if (std::abs(Matrix::productoEscalar(r1n, pn)) > 0.017452406) {
        error = "not coplanar";
    }

    std::vector<double> d = { p[0] + q[0] + w[0], p[1] + q[1] + w[1], p[2] + q[2] + w[2] };
    double magd = Matrix::normaVector(d);
    std::vector<double> n = { magr1 * p[0] + magr2 * q[0] + magr3 * w[0],
                              magr1 * p[1] + magr2 * q[1] + magr3 * w[1],
                              magr1 * p[2] + magr2 * q[2] + magr3 * w[2] };
    double magn = Matrix::normaVector(n);
    std::vector<double> nn = Matrix::unit(n);
    std::vector<double> dn = Matrix::unit(d);

    // -------------------------------------------------------------
    //determine if  the orbit is possible. both d and n must be in
    //the same direction, and non-zero.
    // -------------------------------------------------------------
    if (std::abs(magd) < small || std::abs(magn) < small || Matrix::productoEscalar(nn, dn) < small) {
        error = "  impossible";
    } else {
        theta = Matrix::angl(r1, r2);
        theta1 = Matrix::angl(r2, r3);

        // ----------- perform gibbs method to find v2 -----------
        double r1mr2 = magr1 - magr2;
        double r3mr1 = magr3 - magr1;
        double r2mr3 = magr2 - magr3;
        std::vector<double> s = { r1mr2 * r3[0] + r3mr1 * r2[0] + r2mr3 * r1[0],
                                  r1mr2 * r3[1] + r3mr1 * r2[1] + r2mr3 * r1[1],
                                  r1mr2 * r3[2] + r3mr1 * r2[2] + r2mr3 * r1[2] };
        std::vector<double> b = Matrix::productoVectorial(d, r2);
        double l = sqrt(Const::GM_Earth / (magd * magn));
        double tover2 = l / magr2;
        for (int i = 0; i < 3; ++i) {
            v2[i] = tover2 * b[i] + l * s[i];
        }
    }
}