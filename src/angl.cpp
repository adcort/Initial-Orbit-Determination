#include "angl.h"
#include "Matrix.h"
#include "sign_.h"
#include <cmath>
/**
 * @file angl.cpp
 * @brief Calculates the angle between two vectors.
 *
 * @param[in] vec1 The first vector. It is an Eigen::VectorXd object.
 * @param[in] vec2 The second vector. It is an Eigen::VectorXd object.
 * @return The angle between the two vectors in the range of -pi to pi.
 *
 * @code
 * #include <Eigen/Dense>
 *
 * Eigen::VectorXd vec1(3);
 * Eigen::VectorXd vec2(3);
 *
 * // Initialization of vec1 and vec2
 *
 * double theta = calculateAngle(vec1, vec2);
 */

double angl(const std::vector<double> vec1, const std::vector<double> vec2){
    double small     = 0.00000001;
    double undefined = 999999.1;

    double magv1 = Matrix::normaVector(vec1);
    double magv2 = Matrix::normaVector(vec2);

    if (magv1*magv2 > small^2){
        double temp = Matrix::productoEscalar(vec1, vec2) / (magv1 * magv2);
        if (abs( temp ) > 1.0){
            double temp= sign_(temp, temp) * 1.0;
        }
        return std::round(acos(temp) * 100.0) / 100.0; //Tenemos que devolver el resultado redondeado a las centésimas porque Matlab lo hace.
    }
    else{
        return undefined;
    }
}