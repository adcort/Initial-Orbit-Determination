#include "unit.h"
#include <cmath>
#include "Matrix.h"

/**
 * @file unit.cpp
 * @brief Calculates a unit vector given the original vector
 *
 * @details
 * This function calculates a unit vector given the original vector. If a
 * zero vector is input, the vector is set to zero.
 *
 * @param vec Vector
 * @return std::vector<double> Unit vector
 */


std::vector<double> unit(std::vector<double> vec) {
    double small = 0.000001;
    double magv = Matrix::normaVector(vec);
    std::vector<double> outvec(3, 0.0);
    if (magv > small) {
        for (int i = 0; i < 3; ++i) {
            outvec[i] = vec[i] / magv;
        }
    } else {
        for(int i = 0; i < 3; i++){
            outvec[i] = 0.0;
        }
    }
    return outvec;
}