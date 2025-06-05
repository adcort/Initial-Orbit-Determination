#ifndef PROYECTO_MEASUPDATE_H
#define PROYECTO_MEASUPDATE_H

#include "Matrix.h"
#include <vector>
void MeasUpdate(Matrix& x, Matrix& z, Matrix& g,  std::vector<double>& s,  Matrix& G, Matrix& P, int n, Matrix& K);

#endif //PROYECTO_MEASUPDATE_H
