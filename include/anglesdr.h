#ifndef PROYECTO_ANGLESDR_H
#define PROYECTO_ANGLESDR_H
#include <vector>

void anglesdr(double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3,
              const std::vector<double>& rsite1, const std::vector<double>& rsite2, const std::vector<double>& rsite3,
              std::vector<double>& r2, std::vector<double>& v2);

#endif //PROYECTO_ANGLESDR_H
