#ifndef PROYECTO_HGIBBS_H
#define PROYECTO_HGIBBS_H
#include <vector>
#include <string>

void hgibbs(const std::vector<double>& r1, const std::vector<double>& r2, const std::vector<double>& r3,
            double Mjd1, double Mjd2, double Mjd3,
            std::vector<double>& v2, double& theta, double& theta1, double& copa, std::string& error);

#endif //PROYECTO_HGIBBS_H
