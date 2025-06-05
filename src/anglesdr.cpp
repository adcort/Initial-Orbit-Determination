#include "anglesdr.h"
#include "Const.h"
#include "Geodetic.h"
#include "Matrix.h"
#include "IERS.h"
#include "LTC.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "doubler.h"
#include <string>
#include <cmath>


extern Matrix eopdata;
void anglesdr(double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3,
              std::vector<double>& rsite1, std::vector<double>& rsite2, std::vector<double>& rsite3,
              std::vector<double>& r2, std::vector<double>& v2){

    double magr1in = 1.1*Const::R_Earth;
    double magr2in = 1.11*Const::R_Earth;
    char direct = 'y';

    double tol    = 1e-8*Const::R_Earth;
    double pctchg = 0.005;

    double t1 = (Mjd1 - Mjd2)*86400.0;
    double t3 = (Mjd3 - Mjd2)*86400.0;

    std::vector<double> los1 = {cos(el1) * sin(az1), cos(el1) * cos(az1), sin(el1)};
    std::vector<double> los2 = {cos(el2) * sin(az2), cos(el2) * cos(az2), sin(el2)};
    std::vector<double> los3 = {cos(el3) * sin(az3), cos(el3) * cos(az3), sin(el3)};

    double lon1, lat1, h1, lon2, lat2, h2, lon3, lat3, h3;

    Geodetic(lon1, lat1, h1, rsite1);
    Geodetic(lon2, lat2, h2, rsite2);
    Geodetic(lon3, lat3, h3, rsite3);

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);

    // body-fixed system
    los1 = {M1(0, 0) * los1[0] + M1(0, 1) * los1[1] + M1(0, 2) * los1[2],
            M1(1, 0) * los1[0] + M1(1, 1) * los1[1] + M1(1, 2) * los1[2],
            M1(2, 0) * los1[0] + M1(2, 1) * los1[1] + M1(2, 2) * los1[2]};
    los2 = {M2(0, 0) * los2[0] + M2(0, 1) * los2[1] + M2(0, 2) * los2[2],
            M2(1, 0) * los2[0] + M2(1, 1) * los2[1] + M2(1, 2) * los2[2],
            M2(2, 0) * los2[0] + M2(2, 1) * los2[1] + M2(2, 2) * los2[2]};
    los3 = {M3(0, 0) * los3[0] + M3(0, 1) * los3[1] + M3(0, 2) * los3[2],
            M3(1, 0) * los3[0] + M3(1, 1) * los3[1] + M3(1, 2) * los3[2],
            M3(2, 0) * los3[0] + M3(2, 1) * los3[1] + M3(2, 2) * los3[2]};

    // mean of date system (J2000)
    /*double Mjd_UTC = Mjd1;

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los1 = {E(0, 0) * los1[0] + E(0, 1) * los1[1] + E(0, 2) * los1[2],
            E(1, 0) * los1[0] + E(1, 1) * los1[1] + E(1, 2) * los1[2],
            E(2, 0) * los1[0] + E(2, 1) * los1[1] + E(2, 2) * los1[2]};
    rsite1 = {E(0, 0) * rsite1[0] + E(0, 1) * rsite1[1] + E(0, 2) * rsite1[2],
              E(1, 0) * rsite1[0] + E(1, 1) * rsite1[1] + E(1, 2) * rsite1[2],
              E(2, 0) * rsite1[0] + E(2, 1) * rsite1[1] + E(2, 2) * rsite1[2]};

    Mjd_UTC = Mjd2;
    IERS(eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los2 = {E(0, 0) * los2[0] + E(0, 1) * los2[1] + E(0, 2) * los2[2],
            E(1, 0) * los2[0] + E(1, 1) * los2[1] + E(1, 2) * los2[2],
            E(2, 0) * los2[0] + E(2, 1) * los2[1] + E(2, 2) * los2[2]};
    rsite2 = {E(0, 0) * rsite2[0] + E(0, 1) * rsite2[1] + E(0, 2) * rsite2[2],
              E(1, 0) * rsite2[0] + E(1, 1) * rsite2[1] + E(1, 2) * rsite2[2],
              E(2, 0) * rsite2[0] + E(2, 1) * rsite2[1] + E(2, 2) * rsite2[2]};

    Mjd_UTC = Mjd3;
    IERS(eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los3 = {E(0, 0) * los3[0] + E(0, 1) * los3[1] + E(0, 2) * los3[2],
            E(1, 0) * los3[0] + E(1, 1) * los3[1] + E(1, 2) * los3[2],
            E(2, 0) * los3[0] + E(2, 1) * los3[1] + E(2, 2) * los3[2]};
    rsite3 = {E(0, 0) * rsite3[0] + E(0, 1) * rsite3[1] + E(0, 2) * rsite3[2],
              E(1, 0) * rsite3[0] + E(1, 1) * rsite3[1] + E(1, 2) * rsite3[2],
              E(2, 0) * rsite3[0] + E(2, 1) * rsite3[1] + E(2, 2) * rsite3[2]};


    double magr1old  = 99999999.9;
    double magr2old  = 99999999.9;
    double magrsite1 = Matrix::normaVector(rsite1);
    double magrsite2 = Matrix::normaVector(rsite2);
    double magrsite3 = Matrix::normaVector(rsite3);

    double cc1 = 2.0*Matrix::productoEscalar(los1,rsite1);
    double cc2 = 2.0*Matrix::productoEscalar(los2,rsite2);
    double ktr = 0;

    double r3, f1, f2, q1, magr1, magr2, a, deltae32;
    while (abs(magr1in-magr1old) > tol || abs(magr2in-magr2old) > tol)
        ktr = ktr + 1;
    doubler(  cc1,  cc2,  magrsite1,  magrsite2,  magr1in,  magr2in, los1, los2,  los3,rsite1,  rsite2,  rsite3,t1,  t3,  direct,r2,  r3,  f1,  f2,  q1,magr1,  magr2, a,  deltae32);

    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = t3 - sqrt(a^3/const.GM_Earth)*(deltae32-sin(deltae32));
    v2 = (r3 - f*r2)/g;

    magr1o = magr1in;
    magr1in = (1.0+pctchg)*magr1in;
    deltar1 = pctchg*magr1in;
    [r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
    los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
    pf1pr1 = (f1delr1-f1)/deltar1;
    pf2pr1 = (f2delr1-f2)/deltar1;

    magr1in = magr1o;
    deltar1 = pctchg*magr1in;
    magr2o = magr2in;
    magr2in = (1.0+pctchg)*magr2in;
    deltar2 = pctchg*magr2in;
    [r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
    los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
    pf1pr2 = (f1delr2-f1)/deltar2;
    pf2pr2 = (f2delr2-f2)/deltar2;

    magr2in = magr2o;
    deltar2 = pctchg*magr2in;

    delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
    delta1 = pf2pr2*f1 - pf1pr2*f2;
    delta2 = pf1pr1*f2 - pf2pr1*f1;

    deltar1 = -delta1/delta;
    deltar2 = -delta2/delta;

    magr1old = magr1in;
    magr2old = magr2in;

    magr1in = magr1in + deltar1;
    magr2in = magr2in + deltar2;

    end

    [r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
    los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);

    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = t3 - sqrt(a^3/const.GM_Earth)*(deltae32-sin(deltae32));
    v2 = (r3 - f*r2)/g;*/

}