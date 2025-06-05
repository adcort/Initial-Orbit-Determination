#ifndef PROYECTO_CONST_H
#define PROYECTO_CONST_H

class Const{
public:
    //% Mathematical constants
    constexpr static double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;
    constexpr static double doublepi = 2*pi;
    constexpr static double Rad = pi/180;              // Radians per degree
    constexpr static double Deg = 180/pi;              // Degrees per radian
    constexpr static double Arcs      = 3600*180/pi;         // Arcseconds per radian

    // General
    constexpr static double MJD_J2000 = 51544.5;             // Modified Julian Date of J2000
    constexpr static double T_B1950   = -0.500002108;        // Epoch B1950
    constexpr static double c_light   = 299792458.000000000; // Speed of light  [m/s]; DE430
    constexpr static double AU        = 149597870700.000000; // Astronomical unit [m]; DE430

    // Physical parameters of the Earth, Sun and Moon

    // Equatorial radius and flattening
    constexpr static double R_Earth   = 6378.1363e3;     // Earth's radius [m]; DE430
    constexpr static double f_Earth   = 1/298.257223563; //  Flattening; WGS-84
    constexpr static double R_Sun     = 696000e3;        // Sun's radius [m]; DE430
    constexpr static double R_Moon    = 1738e3;          // Moon's radius [m]; DE430

    // Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
    constexpr static double omega_Earth = 15.04106717866910/3600*Rad;   // [rad/s]; WGS-84

    // Gravitational coefficients
    constexpr static double GM_Earth    = 398600.435436e9;                  // [m^3/s^2]; DE430
    constexpr static double GM_Sun      = 132712440041.939400e9;            // [m^3/s^2]; DE430
    constexpr static double GM_Moon     = GM_Earth/81.30056907419062; // [m^3/s^2]; DE430
    constexpr static double GM_Mercury  = 22031.780000e9;                   // [m^3/s^2]; DE430
    constexpr static double GM_Venus    = 324858.592000e9;                  // [m^3/s^2]; DE430
    constexpr static double GM_Mars     = 42828.375214e9;                   // [m^3/s^2]; DE430
    constexpr static double GM_Jupiter  = 126712764.800000e9;               // [m^3/s^2]; DE430
    constexpr static double GM_Saturn   = 37940585.200000e9;                // [m^3/s^2]; DE430
    constexpr static double GM_Uranus   = 5794548.600000e9;                 // [m^3/s^2]; DE430
    constexpr static double GM_Neptune  = 6836527.100580e9;                 // [m^3/s^2]; DE430
    constexpr static double GM_Pluto    = 977.0000000000009e9;              // [m^3/s^2]; DE430

    // Solar radiation pressure at 1 AU
    constexpr static double P_Sol       = 1367/c_light; // [N/m^2] (~1367 W/m^2); IERS 96

};

#endif //PROYECTO_CONST_H
