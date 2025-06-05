cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(Proyecto
        include/R_x.h src/R_x.cpp
        include/R_y.h src/R_y.cpp
        include/R_z.h src/R_z.cpp
        include/Position.h src/Position.cpp
        include/Frac.h
        include/global.h
        include/Matrix.h
        include/sign_.h
        include/unit.h
        src/Frac.cpp
        src/global.cpp
        src/Matrix.cpp
        src/sign_.cpp
        include/timediff.h src/timediff.cpp
        include/AccelPointMass.h src/AccelPointMass.cpp
        EKF_Test.cpp include/AzElPa.h src/AzElPa.cpp include/Cheb3D.h src/Cheb3D.cpp include/EccAnom.h src/EccAnom.cpp include/Geodetic.h src/Geodetic.cpp include/MeanObliquity.h src/MeanObliquity.cpp include/Mjday.h src/Mjday.cpp include/Mjday_TDB.h src/Mjday_TDB.cpp include/NutAngles.h src/NutAngles.cpp include/Legendre.h src/Legendre.cpp src/unit.cpp include/IERS.h src/IERS.cpp include/EqnEquinox.h src/EqnEquinox.cpp include/gmst.h src/gmst.cpp include/gast.h src/gast.cpp include/GHAMatrix.h src/GHAMatrix.cpp include/TimeUpdate.h src/TimeUpdate.cpp include/angl.h src/angl.cpp include/G_AccelHarmonic.h src/G_AccelHarmonic.cpp include/LTC.h src/LTC.cpp include/PoleMatrix.h src/PoleMatrix.cpp include/NutMatrix.h src/NutMatrix.cpp include/elements.h src/elements.cpp include/PrecMatrix.h src/PrecMatrix.cpp include/gibbs.h src/gibbs.cpp include/hgibbs.h src/hgibbs.cpp include/doubler.h src/doubler.cpp include/anglesdr.h src/anglesdr.cpp include/AccelHarmonic.h src/AccelHarmonic.cpp include/MeasUpdate.h src/MeasUpdate.cpp include/VarEqn.h src/VarEqn.cpp include/VarEqn.h src/VarEqn.cpp)
