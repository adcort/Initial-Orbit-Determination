#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "Matrix.h"
#include <cstdio>
#include <cstdlib>

class global {
public:
    static Matrix *eopdata;
    static void eop19620101(int c);
    static int nargin(int arg1, ...);

};

#endif //PROYECTO_GLOBAL_H
