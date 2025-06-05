#include <cstdarg>
#include <fstream>
#include <sstream>
#include "../include/global.h"


/**
 * @file global.cpp
 * @brief global functions
 *
 * @details
 * Purpose:
 * Contains the global functions that will be used by any other function of the project
 *
 */

Matrix *global::eopdata;

void global::eop19620101(int c) {
    global::eopdata = new Matrix(13, c);
    FILE *fid = fopen("../data/eop19620101.txt", "r");

    if(fid==NULL){
        printf("Error");
        exit(EXIT_FAILURE);
    }
    for(int i=1; i<=c; i++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", (&(*global::eopdata)(1, i))
                , (&(*global::eopdata)(2, i)), (&(*global::eopdata)(3, i)), (&(*global::eopdata)(4, i)), (&(*global::eopdata)(5, i))
                , (&(*global::eopdata)(6, i)), (&(*global::eopdata)(7, i)), (&(*global::eopdata)(8, i)), (&(*global::eopdata)(9, i))
                , (&(*global::eopdata)(10, i)), (&(*global::eopdata)(11, i)), (&(*global::eopdata)(12, i)), (&(*global::eopdata)(13, i)));
    }
    fclose(fid);
}


int global::nargin(int arg1, ...) {
    va_list args;
    va_start(args, arg1);

    int count = 1; // Inicializar el contador con el primer argumento
    while (va_arg(args, int)) {
        count++;
    }

    va_end(args);

    return count;
}