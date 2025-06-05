#ifndef _MATRIX_
#define _MATRIX_

#include <vector>

class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);
        Matrix(const Matrix& m);
        ~Matrix();
        double get(int i, int j);
        void set(int i, int j, double d);
        int getNumCol();
        int getNumFilas();
        std::vector<double> getCol(int col); //Obtenemos el vector columna de una matriz
        std::vector<double> getFila(int col); //Obtenemos el vector fila de una matriz
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator-(const Matrix& matrix2);
        Matrix  operator*(const Matrix& matrix2);
        double& operator()(const int i, const int j) const;
        Matrix  operator+(double c);//Redefinimos para sumar un escalar a toda las entradas
        Matrix  operator*(double c);//Redefinimos para multiplicar por un escalar
        std::vector<double> operator*(std::vector<double>& vec); //Redefinimos para multiplicar por un vector
        Matrix  operator/(double c);//Redefinimos para dividir entre un escalar
        double norm();//Norma de Frobenius. Sustituye el norm() de Matlab.
        static double normaVector(std::vector<double> vec); //Calcula la norma de un vector
        double normaVector(); //Calcula la norma de una matriz como si fuera un vector
        static double productoEscalar(std::vector<double>, std::vector<double>); //Calcula el producto escalar de 2 vectores
        double productoEscalar(Matrix & m2); //Calcula el producto escalar entre dos matrices fila
        static std::vector<double> productoVectorial(std::vector<double> vec1, std::vector<double> vec2); //Calcula el producto vectorial de dos vectores de 3 entradas
        Matrix productoVectorial(Matrix& m2); //Realiza el producto vectorial entre dos matrices fila
        static std::vector<double> unit(const std::vector<double>& vec); //Recibe un vector como parámetro y devuelve el vector unitario correspondiente
        static double angl(const std::vector<double> vec1, const std::vector<double> vec2); //Dados dos vectores, devolvemos el angulo entre ellos
        Matrix transponer(); //Devuelve la matriz traspuesta
        Matrix inversa(); //Devuelve la inversa de la matriz
        static Matrix identidad(int n); //Construye una matriz identidad de tamaño n x n
        double determinante();
        void toString();

        void print();
        bool equalMatrix(const Matrix& m1, const Matrix& m2, double tolerancia);
 
    private:
        void initMatrix();
 
    private:
        int fil;
        int col;
        double **matrix;

};

#endif
