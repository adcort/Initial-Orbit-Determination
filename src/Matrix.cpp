#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @file Matrix.cpp
 * @brief Representates a Matrix. Contains all the necessary functions to operate
 * with them and with any given vector.
 *
 */

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}

Matrix::Matrix(int fil, int col, double v[], int n) : fil(fil), col(col) {
    initMatrix();

    int k = 0;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n) { // Verificar si aún quedan elementos en el vector v
                matrix[i][j] = v[k];
                k++;
            } else {
                // Si se agotan los elementos del vector, salir del bucle
                break;
            }
        }
}

 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}

double Matrix::get(int i, int j) {
    return this->matrix[i][j];
}

void Matrix::set(int i, int j, double d) {
    this->matrix[i][j] = d;
}

int Matrix::getNumCol() {
    return this->col;
}

int Matrix::getNumFilas() {
    return this->fil;
}
 
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];
 
    return *this;
}
 
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}

Matrix Matrix::operator*(const Matrix& matrix2)
{
    // Verificar las dimensiones antes de la multiplicación
    if (this->col != matrix2.fil) {
        throw std::invalid_argument("No se pueden multiplicar las matrices: dimensiones incompatibles.");
    }

    Matrix result(this->fil, matrix2.col); // Matriz resultante de tamaño adecuado

    for (int i = 0; i < this->fil ; ++i) {
        for (int j = 0; j < matrix2.col; ++j) {
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; ++k) {
                result.matrix[i][j] += this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }

    return result;
}
//Redefinimos para multiplicar por un escalar
Matrix  Matrix::operator*(double c)
{
    Matrix result(fil, col);

    for (int i = 0; i < this->fil; i++) {
        for (int j = 0; j < this->col; j++) {
            // Sumar el escalar c a cada entrada de la matriz original
            result.matrix[i][j] = this->matrix[i][j] * c;
        }
    }
    return result;
}

//Redefinimos para sumar un escalar a todas las entradas
Matrix Matrix::operator+(double c) {
    Matrix result(fil, col);

    for (int i = 0; i < this->fil; i++) {
        for (int j = 0; j < this->col; j++) {
            // Sumar el escalar c a cada entrada de la matriz original
            result.matrix[i][j] = this->matrix[i][j] + c;
        }
    }
    return result;
}


Matrix  Matrix::operator/(double c)
{
    Matrix result(fil, col);

    for (int i = 0; i < this->fil; i++) {
        for (int j = 0; j < this->col; j++) {
            // Sumar el escalar c a cada entrada de la matriz original
            result.matrix[i][j] = this->matrix[i][j] / c;
        }
    }
    return result;
}
 
 
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}
 
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

bool Matrix::equalMatrix(const Matrix &m1, const Matrix &m2, double tolerancia) {

    if (m1.fil != m2.fil || m1.col != m2.col) {
        return false;
    }
    //A partir de aquí tienen mismo número de filas y columnas
    for (int i = 0; i < m1.fil; i++) {
        for (int j = 0; j < m1.col; j++){
            if (abs(m1.matrix[i][j] - m2.matrix[i][j]) > tolerancia) {
                return false;
            }
        }
    }
    return true;
}
//Calcula la norma de Frobenius de una matriz.
double Matrix::norm(){
    double suma = 0.0;
    for(int i=0; i<this->fil; i++){
        for(int j=0; j<this->col; j++){
            suma += this->matrix[i][j] * this->matrix[i][j];
        }
    }
    return suma;
}
//Calcula la norma de la matriz como si fuera un vector
double Matrix::normaVector(){
    double suma = 0.0;
    for(int i=0; i<this->fil; i++){
        for(int j=0; j<this->col; j++){
            suma += this->matrix[i][j] * this->matrix[i][j];
        }
    }
    return std::sqrt(suma);
}


//Calcula la norma de un vector
double Matrix::normaVector(std::vector<double> vec){
    double sumaCuadrados = 0.0;
    for (double x : vec) {
        sumaCuadrados += x * x;
    }
    return sqrt(sumaCuadrados);
}

//Calcula el producto escalar de dos vectores
double Matrix::productoEscalar(std::vector<double> vec1, std::vector<double> vec2){
    double result = 0.0;
    if(vec1.size() != vec2.size()){return -1;}//Si tienen diferente tamaño devolvemos -1

    for (int i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

double Matrix::productoEscalar(Matrix & m2){
    double result = 0.0;
    if(fil!=m2.getNumFilas() || col!=m2.getNumCol()){
        throw std::invalid_argument("Las matrices no son del tamaño adecuado para el producto escalar.");
    }

    for (int i = 0; i < col; ++i) {
        result += matrix[0][i] * m2.get(0,i);
    }
    return result;
}

//Calcula el producto vectorial de dos vectores de 3 entradas
std::vector<double> Matrix::productoVectorial(std::vector<double> vec1, std::vector<double> vec2){
    if (vec1.size() != 3 || vec2.size() != 3) {
        throw std::invalid_argument("Ambos vectores deben tener exactamente 3 elementos.");
    }

    std::vector<double> resultado(3);
    resultado[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    resultado[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    resultado[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

    return resultado;
}

Matrix Matrix::productoVectorial(Matrix& m2){
    if (col!=3 || m2.getNumCol() != 3) {
        throw std::invalid_argument("Ambos vectores deben tener exactamente 3 elementos.");
    }

    Matrix resultado(1,3);
    resultado.set(0,0, matrix[0][1]*m2.get(0,2) - matrix[0][2] *m2.get(0,1));
    resultado.set(0,1, matrix[0][2]*m2.get(0,0) - matrix[0][0] *m2.get(0,2));
    resultado.set(0,2, matrix[0][0]*m2.get(0,1) - matrix[0][1] *m2.get(0,0));

    return resultado;
}

//Recibe un vector como parámetro y devuelve el vector unitario correspondiente
std::vector<double> Matrix::unit(const std::vector<double>& vec) {
    double vecNorm = 0.0;
    for (double val : vec) {
        vecNorm += val * val;
    }
    vecNorm = std::sqrt(vecNorm);

    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] / vecNorm;
    }
    return result;
}

//Dados dos vectores, devolvemos el angulo entre ellos
double Matrix::angl(const std::vector<double> vec1, const std::vector<double> vec2){
    double dot_product = 0.0;
    double mag_vec1 = 0.0;
    double mag_vec2 = 0.0;

    // Calculamos el producto punto entre los dos vectores
    for (size_t i = 0; i < vec1.size(); ++i) {
        dot_product += vec1[i] * vec2[i];
        mag_vec1 += vec1[i] * vec1[i];
        mag_vec2 += vec2[i] * vec2[i];
    }

    // Calculamos las magnitudes de los vectores
    mag_vec1 = std::sqrt(mag_vec1);
    mag_vec2 = std::sqrt(mag_vec2);

    // Calculamos el ángulo en radianes
    double angle_rad = std::acos(dot_product / (mag_vec1 * mag_vec2));

    return angle_rad;
}

//Obtenemos el vector columna de una matriz
std::vector<double> Matrix::getCol(int col) {
    if (col < 0 || col >= this->col) {
        // Manejo de error: índice de columna fuera de rango
        throw std::out_of_range("Índice de columna fuera de rango");
    }

    std::vector<double> column;
    for (int i = 0; i < this->fil; ++i) {
        column.push_back(matrix[i][col]);
    }
    return column;
}

//Obtenemos el vector fila de una matriz
std::vector<double> Matrix::getFila(int fila) {
    if (fila < 0 || fila >= this->fil) {
        // Manejo de error: índice de fila fuera de rango
        throw std::out_of_range("Índice de fila fuera de rango");
    }
    return std::vector<double>(matrix[fila], matrix[fila] + this->col);
}

//Calcula la traspuesta de la matriz
Matrix Matrix::transponer() {
    Matrix traspuesta(this->col, this->fil);
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            traspuesta.matrix[j][i] = this->matrix[i][j];
        }
    }
    return traspuesta;
}
//Redefinimos la multiplicación de matriz por un vector
std::vector<double> Matrix::operator*(std::vector<double>& vec) {
    if (vec.size() != col) {//Comprobamos que se pueden multiplicar
        throw std::invalid_argument("El tamaño del vector no coincide con el número de columnas de la matriz.");
    }

    std::vector<double> result(fil, 0.0);
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}

//Construye una matriz identidad de tamaño n x n
Matrix Matrix::identidad(int n) {
    Matrix identidad(n, n); // Creamos una matriz cuadrada de tamaño n x n

    // Rellenamos la diagonal principal con unos
    for (int i = 0; i < n; ++i) {
        identidad.matrix[i][i] = 1.0;
    }
    return identidad;
}

void Matrix::toString(){
    for (int i = 0; i < this->fil; ++i) {
        for (int j = 0; j < this->col; ++j) {
            std::cout << this->matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}


double Matrix::determinante() {
    // Verificar si la matriz es cuadrada
    if (fil != col) {
        throw std::runtime_error("La matriz no es cuadrada, no se puede calcular el determinante.");
    }

    int n = fil; // Número de filas (y columnas)

    if (n == 1) {
        return matrix[0][0]; // Si es una matriz 1x1, el determinante es el único elemento
    }

    // Algoritmo para calcular el determinante usando eliminación gaussiana
    Matrix temp(*this); // Creamos una copia de la matriz para evitar cambios no deseados
    double det = 1.0;

    for (int i = 0; i < n; ++i) {
        // Pivoteo parcial: encontrar la fila con el mayor valor en la columna i y cambiarla con la fila actual
        int maxRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(temp.matrix[j][i]) > std::abs(temp.matrix[maxRow][i])) {
                maxRow = j;
            }
        }
        if (maxRow != i) {
            // Intercambiar filas
            std::swap(temp.matrix[i], temp.matrix[maxRow]);
            // Cambiar el signo del determinante debido al intercambio de filas
            det *= -1.0;
        }

        // Hacer ceros por debajo del pivote
        double pivot = temp.matrix[i][i];
        if (pivot == 0) {
            return 0.0; // Si el pivote es cero, el determinante es cero
        }
        for (int j = i + 1; j < n; ++j) {
            double factor = temp.matrix[j][i] / pivot;
            for (int k = i; k < n; ++k) {
                temp.matrix[j][k] -= factor * temp.matrix[i][k];
            }
        }
    }

    // Calcular el determinante multiplicando los elementos diagonales
    for (int i = 0; i < n; ++i) {
        det *= temp.matrix[i][i];
    }

    return det;
}

//Devuelve la inversa de la matriz
Matrix Matrix::inversa() {

}
