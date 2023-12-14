#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <omp.h>
#include <stdio.h>
#include <ctime>
#include <string>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

double elrand (unsigned /*i*/, unsigned /*j*/)
{
    return 10*(rand()/(double)RAND_MAX)-5;
}

double baggsin (unsigned i, unsigned j)
{
    return sin((i+1)*(j+2));
}

using frodo = double (*) (unsigned, unsigned);


/**
 * This class is what we're supposed to work with
*/

template <typename T>
class matrix
{
private:
    typedef std::vector<std::vector<T>> vec;
    void luDecomposition(const vec& A, vec& L, vec& U, std::vector<unsigned>& pivot);
    vec* data;
    // ---------------------------------------------------------------------------------------------
public:
    // Constructors and Destructor
    matrix ();
    matrix (unsigned R, unsigned C); // R - row count, C - column count
    matrix (unsigned R, unsigned C, const T V); // V - initial value
    ~matrix ();
    
    // Copy Constructor and Assignment Operator
    matrix (const matrix& m);
    matrix& operator= (const matrix& m);
    matrix& operator+= (const matrix& B);
    matrix& operator-= (const matrix& B);
    
    // Accessors
    T& operator() (unsigned r, unsigned c);
    T  operator() (unsigned r, unsigned c) const;

    // Common operations
    void fill (frodo f = elrand);
    void resize (unsigned R, unsigned C);
    void display ();

    // Calculate or do something
    T det (T eig = 0);
    void tr (const matrix& m); // A.tr(B) : transpose B assign to A matrix function

    // Friend Declarations
    template <typename U>
    friend matrix<U> operator+ (const matrix<U>& A, const matrix<U>& B);
    template <typename U>
    friend matrix<U> operator- (const matrix<U>& A, const matrix<U>& B);
    template <typename U>
    friend matrix<U> operator* (const matrix<U>& A, const matrix<U>& B);

    unsigned R,C, V;
};
// ---------------------------------------------------------------------------------------------
template <typename T>
matrix<T>::matrix ()
{
    R = 1; C = 1;
    data = new vec(R, std::vector<T>(C));
}
template <typename T>
matrix<T>::matrix (unsigned _R, unsigned _C) : R(_R), C(_C)
{
    data = new vec(R, std::vector<T>(C));
}
template <typename T>
matrix<T>::matrix (unsigned _R, unsigned _C, T _V) : R(_R), C(_C), V(_V)
{
    data = new vec(R, std::vector<T>(C,V));
}
template <typename T>
matrix<T>::~matrix()
{
    delete data; // is that required with std::vectors? dunno yet
}
template <typename T>
matrix<T>::matrix (const matrix& m) : R(m.R), C(m.C)
{
    data = new vec(R, std::vector<T>(C));
    for (unsigned i=0; i<R; ++i)
        for (unsigned j=0; j<C; ++j)
            data->at(i).at(j) = m(i,j);
}
// ---------------------------------------------------------------------------------------------
// Functions

template <typename T>
void matrix<T>::fill (frodo f)
{
    #pragma omp parallel for
    for (unsigned i=0; i<R; ++i)
        for (unsigned j=0; j<C; ++j)
            data->at(i).at(j) = f(i,j);
            
}
template <typename T>
void matrix<T>::resize (unsigned R, unsigned C)
{
    data->resize(R);
    #pragma omp parallel for
    for (unsigned i=0; i<R; ++i)
        data->at(i).resize(C);
    this->R = R;
    this->C = C;
}

template <typename T>
void matrix<T>::tr(const matrix& m)
{
    this->resize(m.C,m.R);
    #pragma omp parallel for
    for (unsigned i=0; i<R; ++i)
        for (unsigned j=0; j<C; ++j)
            data->at(i).at(j) = m(j,i);
}

template <typename T>
void matrix<T>::display() {
    for (unsigned i = 0; i < R; ++i) {
        for (unsigned j = 0; j < C; ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(3) << data->at(i).at(j);
        }
        std::cout << std::endl;
    }
}

template <typename T>
T matrix<T>::det(T eig) {
    if (R != C) {
        throw std::invalid_argument("To get determinant, matrix must be a square matrix");
    }

    // Create a copy of the matrix
    vec A = *data;

    if (eig != 0){
        for (unsigned i = 0; i < R; i++){
            A.at(i).at(i) -= eig;
        }
    }
    // Initialize matrices for LU decomposition
    vec L(R, std::vector<T>(C, 0));
    vec U = A;
    std::vector<unsigned> pivot(R);

    // Perform LU decomposition
    luDecomposition(A, L, U, pivot);

    // Calculate determinant
    T det = 1;
    for (unsigned i = 0; i < R; i++) {
        det *= U[i][i];
    }

    // Consider the effect of the eigenvalue
    if (eig != 0) {
        det *= eig;
    }
    
    return det;
}

template <typename T>
void matrix<T>::luDecomposition(const vec& A, vec& L, vec& U, std::vector<unsigned>& pivot) {
    unsigned n = A.size();

    for (unsigned i = 0; i < n; i++) {
        pivot[i] = i;
    }

    for (unsigned k = 0; k < n - 1; k++) {
        unsigned p = k;
        #pragma omp parallel for
        for (unsigned i = k + 1; i < n; i++) {
            if (std::abs(U[i][k]) > std::abs(U[p][k])) {
                p = i;
            }
        }

        if (U[p][k] == 0) {
            throw std::runtime_error("LU decomposition failed: Matrix is singular.");
        }

        std::swap(pivot[k], pivot[p]);
        #pragma omp parallel for
        for (unsigned i = k + 1; i < n; i++) {
            L[i][k] = U[i][k] / U[k][k];
            for (unsigned j = k + 1; j < n; j++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }
}





// ---------------------------------------------------------------------------------------------
// Operators
template <typename T>
matrix<T>& matrix<T>::operator= (const matrix& m)
{
    
    if (R != m.R){
        data->resize(m.R);
        R = m.R;
    }
    if (C != m.C){
        #pragma omp parallel for
        for (unsigned i=0; i<m.R; ++i)
            data->at(i).resize(m.C);
        C = m.C;
    }

    #pragma omp parallel for collapse(2)
    for (unsigned i=0; i<m.R; ++i)
        for (unsigned j=0; j<m.C; ++j)
            data->at(i).at(j) = m(i,j);
            
    return *this;
}

template <typename T>
inline T& matrix<T>::operator() (unsigned r, unsigned c)
{
    return data->at(r).at(c);
}

template <typename T>
inline T matrix<T>::operator() (unsigned r, unsigned c) const
{
    return data->at(r).at(c);
}

template <typename T>
matrix<T>& matrix<T>::operator+= (const matrix<T>& B)
{
    if (R != B.R || C != B.C) {
            throw std::invalid_argument("Matrices must have the same dimensions for subtraction");
        }

    #pragma omp parallel for collapse(2)
    for (unsigned i=0; i< R; ++i)
        for (unsigned j=0; j< C; ++j)
            data->at(i).at(j) += B(i,j);

    return *this;
}

template <typename T>
matrix<T>& matrix<T>::operator-= (const matrix<T>& B)
{
    if (R != B.R || C != B.C) {
            throw std::invalid_argument("Matrices must have the same dimensions for subtraction");
        }

    #pragma omp parallel for collapse(2)
    for (unsigned i=0; i<R; ++i)
        for (unsigned j=0; j<C; ++j)
            data->at(i).at(j) -= B(i,j);

    return *this;
}


template <typename T>
matrix<T> operator+ (const matrix<T>& A, const matrix<T>& B)
{
    if (A.R != B.R || A.C != B.C) {
            throw std::invalid_argument("Matrices must have the same dimensions for subtraction");
        }

    matrix<T> data(A.R, B.C);
    #pragma omp parallel for collapse(2)
    for (unsigned i=0; i<data.R; ++i)
        for (unsigned j=0; j<data.C; ++j)
            data(i,j) = A(i,j) + B(i,j);

    return data;
}

template <typename T>
matrix<T> operator- (const matrix<T>& A, const matrix<T>& B)
{
    if (A.R != B.R || A.C != B.C) {
            throw std::invalid_argument("Matrices must have the same dimensions for summation");
        }

    matrix<T> data(A.R, B.C);
    #pragma omp parallel for
    for (unsigned i=0; i<data.R; ++i)
        for (unsigned j=0; j<data.C; ++j)
            data(i,j) = A(i,j) - B(i,j);

    return data;
}

template <typename T>
matrix<T> operator* (const matrix<T>& A, const matrix<T>& B)
{
    if (A.C != B.R) {
        throw std::invalid_argument("Number of columns in the first matrix must be equal to the number of rows in the second matrix for multiplication");
    }

    matrix<T> data(A.R, B.C);
    #pragma omp parallel for
    for (unsigned i=0; i<data.R; ++i)
        for (unsigned j=0; j<data.C; ++j)
            for (unsigned k=0; k<A.C; ++k)
                data(i,j) += A(i,k) * B(k,j);

    return data;
}


#endif