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
    
    // Accessors
    T& operator() (unsigned r, unsigned c);
    T  operator() (unsigned r, unsigned c) const;

    // Common operations
    void fill (frodo f = elrand);
    void resize (unsigned R, unsigned C);

    // Friend Declarations
    template <typename U>
    friend matrix<U> operator+ (const matrix<U>& A, const matrix<U>& B);
    template <typename U>
    friend matrix<U> operator- (const matrix<U>& A, const matrix<U>& B);
    template <typename U>
    friend matrix<U> operator* (const matrix<U>& A, const matrix<U>& B);

    unsigned R,C, V;
};

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

// note: matrices must be of the same dimensions, else - error
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
matrix<T> operator+ (const matrix<T>& A, const matrix<T>& B)
{
    // if (A.R != B.R || A.C != B.C) {
    //         throw std::invalid_argument("Matrices must have the same dimensions for subtraction");
    //     }

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