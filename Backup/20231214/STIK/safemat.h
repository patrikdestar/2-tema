#ifndef SAFEMAT_H
#define SAFEMAT_H

#include <cstdlib>
#include <cmath>


double elrand (unsigned /*i*/, unsigned /*j*/)
{
    return 10*(rand()/(double)RAND_MAX);
}

double elsin (unsigned i, unsigned j)
{
    return sin((i+1)*(j+2));
}

using elf = double (*) (unsigned, unsigned);


class SafeMat
{
    public:
    // C/D
    SafeMat (unsigned M, unsigned N);
    ~SafeMat ();
    SafeMat (const SafeMat& m);
    SafeMat& operator= (const SafeMat& m);
    // Direct access
    double& operator() (unsigned r, unsigned c);
    double  operator() (unsigned r, unsigned c) const;
    // Common operations
    void fill (elf f = elrand);
    friend SafeMat operator+ (const SafeMat& A, const SafeMat& B);
    friend SafeMat operator* (const SafeMat& A, const SafeMat& B);

//  private:
    double** data;
    unsigned M, N;
};

SafeMat::SafeMat (unsigned _M, unsigned _N) : M(_M), N(_N)
{
    data = new double*[M];
    double* mem = new double[M*N];
    for (unsigned i=0; i<M; ++i)
        data[i] = &mem[i*N];
}

SafeMat::~SafeMat ()
{
    delete[] data[0];
    delete[] data;
}

SafeMat::SafeMat (const SafeMat& A): M(A.M), N(A.N)
{
    data = new double*[M];
    double* mem = new double[M*N];
    for (unsigned i=0; i<M; ++i)
        data[i] = &mem[i*N];

    for (unsigned i=0; i<M; ++i)
        for (unsigned j=0; j<N; ++j)
            data[i][j] = A(i,j);
}

SafeMat& SafeMat::operator= (const SafeMat& A)
{
    for (unsigned i=0; i<M; ++i)
        for (unsigned j=0; j<N; ++j)
            data[i][j] = A(i,j);
    return *this;
}

inline
double& SafeMat::operator() (unsigned r, unsigned c)
{
    return data[r][c];
}
inline
double SafeMat::operator() (unsigned r, unsigned c) const
{
    return data[r][c];
}

void SafeMat::fill (elf f)
{
    for (unsigned i=0; i<M; ++i)
        for (unsigned j=0; j<N; ++j)
            data[i][j] = f(i,j);
}

SafeMat operator+ (const SafeMat& A, const SafeMat& B)
{
    SafeMat data(A.M, A.N);
    for (unsigned i=0; i<data.M; ++i)
        for (unsigned j=0; j<data.N; ++j)
            data(i,j) = A(i,j) + B(i,j);
    return data;
}

SafeMat operator* (const SafeMat& A, const SafeMat& B)
{
    SafeMat data(A.M, B.N);
    for (unsigned i=0; i<data.M; ++i)
        for (unsigned j=0; j<data.N; ++j)
            for (unsigned k=0; k<A.N; ++k)
                data(i,j) += A(i,k) * B(k,j);
    return data;
}

#endif // SAFEMAT_H
