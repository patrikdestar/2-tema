#include <iostream>
#include <stdio.h>
#include <iomanip>
#include "matrix.h"
#include "Complex.h"

const int MSIZE = 8192;
const int MSIZE2 = 1024;
int threads;
typedef double type;

void show_time (double t, int W)
{
    printf (" in %6.3f s", t);
    printf (" [t*N: %6.3f] ", t*W);
}
void show_vals (matrix<type>& mat, int i = 0)
{
    unsigned size = (i == 0) ? MSIZE : MSIZE2;
    printf ("* Values: ");
    printf ("[ ");
    for (int n=0; n<4; ++n)
    {
        unsigned r = (size/4 * n) + 4, c = (size/4 * (n+1)) - 3;
        printf ("%+.3f ", mat(r, c));
    }
    printf ("]\n");
}

void run()
{
    // My precious matrces
    matrix<type> R(MSIZE, MSIZE), A(MSIZE, MSIZE), B(MSIZE, MSIZE), C(MSIZE, MSIZE);
    std::cout << R.C;
    for (int a = 0; a < 4; a++){
    if (a == 0) threads = 1;
    else threads = 4*a;
    printf("\nSize:   %i | %i threads\n", MSIZE, threads);
    // TEST 0: filling matrix in a single thread
    omp_set_num_threads(1);
    auto start = std::chrono::high_resolution_clock::now();
    R.fill(baggsin);
    auto duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Serial fill   ");
    show_time(((double)duration.count())/1000000000, 1);
    show_vals(R, 0);

    omp_set_num_threads(threads);
    // TEST 1: implicit fill
    start = std::chrono::high_resolution_clock::now();
    A.fill(baggsin);
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread fill   ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 0);

    // TEST 2: implicit sum
    A = R; B = R;
    start = std::chrono::high_resolution_clock::now();
    int REPS = 1;
    for (int k=0; k<REPS; ++k)
    {
        C = A + B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sum    ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 0);

    }

    C.resize(MSIZE2,MSIZE2);A.resize(MSIZE2,MSIZE2);B.resize(MSIZE2,MSIZE2);
    for (int a = 0; a < 4; a++){
    if (a == 0) threads = 1;
    else threads = 4*a;
    omp_set_num_threads(threads);
    printf("\nSize:   %i | %i threads\n", MSIZE2, threads);
    // TEST 3: implicit product
    printf ("Calculating product...\n");
    auto start = std::chrono::high_resolution_clock::now();
    C = A * B;

    auto duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread product");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 1);
    }
}

int main(){
    printf("Hello world\n");
    omp_set_num_threads(2);
    // Pastabos so far:
    // reikia sukurti klase kuri ir yra matrica, ir veiksmai su ta klase
    // turi but aprašomi zenklais

    run();
    
    return 0;
}