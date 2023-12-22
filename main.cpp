#include <iostream>
#include <stdio.h>
#include <iomanip>
#include "matrix.h"
#include "Complex.h"

const int MSIZE = 8192;
const int MSIZE2 = 1024;
int threads = 6;
typedef double type;

void show_time (double t, int W)
{
    printf (" in %6.3f s", t);
    printf (" [t*N: %6.3f] ", t*W);
}
void show_vals (matrix<type>& mat, int i = 0)
{
    unsigned size = mat.R;
    printf ("* Values: ");
    printf ("[ ");
    for (int n=0; n<4; ++n)
    {
        unsigned r = (size/4 * n) + 4, c = (size/4 * (n+1)) - 3;
        printf ("%+.3f ", mat(r, c));
    }
    printf ("]\n");
}

void run_threads()
{
    // My precious matrces
    matrix<type> R(MSIZE, MSIZE), A(MSIZE, MSIZE), B(MSIZE, MSIZE), C(MSIZE, MSIZE);
    for (int a = 0; a < 4; a++){
    if (a == 0) threads = 1;
    else threads = 4*a;
    printf("\nSize:   %i | %i threads\n", MSIZE, threads);
    // TEST 0: filling matrix in a single thread
    omp_set_num_threads(1);
    auto start = std::chrono::high_resolution_clock::now();
    R.fill(baggsin);
    auto duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Serial fill      ");
    show_time(((double)duration.count())/1000000000, 1);
    show_vals(R, 0);

    omp_set_num_threads(threads);
    // TEST 1: implicit fill
    start = std::chrono::high_resolution_clock::now();
    B.fill(baggsin);
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread fill      ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(B, 0);

    // TEST 2: implicit sum
    A.tr(R);
    start = std::chrono::high_resolution_clock::now();
    int REPS = 1;
    for (int k=0; k<REPS; ++k)
    {
        C = A + B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sum (+)   ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 0);

    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        A += B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sum (+=)  ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 0);

    // TEST 3: implicit sub
    A.tr(R);
    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        C = A - B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sub (-)   ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 0);

    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        A -= B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sub (-=)  ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 0);
    A.tr(R);

    // TEST 4: tr
    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        A.tr(C);
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread tr        ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 0);
    std::cout << "\n";
    }

    C.resize(MSIZE2,MSIZE2);A.resize(MSIZE2,MSIZE2);B.resize(MSIZE2,MSIZE2);R.resize(MSIZE2,MSIZE2);
    unsigned detn=MSIZE2/100;
    matrix<type> G(detn, detn, 1);
    for (int a = 0; a < 4; a++){
    if (a == 0) threads = 1;
    else threads = 4*a;
    omp_set_num_threads(threads);
    printf("\nSize:   %i | %i threads\n", MSIZE2, threads);
    // TEST 4: det


    printf ("Calculating det (%uX%u matrix of ones, eig = 1)...\n", detn, detn);
    auto start = std::chrono::high_resolution_clock::now();
    type det = G.det(1);
    auto duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread det    ");
    show_time(((double)duration.count())/1000000000, threads);
    printf ("* Value: [ %+.3f ]\n", det);


    A.tr(R);
    // TEST 6: implicit product
    printf ("\nCalculating product...\n");
    start = std::chrono::high_resolution_clock::now();
    C = A * B;
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread product (*)");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 1);

    printf ("\nCalculating product...\n");
    start = std::chrono::high_resolution_clock::now();
    A *= B;
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread product (*=)");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 1);
    }
}

void run_size()
{
    unsigned sizes[5] = {1024,2048, 4086, 8192, 16384};
    // My precious matrces
    matrix<type> R(MSIZE, MSIZE), A(MSIZE, MSIZE), B(MSIZE, MSIZE), C(MSIZE, MSIZE);
    for (int a = 0; a < 5; a++){
    int size = sizes[a];
    C.resize(size,size);A.resize(size,size);B.resize(size,size);R.resize(size,size);
    printf("\nSize:   %i | %i threads\n", size, threads);
    R.fill(baggsin);

    omp_set_num_threads(threads);
    // TEST 1: implicit fill
    auto start = std::chrono::high_resolution_clock::now();
    B.fill(baggsin);
    auto duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread fill      ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(B, 0);

    // TEST 2: implicit sum
    A.tr(R);
    start = std::chrono::high_resolution_clock::now();
    int REPS = 1;
    for (int k=0; k<REPS; ++k)
    {
        C = A + B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sum (+)   ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 0);

    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        A += B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sum (+=)  ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 0);

    // TEST 3: implicit sub
    A.tr(R);
    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        C = A - B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sub (-)   ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 0);

    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        A -= B;
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread sub (-=)  ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 0);
    A.tr(R);

    // TEST 4: tr
    start = std::chrono::high_resolution_clock::now();
    for (int k=0; k<REPS; ++k)
    {
        A.tr(C);
    }
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread tr        ");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 0);
    std::cout << "\n";
    }
    
    unsigned sizes_d[5] = {128, 256, 512, 1024, 2048};
    unsigned detn;
    matrix<type> G;
    for (int a = 0; a < 5; a++){
        omp_set_num_threads(12);
        threads = 12;

    unsigned size = sizes_d[a];
    C.resize(size,size);A.resize(size,size);B.resize(size,size);R.resize(size,size);
    printf("\nSize:   %i | %i threads\n", size, threads);
    // TEST 4: det

    C.fill();R.fill(baggsin);B.fill(baggsin);A.tr(R);
    detn=size%100;
    G.resize(detn,detn);
    G.fill(1);
    printf ("Calculating det (%uX%u matrix of ones, eig = 1)...\n", detn, detn);
    auto start = std::chrono::high_resolution_clock::now();
    type det = G.det(1);
    auto duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread det    ");
    show_time(((double)duration.count())/1000000000, threads);
    printf ("* Value: [ %+.3f ]\n", det);
    
    // TEST 6: implicit product
    printf ("\nCalculating product...\n");
    start = std::chrono::high_resolution_clock::now();
    C = A * B;
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread product (*)");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(C, 1);

    printf ("\nCalculating product...\n");
    start = std::chrono::high_resolution_clock::now();
    A *= B;
    duration = std::chrono::high_resolution_clock::now() - start;
    printf ("* Thread product (*=)");
    show_time(((double)duration.count())/1000000000, threads);
    show_vals(A, 1);
    }
}

int main(){
    // Main testas
    if (true){run_threads();
    run_size();}

    // inline operatoriu testas
    if (false){
        omp_set_num_threads(4);
        matrix<type> A(2,2,0);
        #pragma omp parallel
        {
            A(0,0) += 1+1;
        }
        A.display();
    }

    // basic pristatymo testas
    if (false){
        omp_set_num_threads(8);
        matrix<type> A(2,2,2);
        printf("\ndeterminants: %+.3f\n",A.det(1));
        A.display();

        matrix<type> B(3,3,2);
        printf("\ndeterminants: %+.3f\n",B.det(1));
        B.display();

        // matrix<type> C;
        // printf("\n");
        // C.display();

        // // printf("\nC = A + B:\n");
        // // C = A + B;
        // // C.display();

        // A.resize(3, 3);
        // A.fill((type)1);
        // printf("\nC = A + B:\n");
        // C = A + B;
        // C.display();

        // printf("\nC = A - B:\n");
        // C = A - B;
        // C.display();

        // printf("\nC = A * B:\n");
        // C = A * B;
        // C.display();

        // A.resize(6, 3);
        // A.fill((type)1);
        // printf("\nC = A * B:\n");
        // C = A * B;
        // C.display();
    }

    return 0;
}
