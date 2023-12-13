#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <thread>
#include "my_spinbarrier.h"
#include "safemat.h"
using namespace std;

const int MSIZE = 8192;
const int WORKERS = 4;
SpinBarrier gandalf(WORKERS);


void update (int ID, SafeMat& mat, bool add = true)
{
    for (int i=0; i<MSIZE; ++i)
    {
        for (int j=0; j<MSIZE; ++j)
        {
            if (add) mat(i,j) += (ID+1)/10.;
            else     mat(i,j) -= (ID+1)/10.;
        }
    }
}


void show_vals (SafeMat& mat)
{
    printf ("* Values: ");
    printf ("[ ");
    for (int n=0; n<4; ++n)
    {
        unsigned r = (MSIZE/4 * n) + 4, c = (MSIZE/4 * (n+1)) - 3;
        printf ("%+.3f ", mat(r, c));
    }
    printf ("]\n");
}

void show_time (double t, int W)
{
    printf (" in %6.2f s", t);
    printf (" [t*N: %6.2f] ", t*W);
}


void worker (int ID, SafeMat& R, SafeMat& A, SafeMat& B, SafeMat& C)
{
    double ts, tf, tp;
    // TEST 1: implicit fill
    ts = clock();
    A.fill(elsin);
    tf = clock();
    tp = (tf-ts)/CLOCKS_PER_SEC;
    if (ID == 0)
    {
        printf ("* Thread fill   ");
        show_time(tp, WORKERS);
        show_vals(A);
    }
    gandalf.wait();
    // TEST 2: explicit update with data race
    int REPS = 32/WORKERS;
    ts = clock();
    for (int k=0; k<REPS; ++k)
    {
        if ((k+ID)%2) update(ID, A, true);
        else          update(ID, A, false);
    }
    tf = clock();
    tp = (tf-ts)/CLOCKS_PER_SEC;
    if (ID == 0)
    {
        printf ("* Thread update ");
        show_time(tp, WORKERS);
        show_vals(A);
    }
    gandalf.wait();
    // TEST 3: implicit sum
    //if (ID == 0) printf ("Preparing...\n");
    A = R; B = R;
    ts = clock();
    REPS = 8;
    for (int k=0; k<REPS; ++k)
    {
        C = A + B;
    }
    tf = clock();
    tp = (tf-ts)/CLOCKS_PER_SEC;
    if (ID == 0)
    {
        printf ("* Thread sum    ");
        show_time(tp, WORKERS);
        show_vals(C);
    }
    gandalf.wait();
    // TEST 4: implicit product
    if (ID == 0) printf ("Calculating product...\n");
    ts = clock();
    C = A * B;
    tf = clock();
    tp = (tf-ts)/CLOCKS_PER_SEC;
    if (ID == 0)
    {
        printf ("* Thread product");
        show_time(tp, WORKERS);
        show_vals(C);
    }
    gandalf.wait();
}

void run_threads ()
{
    // Shared data
    SafeMat R(MSIZE, MSIZE), A(MSIZE, MSIZE), B(MSIZE, MSIZE), C(MSIZE, MSIZE);

    // TEST 0: filling matrix in a single thread
    double ts = clock();
    R.fill(elsin);
    double tf = clock();
    double tp = (tf-ts)/CLOCKS_PER_SEC;
    printf ("* Serial fill   ");
    show_time(tp, 1);
    show_vals(R);

    thread tid[WORKERS];
    for (int i=0; i<WORKERS; ++i)
    {
        /* Create thread and execute worker */
        tid[i] = thread(worker, i, ref(R), ref(A), ref(B), ref(C));
    }
    /* Wait for all threads to end */
    for (int i=0; i<WORKERS; ++i)
        tid[i].join();
}

int main ()
{
    printf ("Matrix size: %6d | %2d threads\n",
            MSIZE, WORKERS);
	run_threads();
}

