#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <omp.h>
#include <stdio.h>
#include <ctime>
#include <string>
#include <chrono>



/**
 * This class is what we're supposed to work with
*/
template <typename T>
class matrix
{
private:
    typedef std::vector<std::vector<T>> vec;
    vec* matrix_A; // pirma matrica
    vec* matrix_B; // antra matrica
    vec* matrix_C; // rezultatas
    long dim_A_x, dim_A_y, dim_B_x, dim_B_y; // dimensions A[x][y]
    int threads;
    std::chrono::duration<double> duration; // laikas

    /* Gets stuff */
    void get_stuff(vec* A, vec* B, vec* C){
        matrix_A = A;
        matrix_B = B;
        matrix_C = C;
        dim_A_x = matrix_A->size();
        dim_A_y = matrix_A->at(0).size();
        dim_B_x = matrix_B->size();
        dim_B_y = matrix_B->at(0).size();
        omp_set_num_threads(threads);
    }
    // ---------------------------------------------------------------------------------------------
public:
    /** Contructor
     * \param A A matrica
     * \param B B matrica
     * \param C C matrica (not initialised)
    */
    matrix(vec* A, vec* B, vec* C)
    {
        threads = omp_get_max_threads();
        get_stuff(A,B,C);
    }
    // ---------------------------------------------------------------------------------------------
    /** Constructor
     * \param A A matrica
     * \param B B matrica
     * \param C C matrica (not initialised)
     * \param thr threads count (default max)
    */
    matrix(vec* A, vec* B, vec* C,
        int thr)
    {
        threads = thr;
        get_stuff(A,B,C);
    }
    // ---------------------------------------------------------------------------------------------
    /* Printer (tikriausiai siai klasei nepriklausys nes negali printint nezinomu dalyku)*/
    void print(std::string which_one, long xx, long yy)
    {
        if (which_one == "sum"){
            printf("\nMatrcų A ir B suma C;\n\n");
            for (int y = 0; y < yy; ++y) {
                for (int x = 0; x < xx; ++x) {
                    std::cout << matrix_C->at(x).at(y) << " ";
                }
                std::cout << std::endl;
            }
        // } else if (which_one == "sub"){
        //     printf("\nMatrcų A ir B skirtumas C;\n\n");
        //     for (int y = 0; y < yy; ++y) {
        //         for (int x = 0; x < xx; ++x) {
        //             std::cout << matrix_C->at(x).at(y) << " ";
        //         }
        //         std::cout << std::endl;
        //     }
        }
        std::cout << "Elapsed time: " << duration.count() << " seconds\n\n";
    }
    // ---------------------------------------------------------------------------------------------
    /* Allocates */
    void allocate(){
        long xx = dim_A_x, yy = dim_A_y;
        if (matrix_C->size() > 0){
            matrix_C->resize(xx);
            for (long i = 0; i < yy; i++){
                matrix_C->at(i).resize(yy);
            }
            } else {
                matrix_C = new vec(xx, std::vector<T>(yy, 0));
            }
    }
    // ---------------------------------------------------------------------------------------------
    /* Sudeda matricas */
    void summation(){
        if (dim_A_x == dim_B_x && dim_A_y == dim_B_y){
            long xx = dim_A_x, yy = dim_A_y; // dimensions
            allocate();
            printf("Sumavimas prasideda:\n");
            auto start = std::chrono::high_resolution_clock::now();
            #pragma omp parallel for collapse(2)
            for (long x = 0; x < xx; ++x) {
                for (long y = 0; y < yy; ++y) {
                    matrix_C->at(x).at(y) = matrix_A->at(x).at(y) + matrix_B->at(x).at(y);
                }
            }
            duration = std::chrono::high_resolution_clock::now() - start;
            print("sum", xx, yy);
        } else {
            printf("\nNegalimas veiksmas, matricų dimensijos nesutampa!!!\n\n");
        }
    }
    // ---------------------------------------------------------------------------------------------
    /* Atima matricas */
    void subtraction(){
        if (dim_A_x == dim_B_x && dim_A_y == dim_B_y){
            long xx = dim_A_x, yy = dim_A_y; // dimensions
            allocate();
            printf("Atemimas prasideda:\n");
            auto start = std::chrono::high_resolution_clock::now();
            #pragma omp parallel for collapse(2)
            for (long x = 0; x < xx; ++x) {
                for (long y = 0; y < yy; ++y) {
                    matrix_C->at(x).at(y) = matrix_A->at(x).at(y) - matrix_B->at(x).at(y);
                }
            }
            duration = std::chrono::high_resolution_clock::now() - start;
            print("sub", xx, yy);
        } else {
            printf("\nNegalimas veiksmas, matricų dimensijos nesutampa!!!\n\n");
        }
    }
};

#endif