// Spline data
#include <data_for_spline.h>

// Naive bicubic interpolation method
#include <bicubic.hpp>

#include <vector>
#include <stdio.h>
#include <stdlib.h> /*srand, rand*/
#include <time.h> /*For srand seed*/

// Measure time
#include <chrono>

// OMP
#include <omp.h>
#include <random>

// Thread for obtaining max number of CPUs
#include <thread>

int main(){
    int THREAD_NUM;
    THREAD_NUM = std::min(32, static_cast<int>(std::thread::hardware_concurrency()));

    omp_set_num_threads(THREAD_NUM);
    srand(time(NULL));
    // Prepare the naive interpolation.
    int n_cols = eq3_Rs.size();
    int n_rows = eq3_Zs.size();

    std::vector<std::vector<double>> reshaped_Psi;
    std::vector<double> buffer;
    for (int i=0; i<n_rows; i++){
        buffer.clear();
        for (int j=0; j<n_cols; j++){
            buffer.push_back(eq3_Psi[i * n_cols + j]);
        }
        // reshaped_Psi.insert(reshaped_Psi.begin(), buffer);
        reshaped_Psi.push_back(buffer);
    }

    // For Random prepare lower and upper boundaries
    double rmin, rmax, rdiff, zmin, zmax, zdiff;
    rmin = eq3_Rs[2];
    rmax = eq3_Rs[eq3_Rs.size() - 2];
    rdiff = rmax - rmin;
    zmin = eq3_Zs[2];
    zmax = eq3_Zs[eq3_Zs.size() - 2];
    zdiff = zmax - zmin;

    // for (int i=0; i<10;i++){
    //     random_number(rmin, rmax, buff_x);
    //     random_number(zmin, zmax, buff_y);
    //     printf("%f %f\n", buff_x, buff_y);
    // }
    int N = 1000*1000*1000;

    std::vector<double> buff_r;
    std::vector<double> buff_z;

    buff_r.resize(N);
    buff_z.resize(N);

    // Interpolation object

    BICUBIC_INTERP *bicubic_interp = new BICUBIC_INTERP();
    bicubic_interp->setArrays(eq3_Rs, eq3_Zs, reshaped_Psi);
    BI_DATA *context= new BI_DATA();
    bicubic_interp->populateContext(context);

    // Let's just make a large array and put the values there...
    printf("Creating %d random points.\n", N);

    double begin, end, elapsed, serial_elapsed;

    begin = omp_get_wtime();
    for (int i=0; i<N; i++){
        buff_r[i] = rmin + rdiff * rand() / RAND_MAX;
        buff_z[i] = zmin + zdiff * rand() / RAND_MAX;
    }
    end = omp_get_wtime();
    elapsed = end-begin;
    printf("Completed in: %f seconds\n", elapsed);

    // SERIAL
    // begin = std::chrono::high_resolution_clock::now();

    begin = omp_get_wtime();
    for (int i=0; i<N; i++){
        context->r = buff_r[i];
        context->z = buff_z[i];
        bicubic_interp->getValues(context);
        // bicubic_interp->dummy(context);
    }
    end = omp_get_wtime();
    elapsed = end-begin;
    serial_elapsed = elapsed;
    printf("Serial interpolation time for %d iterations: %f\n", N, elapsed);

    // OMP WITH ONE OBJECT
    // begin = std::chrono::high_resolution_clock::now();
    begin = omp_get_wtime();

    #pragma omp parallel
    {
        //Class that holds all the data as copies.
        BI_DATA *omp_context = new BI_DATA();
        bicubic_interp->populateContext(omp_context);
        // double *r = &buff_r[0];
        // double *z = &buff_z[0];
        #pragma omp for
        for (int i=0; i<N; i++){
            omp_context->r = buff_r[i];
            omp_context->z = buff_z[i];
            bicubic_interp->getValues(omp_context);
        }
        // This is still on the parallel block
    }
    end = omp_get_wtime();
    elapsed = end-begin;
    printf("Parallel interpolation time for %d iterations with %d threads: %f\n", N, THREAD_NUM, elapsed);
    printf("Speedup: %.2f\n", serial_elapsed / elapsed);
}