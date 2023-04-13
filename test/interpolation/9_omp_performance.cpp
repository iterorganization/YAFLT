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
#define THREAD_NUM 4

int main(){
    omp_set_num_threads(THREAD_NUM);
    srand(time(NULL));
    // Prepare the naive interpolation.
    BICUBIC_INTERP *naive_interp = new BICUBIC_INTERP();
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
    naive_interp->prepareContainers(THREAD_NUM);
    naive_interp->setArrays(eq3_Rs, eq3_Zs, reshaped_Psi);

    // Now prepare THREAD_NUM interpolators, to see if having multiple
    // instances is faster...
    std::vector<BICUBIC_INTERP *> naive_interpolators;
    naive_interpolators.resize(THREAD_NUM);
    for (int i=0; i<THREAD_NUM; i++){
        BICUBIC_INTERP *interp = new BICUBIC_INTERP();
        naive_interpolators[i] = interp;

        // Populate with data
        naive_interpolators[i]->prepareContainers(1); // Only one!
        naive_interpolators[i]->setArrays(eq3_Rs, eq3_Zs, reshaped_Psi);
    }


    // For Random prepare lower and upper boundaries
    double rmin, rmax, rdiff, zmin, zmax, zdiff;
    rmin = eq3_Rs[2];
    rmax = eq3_Rs[eq3_Rs.size() - 2];
    rdiff = rmax - rmin;
    zmin = eq3_Zs[2];
    zmax = eq3_Zs[eq3_Zs.size() - 2];
    zdiff = zmax - zmin;

    // Output variables
    double fval, fvaldx, fvaldy, fvaldxdy;

    // for (int i=0; i<10;i++){
    //     random_number(rmin, rmax, buff_x);
    //     random_number(zmin, zmax, buff_y);
    //     printf("%f %f\n", buff_x, buff_y);
    // }
    int N = 100*1000*1000;
    int chunk_N = N / THREAD_NUM;
    // OMP related variables
    int tid;
    int offset;
    int index;


    std::vector<double> buff_r;
    std::vector<double> buff_z;

    buff_r.resize(N);
    buff_z.resize(N);


    // Let's just make a large array and put the values there...
    printf("Creating %d random points.\n", N);
    auto begin = std::chrono::high_resolution_clock::now();
    for (int i=0; i<N; i++){
        buff_r[i] = rmin + rdiff * rand() / RAND_MAX;
        buff_z[i] = zmin + zdiff * rand() / RAND_MAX;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Completed in: %f seconds\n", elapsed.count() * 1e-9);

    // SERIAL
    begin = std::chrono::high_resolution_clock::now();
    tid = 0;
    for (int i=0; i<N; i++){
        naive_interp->getAllValues(buff_r[i], buff_z[i], fval, fvaldx, fvaldy, fvaldxdy);
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Serial interpolation time for %d iterations: %f\n", N, elapsed.count() * 1e-9);


    // DEPRECATED
    // // OMP WITH ONE OBJECT
    // begin = std::chrono::high_resolution_clock::now();
    // #pragma omp parallel private(tid, offset, index, fval, fvaldx, fvaldy, fvaldxdy)
    // {
    //     tid = omp_get_thread_num();
    //     offset = tid * chunk_N;
    //     for (int i=0; i<chunk_N; i++){
    //         index = offset + i;
    //         naive_interp->getAllValues(buff_r[index], buff_z[index], fval, fvaldx, fvaldy, fvaldxdy, tid);
    //     }
    // }

    // end = std::chrono::high_resolution_clock::now();
    // elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    // printf("Parallel interpolation time for %d iterations with %d threads: %f\n", N, THREAD_NUM, elapsed.count() * 1e-9);

    // OMP WITH THREAD_NUM OBJECTS.
    begin = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid, offset, index, fval, fvaldx, fvaldy, fvaldxdy)
    {
        tid = omp_get_thread_num();
        offset = tid * chunk_N;
        for (int i=0; i<chunk_N; i++){
            index = offset + i;
            naive_interpolators[tid]->getAllValues(buff_r[index], buff_z[index], fval, fvaldx, fvaldy, fvaldxdy);
        }
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Parallel interpolation time for %d iterations with %d threads and %d separate interpolator objects: %f\n", N, THREAD_NUM, THREAD_NUM, elapsed.count() * 1e-9);


    return 0;
}