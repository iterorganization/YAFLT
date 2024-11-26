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

#include <random>

int main(){
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
    naive_interp->setArrays(eq3_Rs, eq3_Zs, reshaped_Psi);

    BI_DATA *context = new BI_DATA();

    naive_interp->populateContext(context);

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


    begin = std::chrono::high_resolution_clock::now();
    for (int i=0; i<N; i++){
        context->r = buff_r[i];
        context->z = buff_z[i];

        naive_interp->getValues(context);
        fval = context->val;
        fvaldx = context->valdx;
        fvaldy = context->valdy;
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Naive interpolation time for %d iterations: %f\n", N, elapsed.count() * 1e-9);
    printf("Or number of evaluations per second: %.2f\ M/s \n", N * 1e3 / (elapsed.count()));
    return 0;
}