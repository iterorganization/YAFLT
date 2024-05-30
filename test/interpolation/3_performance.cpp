// Spline data
#include <data_for_spline.h>

// Naive bicubic interpolation method
#include <bicubic.hpp>

// Alglib implementation
#include <interpolation.h>

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

    // Alglib interpolator
    alglib::real_1d_array R, Z, Psi;
    alglib::spline2dinterpolant alglib_interp;
    R.setcontent(eq3_Rs.size(), &(eq3_Rs[0]));
    Z.setcontent(eq3_Zs.size(), &(eq3_Zs[0]));
    Psi.setcontent(eq3_Psi.size(), &(eq3_Psi[0]));
    spline2dbuildbicubicv(R, eq3_Rs.size(), Z, eq3_Zs.size(),
                          Psi, 1, alglib_interp);

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
    // OMP related variables
    int tid;
    int offset;

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

    begin = std::chrono::high_resolution_clock::now();
    for (int i=0; i<N; i++){
        spline2ddiff(alglib_interp, buff_r[i], buff_z[i], fval, fvaldx, fvaldy, fvaldxdy);
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Alglib interpolation time for %d iterations: %f\n", N, elapsed.count() * 1e-9);


    return 0;
}