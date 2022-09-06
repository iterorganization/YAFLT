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

#include <cmath>

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
    naive_interp->prepareContainers();
    naive_interp->setArrays(eq3_Rs, eq3_Zs, reshaped_Psi);

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
    double naive_fval, naive_fvaldx, naive_fvaldy;//, naive_fvaldxdy;
    double alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy;

    //Random check
    double buff_x, buff_y;

    // for (int i=0; i<10;i++){
    //     random_number(rmin, rmax, buff_x);
    //     random_number(zmin, zmax, buff_y);
    //     printf("%f %f\n", buff_x, buff_y);
    // }
    int N = 1000*1000;
    printf("Comparing %d values\n", N);
    bool print;
    int counts=0;
    double abs_max_psi=-1, abs_max_psidx=-1, abs_max_psidy=-1, buff;
    for (int i=0; i<N; i++){
        print=false;
        buff_x = rmin + rdiff * rand() / RAND_MAX;
        buff_y = zmin + zdiff * rand() / RAND_MAX;
        naive_interp->getValues(buff_x, buff_y, naive_fval, naive_fvaldx, naive_fvaldy);
        spline2ddiff(alglib_interp, buff_x, buff_y, alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy);
        buff = std::fabs(naive_fval - alglib_fval) / alglib_fval;
        if (buff > 0.03) {
            print=true;
        }
        if (buff > abs_max_psi) {
            abs_max_psi = buff;
        }

        buff = std::fabs(naive_fvaldx - alglib_fvaldx) / alglib_fvaldx;
        if (buff > 0.03) {
            print=true;
        }
        if (buff > abs_max_psidx) {
            abs_max_psidx = buff;
        }

        buff = std::fabs(naive_fvaldy - alglib_fvaldy) / alglib_fvaldy;
        if (buff > 0.03) {
            print=true;
        }
        if (buff > abs_max_psidy) {
            abs_max_psidy = buff;
        }

        if (print) {
            counts = counts + 1;
            printf("Point %f %f\n", buff_x, buff_y);
            printf("Naive  %f %f %f\n", naive_fval, naive_fvaldx, naive_fvaldy);
            printf("Alglib %f %f %f\n", alglib_fval, alglib_fvaldx, alglib_fvaldy);
        }
    }
    printf("Total tries: %d\n", N);
    printf("Number of high error points: %d\n", counts);
    printf("Maximum relative error psi: %f\n", abs_max_psi);
    printf("Maximum relative error psidx: %f\n", abs_max_psidx);
    printf("Maximum relative error psidy: %f\n", abs_max_psidy);
    return 0;
}