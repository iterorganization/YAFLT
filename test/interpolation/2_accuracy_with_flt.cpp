// Spline data
// #include <data_for_spline.h>
#include <eqdsk_drsep10_data.hpp>

// Naive bicubic interpolation method
#include <bicubic.hpp>
#include <flt.hpp>

// Alglib implementation
#include <interpolation.h>

#include <vector>
#include <stdio.h>
#include <stdlib.h> /*srand, rand*/
#include <time.h> /*For srand seed*/

#include <cmath>

int main(){
    srand(time(NULL));
    FLT *obj = new FLT();
    // Set input data
    obj->setNDIM(drsep10_Rs.size(), drsep10_Zs.size());

    obj->setRARR(drsep10_Rs);
    obj->setZARR(drsep10_Zs);
    obj->setPSI(drsep10_Psi);

    // The following is useless...
    obj->setFARR(drsep10_Rs);
    obj->setFPOL(drsep10_Rs);
    obj->setVacuumFPOL(drsep10_Zs[drsep10_Zs.size() - 1]);

    obj->prepareInterpolation();
    obj->prepareThreadContainers(5);

    // Alglib interpolator
    alglib::real_1d_array R, Z, Psi;
    alglib::spline2dinterpolant alglib_interp;
    R.setcontent(drsep10_Rs.size(), &(drsep10_Rs[0]));
    Z.setcontent(drsep10_Zs.size(), &(drsep10_Zs[0]));
    Psi.setcontent(drsep10_Psi.size(), &(drsep10_Psi[0]));
    spline2dbuildbicubicv(R, drsep10_Rs.size(), Z, drsep10_Zs.size(),
                          Psi, 1, alglib_interp);

    // For Random prepare lower and upper boundaries
    double rmin, rmax, rdiff, zmin, zmax, zdiff;
    rmin = drsep10_Rs[1];
    rmax = drsep10_Rs[drsep10_Rs.size() - 1];
    rdiff = rmax - rmin;
    zmin = drsep10_Zs[1];
    zmax = drsep10_Zs[drsep10_Zs.size() - 1];
    zdiff = zmax - zmin;

    // Output variables
    double naive_fval, naive_fvaldx, naive_fvaldy;//, naive_fvaldxdy;
    double alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy;

    //Random check
    double buff_x, buff_y;

    int N = 100*1000*1000;
    printf("Comparing %d values", N);
    bool print;
    int counts=0;
    double diff=0.01;
    double abs_max_psi=-1, abs_max_psidx=-1, abs_max_psidy=-1, buff;
    for (int i=0; i<N; i++){
        print=false;
        buff_x = rmin + rdiff * rand() / RAND_MAX;
        buff_y = zmin + zdiff * rand() / RAND_MAX;
        obj->debug_getValues(buff_x, buff_y, naive_fval, naive_fvaldx, naive_fvaldy, 3);
        spline2ddiff(alglib_interp, buff_x, buff_y, alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy);
        buff = std::fabs(naive_fval - alglib_fval) / alglib_fval;
        if (buff > diff) {
            print=true;
        }
        if (buff > abs_max_psi) {
            abs_max_psi = buff;
        }

        buff = std::fabs(naive_fvaldx - alglib_fvaldx) / alglib_fvaldx;
        if (buff > diff) {
            print=true;
        }
        if (buff > abs_max_psidx) {
            abs_max_psidx = buff;
        }

        buff = std::fabs(naive_fvaldy - alglib_fvaldy) / alglib_fvaldy;
        if (buff > diff) {
            print=true;
        }
        if (buff > abs_max_psidy) {
            abs_max_psidy = buff;
        }

        // Now check if it is inside the vessell!
        if (print && 4.0 <= buff_x && buff_x <= 8.5 && -4.3 <= buff_y && buff_y <= 4.8) {
            counts = counts + 1;
            printf("Point %f %f\n", buff_x, buff_y);
            printf("Naive  %f %f %f\n", naive_fval, naive_fvaldx, naive_fvaldy);
            printf("Alglib %f %f %f\n", alglib_fval, alglib_fvaldx, alglib_fvaldy);
        }
    }
    printf("Total tries: %d\n", N);
    printf("Number of high error points: %d\n", counts);
    printf("Maximum error psi: %f\n", abs_max_psi);
    printf("Maximum error psidx: %f\n", abs_max_psidx);
    printf("Maximum error psidy: %f\n", abs_max_psidy);
    return 0;
}