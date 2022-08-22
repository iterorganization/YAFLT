// Spline data
#include <eqdsk_drsep6.65_data.hpp>
#include <fwp17_baryCent_data.hpp>

// Naive bicubic interpolation method
#include <bicubic.hpp>

// Alglib implementation
#include <interpolation.h>

#include <vector>
#include <stdio.h>
#include <stdlib.h> /*srand, rand*/
#include <time.h> /*For srand seed*/

int main(){
    // Prepare the naive interpolation.
    BICUBIC_INTERP *naive_interp = new BICUBIC_INTERP();
    int n_cols = baseline_r.size();
    int n_rows = baseline_z.size();

    std::vector<std::vector<double>> reshaped_Psi;
    std::vector<double> buffer;
    for (int i=0; i<n_rows; i++){
        buffer.clear();
        for (int j=0; j<n_cols; j++){
            buffer.push_back(baseline_psi[i * n_cols + j]);
        }
        // reshaped_Psi.insert(reshaped_Psi.begin(), buffer);
        reshaped_Psi.push_back(buffer);
    }
    naive_interp->prepareContainers();
    naive_interp->setArrays(baseline_r, baseline_z, reshaped_Psi);

    // Alglib interpolator
    alglib::real_1d_array R, Z, Psi;
    alglib::spline2dinterpolant alglib_interp;
    R.setcontent(baseline_r.size(), &(baseline_r[0]));
    Z.setcontent(baseline_z.size(), &(baseline_z[0]));
    Psi.setcontent(baseline_psi.size(), &(baseline_psi[0]));
    spline2dbuildbicubicv(R, baseline_r.size(), Z, baseline_z.size(),
                          Psi, 1, alglib_interp);

    // For Random prepare lower and upper boundaries
    double rmin, rmax, rdiff, zmin, zmax, zdiff;
    rmin = baseline_r[2];
    rmax = baseline_r[baseline_r.size() - 2];
    rdiff = rmax - rmin;
    zmin = baseline_z[2];
    zmax = baseline_z[baseline_z.size() - 2];
    zdiff = zmax - zmin;

    // Output variables
    double naive_fval, naive_fvaldx, naive_fvaldy;//, naive_fvaldxdy;
    double alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy;

    //Random check
    double buff_r, buff_z;

    // for (int i=0; i<10;i++){
    //     random_number(rmin, rmax, buff_r);
    //     random_number(zmin, zmax, buff_z);
    //     printf("%f %f\n", buff_r, buff_z);
    // }
    int N = fwp17_r.size();
    for (int i=0; i<N; i++){

        buff_r = fwp17_r[i];
        buff_z = fwp17_z[i];
        naive_interp->getValues(buff_r, buff_z, naive_fval, naive_fvaldx, naive_fvaldy);
        spline2ddiff(alglib_interp, buff_r, buff_z, alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy);

        printf("Point %f %f\n", buff_r, buff_z);
        printf("Naive  %f %f %f\n", naive_fval, naive_fvaldx, naive_fvaldy);
        printf("Alglib %f %f %f\n", alglib_fval, alglib_fvaldx, alglib_fvaldy);

    }
    return 0;
}