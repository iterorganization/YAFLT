// Spline data
#include <eqdsk_drsep6.65_data.hpp>
#include <fwp17_baryCent_data.hpp>

// Naive bicubic interpolation method
#include <bicubic.hpp>
#include <flt.hpp>

// Alglib implementation
#include <interpolation.h>

#include <vector>
#include <stdio.h>

int main(){
    FLT *obj = new FLT();
    // Set input data
    obj->setNDIM(baseline_r.size(), baseline_z.size());

    obj->setRARR(baseline_r);
    obj->setZARR(baseline_z);
    obj->setPSI(baseline_psi);

    // The following is useless...
    obj->setFARR(baseline_r);
    obj->setFPOL(baseline_r);
    obj->setVacuumFPOL(baseline_z[baseline_z.size() - 1]);

    obj->prepareInterpolation();
    obj->prepareThreadContainers(5);

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
    double naive_fval, naive_fvaldx, naive_fvaldy, naive_fvaldxdy;
    double alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy;

    // For bcart
    std::vector<double> buff_out;
    buff_out.resize(3);

    //Random check
    double buff_r, buff_z, buff_phi;

    // for (int i=0; i<10;i++){
    //     random_number(rmin, rmax, buff_r);
    //     random_number(zmin, zmax, buff_z);
    //     printf("%f %f\n", buff_r, buff_z);
    // }
    int N = fwp17_r.size();
    for (int i=0; i<N; i++){

        buff_r = fwp17_r[i];
        buff_z = fwp17_z[i];
        buff_phi = fwp17_phi[i];

        obj->debug_getValues(buff_r, buff_z, naive_fval, naive_fvaldx, naive_fvaldy, naive_fvaldxdy, 3);
        spline2ddiff(alglib_interp, buff_r, buff_z, alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy);

        printf("Point %f %f\n", buff_r, buff_z);
        printf("Naive  %f %f %f %f\n", naive_fval, naive_fvaldx, naive_fvaldy, naive_fvaldxdy);
        printf("Alglib %f %f %f %f\n", alglib_fval, alglib_fvaldx, alglib_fvaldy, alglib_fvaldxdy);

        // Now the BCart
        obj->getBCart(buff_r, buff_z, buff_phi, buff_out);
        printf("Naive %f %f %f\n", buff_out[0], buff_out[1], buff_out[2]);

    }
    return 0;
}