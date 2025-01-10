/*
  Software for testing RKF4(5) accuracy by tracing a line on a contour in a
  cylindrical coordinate system.

  For input we take a simple analytical function:

  f(x, y) = 0.5*(x - 3)**2 + 2.3*(y - 5)**2

  The gradient:

  df/dx = 0.5*2(x-5)
  df/dy = 2.3*2(y-5)

  The pde we actually solve is:
  ds/dx = -df/dy
  ds/dy =  df/dx


 The problem is time independant. That is, the field does not change over time.

 */

// PDE derivative function

#include <cmath>
#include <cstdio>

double function(double y[2]){
    return 0.5*(y[0] - 3.0) * (y[0] - 3.0) + 2.3*(y[1] - 5.0) * (y[1] - 5.0);
}

void pde(double y[2], double yp[2]){

    yp[0] = -2.3*2.0 * (y[1] - 5.0);
    yp[1] =  0.5*2.0 * (y[0] - 3.0);
}

// Fehlberg step

void fehl_step(double y[2], double h, double yp[2],
               double k2[2], double k3[2], double k4[2], double k5[2],
               double k6[2], double s[2]){
    // Calculate the k1, k2, k3, k4, k5, k6 factors of the RKF4(5) method and
    // the approximate solution. The k* factors are used for error estimation
    // while the s holds the approximate solution. If the error is acceptable
    // then the approximate solution is the new solution.

    // The factors are calculated based on Table II in Fehlberg, NASA Technical
    // Report 287.

    // In this case we have a non-time dependant PDE

    // Use k6 as container for  calculating other k factors (used as argument
    // for the pde method).
    double ch;

    // k1 is h pde(t, y)

    // k2
    ch = h * 0.25;
    k6[0] = y[0] + ch * yp[0];
    k6[1] = y[1] + ch * yp[1];

    // Evaluate the function for k2
    // pde(t + ch, k6, k2);
    pde(k6, k2);

    // k3
    ch = 3.0 * h / 32.0;
    k6[0] = y[0] + ch * (yp[0] + 3.0 * k2[0]);
    k6[1] = y[1] + ch * (yp[1] + 3.0 * k2[1]);

    // Evaluate the function for k3
    // pde(t + 3.0*h/8.0, k6, k3);
    pde(k6, k3);

    // k4
    ch = h / 2197.0;
    k6[0] = y[0] + ch * (1932.0 * yp[0] + (7296.0 * k3[0] - 7200.0 * k2[0]));
    k6[1] = y[1] + ch * (1932.0 * yp[1] + (7296.0 * k3[1] - 7200.0 * k2[1]));

    // Evaluate the function for k5
    // pde(t + 12.0 * h / 13.0, k6, k4);
    pde(k6, k4);

    // k5
    ch = h / 4104.0;
    k6[0] = y[0] + ch * ((8341.0 * yp[0] - 845.0 * k4[0]) + (29440.0 * k3[0] - 32832.0 * k2[0]));
    k6[1] = y[1] + ch * ((8341.0 * yp[1] - 845.0 * k4[1]) + (29440.0 * k3[1] - 32832.0 * k2[1]));

    // Evaluate the function at the new position
    // pde(t + h, k6, k5);
    pde(k6, k5);

    // Finally k6
    // The k2 factor is not used in the last evaluations, so we will use k2
    // here for storing result for k6.
    ch = h / 20520.0;
    k2[0] = y[0] + ch * (-6080.0 * yp[0] + (9295.0 * k4[0] - 5643.0 * k5[0]) + (41040.0 * k2[0] - 28352.0 * k3[0]));
    k2[1] = y[1] + ch * (-6080.0 * yp[1] + (9295.0 * k4[1] - 5643.0 * k5[1]) + (41040.0 * k2[1] - 28352.0 * k3[1]));

    // Evaluate the function for k6
    // pde(t + 0.5 * h, k2, k6);
    pde(k2, k6);

    // Approximate solution using up to order 6
    ch = h / 7618050.0;
    s[0] = y[0] + ch * (902880.0 * yp[0] + (3855735.0 * k4[0] - 1371249.0 * k5[0]) + (3953664.0 * k3[0] + 277020.0 * k6[0]));
    s[1] = y[1] + ch * (902880.0 * yp[1] + (3855735.0 * k4[1] - 1371249.0 * k5[1]) + (3953664.0 * k3[1] + 277020.0 * k6[1]));
}

int main()
{
    /* code */
    double relerr = 1e-4;
    double abserr = 1e-4;
    double dt = 0.01;

    double y[2], yp[2], orig_contour_val, contour_rel_error, max_error;

    max_error = 0.0;
    y[0] = 2.0;
    y[1] = 4.0;
    pde(y, yp); // Initial derivative value

    // Value of the contour we will compare too
    orig_contour_val = function(y);

    // RKF4(5) variables
    double h, h_scale; // Actual step performed if function starts to act rowdy.

    // Step increment arrays
    double k2[2], k3[2], k4[2], k5[2], k6[2];
    // Solution array
    double s[2];
    // Error values
    const double error_scale = 2.0/relerr;
    const double error_absolute = error_scale * abserr;
    double error_et;
    double error_ete;
    double error_eot;
    double error_tolerance;

    int N = 10000000;


    for(int i=0; i<N; i++){

        // Start the integration
        h = dt;
        while (true){
            fehl_step(y, h, yp, k2, k3, k4, k5, k6, s);
            // Compute and test allowable tolerance.
            // Get average of the function magnitudes at t and t + h
            error_et = std::abs(y[0]) + std::abs(k2[0]) + error_absolute;
            // We assume that the errors set are positive!

            // Truncation error
            error_ete = std::abs(-2090.0*yp[0] + 21970.0 * k4[0] - 15048.0 * k5[0] + 22528.0 * k3[0] - 27360.0 * k6[0]);

            // Get the relative error
            error_eot = error_ete/error_et;

            // Do the same at second dimension
            error_et = std::abs(y[1]) + std::abs(k2[1]) + error_absolute;
            // We assume that the errors set are positive!

            // Truncation error
            error_ete = std::abs(-2090.0*yp[1] + 21970.0 * k4[1] - 15048.0 * k5[1] + 22528.0 * k3[1] - 27360.0 * k6[1]);

            // Take the max relative error of the two
            error_eot = std::fmax(error_eot, error_ete/error_et);

            // Final error tolerance evaluation.
            error_tolerance = std::abs(h) * error_eot * error_scale / 752400.0;

            if (error_tolerance <= 1.0){
                // Successful step
                break;
            }
            // Unsuccessful step
            if (error_tolerance < 59049.0){
                h_scale = 0.9 / std::pow(error_tolerance, 0.2);
            }
            else {
                h_scale = 0.1;
            }
            printf("scalling\n");
            h = h_scale * h;
        } // END SUBSTEPPING. Successful integration
        // Here you would increase the parametric time. But it doesn't matter
        // since we will just follow the value of the contour in the current
        // position.
        y[0] = s[0];
        y[1] = s[1];
        pde(y, yp); // Renew the derivative values

        // Increase step if possible
        if ( 0.0001889568 < error_tolerance )
        {
          h_scale = 0.9 / std::pow ( error_tolerance, 0.2 );
        }
        else
        {
          h_scale = 5.0;
        }

        contour_rel_error = std::abs(function(y) - orig_contour_val) / orig_contour_val;
        printf("%.8f %.8f\n", y[0], y[1]);
        max_error = std::fmax(max_error, contour_rel_error);
        if (contour_rel_error > relerr) {
            printf("%.4e is larger than the prescribed %.4e. After %d steps!\n", contour_rel_error, relerr, i);
            break;
        }
    }
    printf("Finished successfully with maximum relative error: %e.\n", max_error);
    return 0;
}
