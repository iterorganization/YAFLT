/*
  Software for testing RKF4(5) accuracy by tracing a line on a contour in a
  cylindrical coordinate system.

  For input we take a simple analytical function:

  f(x, y) = 0.5*(x - 3)**2 + 2.3*(y - 5)**2

  The gradient:

  df/dx = 0.5*2(x-5)
  df/dy = 2.3*2(y-5)

  We do not directly solve this. What we do is we actually search for s(t) of
  which the gradient is:

  ds/dx = -df/dy
  ds/dy =  df/dx

  And then we integrate over the parametric time t.

Table III from NASA TR R-315


+---+-------+----------------------------------------------------------+-----------+-------------+------------+
| L | Ak    |                        Bkl                               | Ck        | Ch          | Ce         |
+===+=======+===========+============+============+===========+========+===========+=============+============+
| k | /     | 0         | 1          | 2          | 3         | 4      |           |             | 1          |
+---+-------+-----------+------------+------------+-----------+--------+-----------+-------------+------------+
| 0 | 0     | 0         | /          | /          | /         | /      | 25/216    | 16/135      | -1/360.0   |
+---+-------+-----------+------------+------------+-----------+--------+-----------+-------------+------------+
| 1 | 1/4   | 1/4       | /          | /          | /         | /      | 0         | 0           | 1          |
+---+-------+-----------+------------+------------+-----------+--------+-----------+-------------+------------+
| 2 | 3/8   | 3/32      | 9/32       | /          | /         | /      | 1408/2565 | 6656/12825  | 128/4275.0 |
+---+-------+-----------+------------+------------+-----------+--------+-----------+-------------+------------+
| 3 | 12/13 | 1932/2197 | -7200/2197 | 7296/2197  | /         | /      | 2197/4104 | 28561/56430 | 2197/75240 |
+---+-------+-----------+------------+------------+-----------+--------+-----------+-------------+------------+
| 4 | 1     | 439/216   | -8         | 3680/513   | -845/4104 | /      | -1/5      | -9/50       | -1/50      |
+---+-------+-----------+------------+------------+-----------+--------+-----------+-------------+------------+
| 5 | 1/2   | -8/27     | 2          | -3544/2565 | 1859/4104 | -11/40 | /         | 2/55        | -2/55      |
+---+-------+-----------+------------+------------+-----------+--------+-----------+-------------+------------+
                                                              |   LCM  | 20520     |  282150     |  376200    |
                                                              +--------+-----------+-------------+------------+


f[i] = h Fun(x + Ak * h, sum(l=0, i-1) Bil fl )
4th order
y = y0 + h * sum(i=0, 4) Ck[i] f[i]
5th order
y = y0 + h * sum(i=0, 5) Ch[i] f[i]


Local truncation error
lte = |y_5th - y_4th|

 */

// PDE derivative function

#include <cmath>
#include <cstdio>
#include <chrono>

inline double function(double y[2]){
    return 0.5*(y[0] - 3.0) * (y[0] - 3.0) + 2.3*(y[1] - 5.0) * (y[1] - 5.0);
}

inline void pde(double y[2], double yp[2]){

    yp[0] = -4.6 * (y[1] - 5.0);
    yp[1] =  (y[0] - 3.0);
}

// Fehlberg step

void fehl_step(double y[2], double h, double yp[2], double k1[2],
               double k2[2], double k3[2], double k4[2], double k5[2],
               double k6[2], double s[2]){
    // Calculate the k1, k2, k3, k4, k5, k6 factors of the RKF4(5) method and
    // the approximate solution. The k* factors are used for error estimation
    // while the s holds the approximate solution. If the error is acceptable
    // then the approximate solution is the new solution.

    // The factors are calculated based on Table II in Fehlberg, NASA Technical
    // Report 315.

    // We use the re-parametrized mode, meaning we actually integrate over
    // parametric time t.

    // Use k6 as container for  calculating other k factors (used as argument
    // for the pde method).
    // k1 is h pde(t, y)
    double ch;
    k1[0] = yp[0];
    k1[1] = yp[1];

    // k2 = y + h * yp / 4
    ch = 0.25 * h;
    k6[0] = y[0] + ch * yp[0];
    k6[1] = y[1] + ch * yp[1];

    // Evaluate the function for k2
    // pde(t + ch, k6, k2);
    pde(k6, k2);

    // k3 = y + h * (3/32 yp + 9/32 k2)
    ch = 3.0 * h / 32.0;
    k6[0] = y[0] + ch * (yp[0] + 3.0 * k2[0]);
    k6[1] = y[1] + ch * (yp[1] + 3.0 * k2[1]);

    // Evaluate the function for k3
    pde(k6, k3);

    // k4 = y + h * (1932/2197 yp - 7200 / 2197 k2 + 7296/2197 k3)
    ch = h / 2197.0;
    k6[0] = y[0] + ch * (1932.0 * yp[0] + (7296.0 * k3[0] - 7200.0 * k2[0]));
    k6[1] = y[1] + ch * (1932.0 * yp[1] + (7296.0 * k3[1] - 7200.0 * k2[1]));

    // Evaluate the function for k4
    // pde(t + 12.0 * h / 13.0, k6, k4);
    pde(k6, k4);

    // k5 = y + h * (439/216 yp - 8 k2 + 3680/513 k3 - 845/4104 k4)
    ch = h/ 4104.0;
    k6[0] = y[0] + ch * (8341.0 * yp[0] - 32832.0 * k2[0] + 29440.0 * k3[0] - 845.0 * k4[0]);
    k6[1] = y[1] + ch * (8341.0 * yp[1] - 32832.0 * k2[1] + 29440.0 * k3[1] - 845.0 * k4[1]);

    // Evaluate the function at the new position
    // pde(t + h, k6, k5);
    pde(k6, k5);

    // Finally k6 = y + h * (-8/27 yp + 2 k2 - 3544/2565 k3 + 1850/4104 k4 - 11/40 k5)
    // The k2 factor is not used in the last evaluations, so we will use k2
    // here for storing result for k6.
    ch = h / 20520;
    k2[0] = y[0] + ch * (-6080.0 * yp[0] + 41040 * k2[0] - 28352.0 * k3[0] + 9295.0 * k4[0] - 5643.0 * k5[0]);
    k2[1] = y[1] + ch * (-6080.0 * yp[1] + 41040 * k2[1] - 28352.0 * k3[1] + 9295.0 * k4[1] - 5643.0 * k5[1]);

    // Evaluate the function for k6
    // pde(t + 0.5 * h, k2, k6);
    pde(k2, k6);

    // Approximate solution 4th order.
    // ch = h / 20520;
    // s[0] = y[0] + ch * (k1[0] * 2375.0 + k3[0] * 11264.0 + (k4[0] * 10985.0 - k5[0] * 4104.0));
    // s[1] = y[1] + ch * (k1[1] * 2375.0 + k3[1] * 11264.0 + (k4[1] * 10985.0 - k5[1] * 4104.0));

    // // Approximate solution 5th order.
    ch = h / 282150.0;
    s[0] = y[0] + ch * (k1[0] * 33440.0 + k3[0] * 146432.0 + k4[0] * 142805.0 - k5[0] * 50787.0 + k6[0] * 10260.0);
    s[1] = y[1] + ch * (k1[1] * 33440.0 + k3[1] * 146432.0 + k4[1] * 142805.0 - k5[1] * 50787.0 + k6[1] * 10260.0);
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
    bool failed_step;

    // Step increment arrays
    double k1[2], k2[2], k3[2], k4[2], k5[2], k6[2];
    // Solution array
    double s[2];
    // Error values
    double error_et;
    double error_ete;
    double error_lte_relative;

    int N = 1000000;
    int count;
    bool ok = true;
    printf("Running with desired relerr %e and step %f.\n", relerr, dt);
    h = dt;
    auto begin = std::chrono::high_resolution_clock::now();
    for(int i=0; i<N; i++){

        // Start the integration
        failed_step = false;
        while (1){
            fehl_step(y, h, yp, k1, k2, k3, k4, k5, k6, s);
            // Compute and test allowable tolerance.
            // Get average of the function magnitudes at t and t + h
            error_et = std::abs(s[0]) * relerr + abserr;
            // We assume that the errors set are positive!

            // Truncation error
            error_ete = std::abs((-1.0 * k1[0] / 360.0 + 128.0 * k3[0] / 4275.0 + 2197.0 * k4[0] / 75240.0 - 1.0 * k5[0] / 50.0 - 2.0 * k6[0] / 55.0) * h);

            // Get the relative error
            error_lte_relative = error_ete/error_et;

            // Do the same at second dimension
            error_et = std::abs(s[1]) * relerr + abserr;
            // We assume that the errors set are positive!
            error_ete = std::abs((-1.0 * k1[1] / 360.0 + 128.0 * k3[1] / 4275.0 + 2197.0 * k4[1] / 75240.0 - 1.0 * k5[1] / 50.0 - 2.0 * k6[1] / 55.0) * h);

            // Truncation max error
            error_lte_relative = std::fmax(error_lte_relative, error_ete/error_et);

            // Final error tolerance evaluation.
            if (error_lte_relative <= 1.0){
                // Successful step
                break;
            }
            failed_step=true;
            // Unsuccessful step

            // Lower it by maximum one order of magnitude
            if (error_lte_relative < 59049.0){
                h_scale = 0.9 / std::pow(error_lte_relative, 0.2);
            }
            else {
                h_scale = 0.1; // Limit lowering by one order of magnitudes
            }
            h = h_scale * h;
        } // END SUBSTEPPING. Successful integration
        // Here you would increase the parametric time. But it doesn't matter
        // since we will just follow the value of the contour in the current
        // position.
        y[0] = s[0];
        y[1] = s[1];
        pde(y, yp); // Renew the derivative values

        // Increase step if possible. Have a maximum increase of 5 times. Limit
        // it by the user setting.
        if (h < dt){
            if ( 0.0001889568 < error_lte_relative )
            {
              h_scale = 0.9 / std::pow ( error_lte_relative, 0.2 );
            }
            else
            {
              h_scale = 5.0;
            }
            if (failed_step){
                h_scale = std::fmin(1.0, h_scale);
            }
            h = std::fmin(h_scale * h, dt);
        }
        // printf("%f\n", h_scale);
        contour_rel_error = std::abs(function(y) - orig_contour_val) / orig_contour_val;
        max_error = std::fmax(max_error, contour_rel_error);
        count = i;
        // if (contour_rel_error > relerr) {
        //     printf("%.4e is larger than the prescribed %.4e. After %d steps!\n", contour_rel_error, relerr, i);
        //     ok = false;
        //     break;
        // }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    float performance = count * 1e3 / elapsed.count();
    printf("Performance of steps per second: %f M/s\n", performance);
    if (ok) {
        printf("Finished after %d steps with maximum relative error: %e.\n", N, max_error);
        return 0;
    }
    else {
        return 1;
    }
}
