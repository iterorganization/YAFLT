/*
Example showing how to use the RKF45 solver. Also it shows the stability by
following a single magnetic surface for quite some time.

It prints the Phi, R, Z values as well as the relative error of the flux value,
where the initial flux value is compared as reference.

One can plot the contents simply with

.. code-block:: python

   phi = []
   r = []
   z = []
   with open("OUT", "r") as f:
      for line in f:
          sline = line.split()
          phi.append(float[sline[0]])
          r.append(float[sline[1]])
          z.append(float[sline[2]])

   import matplotlib.pyplot as plt
   plt.plot(r, z)
   plt.show()
*/

#include <cmath>
#include <eqdsk_drsep10_data.hpp>
#include <bicubic.hpp>
#include <cstdio>
#include <chrono>

class TEST_SOLVER
{
    public:
        TEST_SOLVER();
        ~TEST_SOLVER();
        BICUBIC_INTERP *m_interp_psi;
        BI_DATA *m_context;
        double m_vacuum_fpol;
        void flt_pde(double y[2], double yp[2]);
        void fehl_step(double y[2], double h, double yp[2], double k1[2],
               double k2[2], double k3[2], double k4[2], double k5[2],
               double k6[2], double s[2]);
};

TEST_SOLVER::TEST_SOLVER(){
    m_interp_psi = new BICUBIC_INTERP();
    m_context = new BI_DATA();
}

TEST_SOLVER::~TEST_SOLVER(){
    delete m_interp_psi;
    delete m_context;
}

void TEST_SOLVER::flt_pde(double y[2], double yp[2]){
    double factor;

    m_context->r = y[0];
    m_context->z = y[1];
    m_interp_psi->getValues(m_context);
    factor = y[0] / m_vacuum_fpol;

    yp[0] = -m_context->valdy * factor;
    yp[1] =  m_context->valdx * factor;
}

void TEST_SOLVER::fehl_step(double y[2], double h, double yp[2], double k1[2],
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
    // flt_pde(t + ch, k6, k2);
    flt_pde(k6, k2);

    // k3 = y + h * (3/32 yp + 9/32 k2)
    ch = 3.0 * h / 32.0;
    k6[0] = y[0] + ch * (yp[0] + 3.0 * k2[0]);
    k6[1] = y[1] + ch * (yp[1] + 3.0 * k2[1]);

    // Evaluate the function for k3
    flt_pde(k6, k3);

    // k4 = y + h * (1932/2197 yp - 7200 / 2197 k2 + 7296/2197 k3)
    ch = h / 2197.0;
    k6[0] = y[0] + ch * (1932.0 * yp[0] + (7296.0 * k3[0] - 7200.0 * k2[0]));
    k6[1] = y[1] + ch * (1932.0 * yp[1] + (7296.0 * k3[1] - 7200.0 * k2[1]));

    // Evaluate the function for k4
    // flt_pde(t + 12.0 * h / 13.0, k6, k4);
    flt_pde(k6, k4);

    // k5 = y + h * (439/216 yp - 8 k2 + 3680/513 k3 - 845/4104 k4)
    ch = h/ 4104.0;
    k6[0] = y[0] + ch * ((8341.0 * yp[0] - 845.0 * k4[0]) + (29440.0 * k3[0] - 32832.0 * k2[0]));
    k6[1] = y[1] + ch * ((8341.0 * yp[1] - 845.0 * k4[1]) + (29440.0 * k3[1] - 32832.0 * k2[1]));

    // Evaluate the function at the new position
    // flt_pde(t + h, k6, k5);
    flt_pde(k6, k5);

    // Finally k6 = y + h * (-8/27 yp + 2 k2 - 3544/2565 k3 + 1850/4104 k4 - 11/40 k5)
    // The k2 factor is not used in the last evaluations, so we will use k2
    // here for storing result for k6.
    ch = h / 20520.0;
    k2[0] = y[0] + ch * ((-6080.0 * yp[0] + (9295.0 * k4[0] - 5643.0 * k5[0])) + (41040.0 * k2[0] - 28352.0 * k3[0]));
    k2[1] = y[1] + ch * ((-6080.0 * yp[1] + (9295.0 * k4[1] - 5643.0 * k5[1])) + (41040.0 * k2[1] - 28352.0 * k3[1]));

    // Evaluate the function for k6
    // flt_pde(t + 0.5 * h, k2, k6);
    flt_pde(k2, k6);

    // Approximate solution 4th order.
    // ch = h / 20520;
    // s[0] = y[0] + ch * (yp[0] * 2375.0 + k3[0] * 11264.0 + (k4[0] * 10985.0 - k5[0] * 4104.0));
    // s[1] = y[1] + ch * (yp[1] * 2375.0 + k3[1] * 11264.0 + (k4[1] * 10985.0 - k5[1] * 4104.0));

    // // Approximate solution 5th order.
    ch = h / 282150.0;
    s[0] = y[0] + ch * (33440.0 * yp[0] + 142805.0 * k4[0] - 50787.0 * k5[0] + 146432.0 * k3[0] + 10260.0 * k6[0]);
    s[1] = y[1] + ch * (33440.0 * yp[1] + 142805.0 * k4[1] - 50787.0 * k5[1] + 146432.0 * k3[1] + 10260.0 * k6[1]);
}

int main(){

    TEST_SOLVER *solver = new TEST_SOLVER();

    int n_cols = drsep10_Rs.size();
    int n_rows = drsep10_Zs.size();

    std::vector<std::vector<double>> reshaped_Psi;
    std::vector<double> buffer;
    for (int i=0; i<n_rows; i++){
        buffer.clear();
        for (int j=0; j<n_cols; j++){
            buffer.push_back(drsep10_Psi[i * n_cols + j]);
        }
        // reshaped_Psi.insert(reshaped_Psi.begin(), buffer);
        reshaped_Psi.push_back(buffer);
    }
    solver->m_interp_psi->setArrays(drsep10_Rs, drsep10_Zs, reshaped_Psi);
    double drsep10_fpolvacuum = -32.86;
    solver->m_vacuum_fpol = drsep10_fpolvacuum;

    // Set the initial parameters. And follow a fieldline for absurd amount
    // of distance

    // Containers for the position and derivative
    double y[2], yp[2];
    double time=0, time_step=.01;
    double h, h_scale;
    bool failed_step;
    // Flag=1, 2 sets the solver in a single step-mode. Use 1 for first
    // initialization. Required!

    // Desired accuracy.
    double relerr=1e-4, abserr=1e-4;
    double desired_relerr=1e-4;

    BI_DATA *context = new BI_DATA();
    solver->m_interp_psi->populateContext(solver->m_context);
    solver->m_interp_psi->populateContext(context);
    y[0] = 6.4;
    y[1] = 0.0;
    solver->flt_pde(y, yp);
    double rel_error, flux;
    // Get the initial flux value to compare on the new locations.
    context->r = y[0];
    context->z = y[1];
    solver->m_interp_psi->getValues(context);
    flux = context->val;
    unsigned int counter=0;
    unsigned int N=100000000;
    unsigned int ten_percentile = N/10;
    double percentage;

    bool ok = true;
    auto begin = std::chrono::high_resolution_clock::now();
    double max_error=0.0;
    double k1[2], k2[2], k3[2], k4[2], k5[2], k6[2];
    // Solution array
    double s[2];
    // Error values
    double last_error;
    double error_et;
    double error_ete;
    double error_lte_relative;
    h = time_step;
    while (counter < N){
        failed_step=false;
        while (1){
            solver->fehl_step(y, h, yp, k1, k2, k3, k4, k5, k6, s);
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
            if (error_lte_relative < 1.0){
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
        }

        if (h < time_step){
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
            h = std::fmin(h_scale * h, time_step);
        }

        if (counter % ten_percentile == 0){
            percentage = 10*static_cast<double>(counter / ten_percentile);
            printf("%.2f\n", percentage);
        }
        y[0] = s[0];
        y[1] = s[1];
        solver->flt_pde(y, yp);
        time = time + h;

        context->r = y[0];
        context->z = y[1];
        solver->m_interp_psi->getValues(context);
        rel_error = std::fabs(std::fabs(context->val - flux) / flux);
        max_error=std::fmax(rel_error, max_error);
        last_error = rel_error;
        if (rel_error > desired_relerr){
            printf("%d %f %f %f %e\n", counter, time / (2 * 3.14), y[0], y[1], rel_error);
            ok = false;
            break;
        }
        // Set the flag back to 2 and ignore other messages from the solver.
        counter = counter + 1;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    float performance = counter *1e3 / elapsed.count();
    printf("Performance of steps per second: %f M/s\n", performance);
    printf("Max error: %e\n", max_error);
    printf("Last error: %e\n", last_error);
    printf("Final time: %f\n", time);
    // Free objects
    delete solver; // Also deletes the interpolation object tied to it.
    delete context;
    if (!ok){
        return 1;
    }
    return 0;
}