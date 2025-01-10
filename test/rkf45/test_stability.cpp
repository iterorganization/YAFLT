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
#include <rkf45.hpp>
#include <bicubic.hpp>
#include <cstdio>
int main(){

    RKF45 *solver = new RKF45();
    BICUBIC_INTERP *naive_interp = new BICUBIC_INTERP();

    solver->set_interpolator(naive_interp);

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
    naive_interp->setArrays(drsep10_Rs, drsep10_Zs, reshaped_Psi);
    double drsep10_fpolvacuum = -32.86;
    solver->set_vacuum_fpol(drsep10_fpolvacuum);

    // Set the initial parameters. And follow a fieldline for absurd amount
    // of distance

    // Containers for the position and derivative
    double y[2], yp[2];
    double time=0, new_time=0, time_step=.01;
    // Flag=1, 2 sets the solver in a single step-mode. Use 1 for first 
    // initialization. Required!
    int flag=1;

    // Desired accuracy.
    double relerr=1e-4, abserr=1e-4;
    double desired_relerr=1e-4;
    y[0] = 6.4;
    y[1] = 0.0;

    BI_DATA *context = new BI_DATA();
    naive_interp->populateContext(context);
    double rel_error, flux;
    // Get the initial flux value to compare on the new locations.
    context->r = y[0];
    context->z = y[1];
    naive_interp->getValues(context);
    flux = context->val;
    unsigned int counter=0;
    unsigned int N=1000000;
    unsigned int ten_percentile = N/10;
    double percentage;

    bool ok = true;

    while (counter < N){
        if (counter % ten_percentile == 0){
            percentage = 10*static_cast<double>(counter / ten_percentile);
            printf("%.2f\n", percentage);
        }
        new_time = time + time_step;
        flag = solver->r8_rkf45(y, yp, &time, new_time, &relerr, abserr, flag);
        // printf("%f %f %d\n", y[0], y[1], flag);
        context->r = y[0];
        context->z = y[1];
        naive_interp->getValues(context);
        rel_error = std::fabs(std::fabs(context->val - flux) / flux);
        if (rel_error > desired_relerr){
            printf("%f %f %f %e\n", time / (2 * 3.14), y[0], y[1], rel_error);
            ok = false;
            break;
        }
        // Set the flag back to 2 and ignore other messages from the solver.
        flag = 2;
        counter = counter + 1;
    }

    // Free objects
    delete solver; // Also deletes the interpolation object tied to it.
    delete context;
    if (!ok){
        return 1;
    }
    return 0;
}