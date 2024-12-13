/// Testing the accuracy of the naive definition of the code. Selected a
/// function from which we can easily calculate the value of the function and
/// it's derivative

#include <bicubic.hpp>

#include <vector>
#include <stdio.h>
#include <stdlib.h> /*srand, rand*/
#include <time.h> /*For srand seed*/
#include <cmath>

double fun(double x, double y, double a, double b){
    return pow(x - a, 2) * pow(y - b, 2);
}

int main(){
    srand(time(NULL));
    BICUBIC_INTERP *obj = new BICUBIC_INTERP();
    BI_DATA *context = new BI_DATA();
    int Nx = 50;
    int Ny = 100;
    double Dx = 5.0/Nx;
    double Dy = 10.0/Ny;

    /// The domain is defined over [0, 5]x[0, 10], with the following function
    /// F(x, y) = exp(-(x - alpha)**2) exp(-(y - 5)**2)

    double alpha, beta;
    alpha=2.5;
    beta=2.5;

    std::vector<double> X, Y;
    X.resize(Nx);
    Y.resize(Ny);

    for (int i=0; i<Nx; i++){
        X[i] = i * Dx;
    }
    for (int j=0; j<Ny; j++){
        Y[j] = j * Dy;
    }

    std::vector<std::vector<double>> F;
    std::vector<double> buffer;

    for (int j=0; j<Ny; j++){ // Rows
        buffer.clear();
        for (int i=0; i<Nx; i++){ // Columns
            buffer.push_back(fun(i*Dx, j*Dy, alpha, beta));
        }
        F.push_back(buffer);
    }

    obj->setArrays(X, Y, F);
    obj->populateContext(context);

    // For Random prepare lower and upper boundaries
    double xmin, xmax, xdiff, ymin, ymax, ydiff;
    xmin = X[0];
    xmax = X[X.size() - 2];
    xdiff = xmax - xmin;
    ymin = Y[0];
    ymax = Y[Y.size() - 2];
    ydiff = ymax - ymin;

    //Random check
    double buff_x, buff_y, buff;
    double exact_f, exact_fdx, exact_fdy;
    int N = 10*1000*1000;

    bool print;
    double abs_max_f=-1, abs_max_fdx=-1, abs_max_fdy=-1;
    int counts=0;

    for (int i=0; i<N; i++){
        //
        print=false;
        buff_x = xmin + xdiff * rand() / RAND_MAX;
        buff_y = ymin + ydiff * rand() / RAND_MAX;

        // exact_f = exp(-pow(buff_x - alpha, 2)) * exp(-pow(buff_y - beta, 2));
        exact_f = fun(buff_x, buff_y, alpha, beta);

        context->r = buff_x;
        context->z = buff_y;
        obj->getValues(context);

        // Compare values
        buff = std::fabs(context->val - exact_f) / exact_f;
        if (buff > 0.0001) {
            print=true;
        }
        if (buff > abs_max_f) {
            abs_max_f = buff;
        }

        // Compare x derivative
        exact_fdx = 2.0 * (buff_x - alpha) * pow(buff_y - beta, 2.0);
        // Compare values
        buff = std::fabs(context->valdx - exact_fdx) / exact_fdx;
        if (buff > 0.0001) {
            print=true;
        }
        if (buff > abs_max_fdx) {
            abs_max_fdx = buff;
        }

        // Compare y derivative
        exact_fdy = 2.0 * (buff_y - beta) * pow(buff_x - alpha, 2.0);
        // Compare values
        buff = std::fabs(context->valdy - exact_fdy) / exact_fdy;
        if (buff > 0.0001) {
            print=true;
        }
        if (buff > abs_max_fdy) {
            abs_max_fdy = buff;
        }
        if (print) {
            counts = counts + 1;
            // printf("Point %f %f\n", buff_x, buff_y);
            // printf("Type:\tVal\tValdx\tValdy\n");
            // printf("Naive\t%f\t%f\t%f\n", context->val, context->valdx, context->valdy);
            // printf("Exact\t%f\t%f\t%f\n", exact_f, exact_fdx, exact_fdy);
            // printf("Error\t%f\t%f\t%f\n", abs_max_f, abs_max_fdx, abs_max_fdy);
        }
    }
    double ratio=(double)counts/N;
    // printf("Total tries: %d\n", N);
    // printf("Number of high error points: %d (%.2f)\n", counts, ratio*100);
    // printf("Maximum relative error f: %f\n", abs_max_f);
    // printf("Maximum relative error fdx: %f\n", abs_max_fdx);
    // printf("Maximum relative error fdy: %f\n", abs_max_fdy);
    // printf("Location %f %f\n", abs_max_r, abs_max_z);
    // printf("R goes from %f to %f\n", xmin, xmax);
    // printf("Z goes from %f to %f\n", ymin, ymax);

    if (ratio > 0.05){
        // If the error ratio is more than 5% return 1
        return 1;
    }
    return 0;
}