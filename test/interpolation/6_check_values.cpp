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

bool cmpd(double a, double b, double eps=1e-3){
    return (fabs(a - b) < eps);
}

void process_input_arrays(std::vector<double> x, std::vector<double> y,
                         std::vector<std::vector<double>> f, double *out_x,
                        double *out_y, double *out_f, double *out_fdx,
                        double *out_fdy, double *out_fdxdy){
    // Prepare the arrays for interpolation

    int n_cols, n_rows, n_xy;
    n_cols = x.size(); // Number of columns
    n_rows = y.size(); // Number of rows
    n_xy = n_rows * n_cols;



    int i, j, k;
    // Simply copy data
    for (i=0; i<n_cols; i++){
        out_x[i] = x[i];
    }

    for (j=0; j<n_rows; j++){
        out_y[j] = y[j];
    }
    // n_rows is the number of rows. n_cols is the number of columns.
    for (i=0; i<n_rows;i++){
        for(j=0; j<n_cols;j++){
            k = i * n_cols + j;
            out_f[k] = f[i][j];
        }
    }

    // Now the derivatives.
    // Upper left border values
    out_fdx[0] = 0.5 * (f[0][1] - f[0][0]);
    out_fdy[0] = 0.5 * (f[1][0] - f[0][0]);
    out_fdxdy[0] = 0.25 * (f[1][1] - f[0][1] - f[1][0] + f[0][0]);

    // Upper right border values
    k = n_cols - 1;
    out_fdy[k] = 0.5 * (f[1][n_cols - 1] - f[0][n_cols - 1]);
    out_fdx[k] = 0.5 * (f[0][n_cols - 1] - f[0][n_cols - 2]);
    out_fdxdy[k] = 0.25 * (f[1][n_cols - 1] - f[0][n_cols - 1] - f[1][n_cols - 2] + f[0][n_cols - 2]);

    // Lower left border values
    k = (n_rows - 1) * n_cols;
    out_fdx[k] = 0.5 * (f[n_rows - 1][1] - f[n_rows - 1][0]);
    out_fdy[k] = 0.5 * (f[n_rows - 1][0] - f[n_rows - 2][0]);
    out_fdxdy[k] = 0.25 * (f[n_rows - 1][1] - f[n_rows - 2][1] - f[n_rows - 1][0] + f[n_rows - 2][0]);

    // Lower right border values
    k = n_xy - 1;
    out_fdx[k] = 0.5 * (f[n_rows - 1][n_cols - 1] - f[n_rows - 1][n_cols - 2]);
    out_fdy[k] = 0.5 * (f[n_rows - 1][n_cols - 1] - f[n_rows - 2][n_cols - 1]);
    out_fdxdy[k] = 0.25 * (f[n_rows - 1][n_cols - 1] - f[n_rows - 2][n_cols - 1] - f[n_rows - 1][n_cols - 2] + f[n_rows - 2][n_cols - 2]);


    // Left and right border lines
    for(int i=1; i<n_rows-1; i++){
        k = i * n_cols;
        // Forward derivative - left side
        out_fdx[k]= 0.5 * (f[i][1] - f[i][0]);
        out_fdy[k]= 0.5 * (f[i+1][0] - f[i-1][0]);
        out_fdxdy[k]= 0.25 * (f[i+1][1] - f[i-1][1] - f[i+1][0] + f[i-1][0]);

        // Backward derivative - right side
        k = (i+1) * n_cols - 1;
        out_fdx[k]  = 0.5 * (f[i][n_cols - 1] - f[i][n_cols - 2]);
        out_fdy[k]  = 0.5 * (f[i+1][n_cols - 1] - f[i-1][n_cols - 1]);
        out_fdxdy[k]= 0.25 * (f[i+1][n_cols - 1] - f[i-1][n_cols - 1] - f[i+1][n_cols - 2] + f[i-1][n_cols - 2]);
    }

    // Upper and bottom border lines
    for(int j=1; j<n_cols-1;j++){
        // Upper side
        k = j;
        out_fdx[k]= 0.5 * (f[0][j+1] - f[0][j-1]);
        out_fdy[k]= 0.5 * (f[1][j] - f[0][j]);
        out_fdxdy[k]= 0.25 * (f[1][j+1] - f[0][j+1] - f[1][j-1] + f[0][j-1]);

        // Lower side
        k = (n_rows - 1) * n_cols + j;
        out_fdx[k] = 0.5 * (f[n_rows-1][j+1] - f[n_rows-1][j-1]);
        out_fdy[k] = 0.5 * (f[n_rows-1][j] - f[n_rows-2][j]);
        out_fdxdy[k] = 0.25 * (f[n_rows-1][j+1] - f[n_rows-2][j+1] - f[n_rows-1][j-1] + f[n_rows-2][j-1]);
    }

    // Centered finite differences.
    // Now for the rest of the space
    for(int i=1; i<n_rows-1;i++){
        for(int j=1; j<n_cols-1;j++){
            // fdx[i][j] = 0.5 *(f[i][j+1] - f[i][j-1]);
            k = i * n_cols + j;
            out_fdx[k] = 0.5 *(f[i][j+1] - f[i][j-1]);
            out_fdy[k] = 0.5 *(f[i+1][j] - f[i-1][j]);
            out_fdxdy[k] = 0.25 * (f[i+1][j+1] - f[i-1][j+1] - f[i+1][j-1] + f[i-1][j-1]);
        }
    }

    // Higher order finite difference
    for(int i=2; i<n_rows-2;i++){
        for(int j=2; j<n_cols-2;j++){
            k = i * n_cols + j;
            out_fdx[k] = (-f[i][j+2] + 8.0 * f[i][j+1] - 8.0 * f[i][j-1] + f[i][j-2]) / 12.0;
            out_fdy[k] = (-f[i+2][j] + 8.0 * f[i+1][j] - 8.0 * f[i-1][j] + f[i-2][j]) / 12.0;
            out_fdxdy[k] = (f[i+2][j+2] - 8.0 * f[i+1][j+2] + 8.0 * f[i-1][j+2] - f[i-2][j+2] + \
                             8.0 * (-f[i+2][j+1] + 8.0 * f[i+1][j+1] - 8.0 * f[i-1][j+1] + f[i-2][j+1]) - \
                             8.0 * (-f[i+2][j-1] + 8.0 * f[i+1][j-1] - 8.0 * f[i-1][j-1] + f[i-2][j-1]) - \
                             f[i+2][j-2] + 8.0 * f[i+1][j-2] - 8.0 * f[i-1][j-2] + f[i-2][j-2]) / 144.0;
        }
    }

    // 6th order centered difference approximations
    for(int i=3; i<n_rows-3;i++){
        for(int j=3; j<n_cols-3;j++){
            k = i * n_cols + j;
            out_fdx[k] = (-f[i][j-3] + 9.0 * f[i][j-2] - 45.0 * f[i][j-1] + 45.0 * f[i][j+1] - 9.0 * f[i][j+2] + f[i][j+3]) / 60.0;
            out_fdy[k] = (-f[i-3][j] + 9.0 * f[i-2][j] - 45.0 * f[i-1][j] + 45.0 * f[i+1][j] - 9.0 * f[i+2][j] + f[i+3][j]) / 60.0;

            out_fdxdy[k] =    (-1.0 * (-f[i-3][j-3] + 9.0 * f[i-2][j-3] - 45.0 * f[i-1][j-3] + 45.0 * f[i+1][j-3] - 9.0 * f[i+2][j-3] + f[i+3][j-3]) + \
                                9.0 * (-f[i-3][j-2] + 9.0 * f[i-2][j-2] - 45.0 * f[i-1][j-2] + 45.0 * f[i+1][j-2] - 9.0 * f[i+2][j-2] + f[i+3][j-2]) - \
                               45.0 * (-f[i-3][j-1] + 9.0 * f[i-2][j-1] - 45.0 * f[i-1][j-1] + 45.0 * f[i+1][j-1] - 9.0 * f[i+2][j-1] + f[i+3][j-1]) + \
                               45.0 * (-f[i-3][j+1] + 9.0 * f[i-2][j+1] - 45.0 * f[i-1][j+1] + 45.0 * f[i+1][j+1] - 9.0 * f[i+2][j+1] + f[i+3][j+1]) - \
                                9.0 * (-f[i-3][j+2] + 9.0 * f[i-2][j+2] - 45.0 * f[i-1][j+2] + 45.0 * f[i+1][j+2] - 9.0 * f[i+2][j+2] + f[i+3][j+2]) + \
                                1.0 * (-f[i-3][j+3] + 9.0 * f[i-2][j+3] - 45.0 * f[i-1][j+3] + 45.0 * f[i+1][j+3] - 9.0 * f[i+2][j+3] + f[i+3][j+3])) / 3600.0;
        }
    }

    // 8th order centered difference approximation
    for(int i=4; i<n_rows-4;i++){
        for(int j=4; j<n_cols-4;j++){
            k = i * n_cols + j;
            out_fdx[k] = (3.0 * f[i][j-4] - 32.0 * f[i][j-3] + 168.0 * f[i][j-2] - 672.0 * f[i][j-1] + 672.0 * f[i][j+1] - 168.0 * f[i][j+2] + 32.0 * f[i][j+3] - 3.0 * f[i][j+4]) / 840.0;
            out_fdy[k] = (3.0 * f[i-4][j] - 32.0 * f[i-3][j] + 168.0 * f[i-2][j] - 672.0 * f[i-1][j] + 672.0 * f[i+1][j] - 168.0 * f[i+2][j] + 32.0 * f[i+3][j] - 3.0 * f[i+4][j]) / 840.0;

            out_fdxdy[k] =    (3.0 * (3.0 * f[i-4][j-4] - 32.0 * f[i-3][j-4] + 168.0 * f[i-2][j-4] - 672.0 * f[i-1][j-4] + 672.0 * f[i+1][j-4] - 168.0 * f[i+2][j-4] + 32.0 * f[i+3][j-4] - 3.0 * f[i+4][j-4]) - \
                              32.0 * (3.0 * f[i-4][j-3] - 32.0 * f[i-3][j-3] + 168.0 * f[i-2][j-3] - 672.0 * f[i-1][j-3] + 672.0 * f[i+1][j-3] - 168.0 * f[i+2][j-3] + 32.0 * f[i+3][j-3] - 3.0 * f[i+4][j-3]) + \
                             168.0 * (3.0 * f[i-4][j-2] - 32.0 * f[i-3][j-2] + 168.0 * f[i-2][j-2] - 672.0 * f[i-1][j-2] + 672.0 * f[i+1][j-2] - 168.0 * f[i+2][j-2] + 32.0 * f[i+3][j-2] - 3.0 * f[i+4][j-2]) - \
                             672.0 * (3.0 * f[i-4][j-1] - 32.0 * f[i-3][j-1] + 168.0 * f[i-2][j-1] - 672.0 * f[i-1][j-1] + 672.0 * f[i+1][j-1] - 168.0 * f[i+2][j-1] + 32.0 * f[i+3][j-1] - 3.0 * f[i+4][j-1]) + \
                             672.0 * (3.0 * f[i-4][j+1] - 32.0 * f[i-3][j+1] + 168.0 * f[i-2][j+1] - 672.0 * f[i-1][j+1] + 672.0 * f[i+1][j+1] - 168.0 * f[i+2][j+1] + 32.0 * f[i+3][j+1] - 3.0 * f[i+4][j+1]) - \
                             168.0 * (3.0 * f[i-4][j+2] - 32.0 * f[i-3][j+2] + 168.0 * f[i-2][j+2] - 672.0 * f[i-1][j+2] + 672.0 * f[i+1][j+2] - 168.0 * f[i+2][j+2] + 32.0 * f[i+3][j+2] - 3.0 * f[i+4][j+2]) + \
                              32.0 * (3.0 * f[i-4][j+3] - 32.0 * f[i-3][j+3] + 168.0 * f[i-2][j+3] - 672.0 * f[i-1][j+3] + 672.0 * f[i+1][j+3] - 168.0 * f[i+2][j+3] + 32.0 * f[i+3][j+3] - 3.0 * f[i+4][j+3]) - \
                               3.0 * (3.0 * f[i-4][j+4] - 32.0 * f[i-3][j+4] + 168.0 * f[i-2][j+4] - 672.0 * f[i-1][j+4] + 672.0 * f[i+1][j+4] - 168.0 * f[i+2][j+4] + 32.0 * f[i+3][j+4] - 3.0 * f[i+4][j+4])) / 705600.0;
        }
    }
}

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


    double *arr_x, *arr_y, *arr_f, *arr_fdx, *arr_fdy, *arr_fdxdy;
    arr_x = new double [n_cols];
    arr_y = new double [n_rows];
    arr_f = new double [n_cols * n_rows];
    arr_fdx = new double [n_cols * n_rows];
    arr_fdy = new double [n_cols * n_rows];
    arr_fdxdy = new double [n_cols * n_rows];
    process_input_arrays(eq3_Rs, eq3_Zs, reshaped_Psi, arr_x, arr_y, arr_f,
                        arr_fdx, arr_fdy, arr_fdxdy);
    printf("Checking array alloc\n");
    printf("arr_x[0]=%f\n", arr_x[0]);
    bool ok=true;
    int i,j,k;
    // Checking values of psi
    for (i=0; i<n_rows; i++){
        for (j=0; j<n_cols; j++){
            k = i * n_cols + j;
            if (!cmpd(naive_interp->m_f[i][j], arr_f[k])){
                ok = false;
            }
        }
    }
    if (!ok){
        printf("Psi values differ!\n");
    }
    else{
        printf("Psi values are the same!\n");
    }


    // Checking values of psi dx
    for (i=0; i<n_rows; i++){
        for (j=0; j<n_cols; j++){
            k = i * n_cols + j;
            if (!cmpd(naive_interp->m_fdx[i][j], arr_fdx[k])){
                ok = false;
            }
        }
    }
    if (!ok){
        printf("Psidx values differ!\n");
    }
   else{
        printf("Psidx values are the same!\n");
    }
    ok = true;

    // Checking values of psi
    for (i=0; i<n_rows; i++){
        for (j=0; j<n_cols; j++){
            k = i * n_cols + j;
            if (!cmpd(naive_interp->m_fdy[i][j], arr_fdy[k])){
                ok = false;
            }
        }
    }

    if (!ok){
        printf("Psidy values differ!\n");
    }
    else{
        printf("Psidy values are the same!\n");
    }
    ok = true;

    // Checking values of psi
    for (i=0; i<n_rows; i++){
        for (j=0; j<n_cols; j++){
            k = i * n_cols + j;
            if (!cmpd(naive_interp->m_fdxdy[i][j], arr_fdxdy[k])){
                ok = false;
            }
        }
    }
    if (!ok){
        printf("Psidxdy values differ!\n");
    }
    else{
        printf("Psidxdy values are the same!\n");
    }

    delete [] arr_x;
    delete [] arr_y;
    delete [] arr_f;
    delete [] arr_fdx;
    delete [] arr_fdy;
    delete [] arr_fdxdy;

}