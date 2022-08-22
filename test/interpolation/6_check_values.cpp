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