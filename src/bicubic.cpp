#include "bicubic.hpp"
#include <stdio.h>
#include <cmath>

double __clip(double x, double lower, double upper){
    const double t = x < lower ? lower : x;
    return t > upper ? upper : t;
}


void BICUBIC_INTERP::setArrays(std::vector<double> x, std::vector<double> y,
                               std::vector<std::vector<double>> f){
    /*Copy the input data for interpolation*/

    /*Reset the vectors*/
    m_x.clear();
    m_y.clear();
    m_f.clear();
    m_fdx.clear();
    m_fdy.clear();
    m_fdxdy.clear();
    m_a.clear();

    /*Copy the vectors.*/
    m_x = x;
    m_y = y;
    m_f = f;

    /*Now set the mins, maxes and deltas*/
    m_nx = m_x.size();
    m_ny = m_y.size();

    m_minx = m_x[0];
    m_maxx = m_x[m_nx-1];
    // m_dx = (m_maxx - m_minx) / m_nx;
    m_dx = m_x[1] - m_x[0];

    m_miny = m_y[0];
    m_maxy = m_y[m_ny-1];
    // m_dy = (m_maxy - m_miny) / m_ny;
    m_dy = m_y[1] - m_y[0];


    int n_rows, n_cols;
    n_rows = m_ny;
    n_cols = m_nx;

    // Now calcuate the derivatives
    // The following is the direction of Rows and Columns!!!
    // o---------> X (R)
    // |
    // |
    // |
    // |
    // v
    // Y (Z)
    // m_fdx is the partial derivative in R direciton
    m_fdx.resize(n_rows, std::vector<double>(n_cols));
    // m_fdy is the partial derivative in Z direction
    m_fdy.resize(n_rows, std::vector<double>(n_cols));
    // m_fdxdy is the partial mixed derivative
    m_fdxdy.resize(n_rows, std::vector<double>(n_cols));

    // m_a holds the interpolation constants
    m_a.resize(n_rows*n_cols*16);

    // Upper left border values - forward finite difference
    m_fdx[0][0] = 0.5 * (m_f[0][1] - m_f[0][0]);
    m_fdy[0][0] = 0.5 * (m_f[1][0] - m_f[0][0]);
    m_fdxdy[0][0] = 0.25 * (m_f[1][1] - m_f[0][1] - m_f[1][0] + m_f[0][0]);

    // Upper right border values - backward finite difference
    m_fdx[0][n_cols - 1] = 0.5 * (m_f[0][n_cols - 1] - m_f[0][n_cols - 2]);
    m_fdy[0][n_cols - 1] = 0.5 * (m_f[1][n_cols - 1] - m_f[0][n_cols - 1]);
    m_fdxdy[0][n_cols - 1] = 0.25 * (m_f[1][n_cols - 1] - m_f[0][n_cols - 1] - m_f[1][n_cols - 2] + m_f[0][n_cols - 2]);

    // Lower left border values - forward finite difference
    m_fdx[n_rows - 1][0] = 0.5 * (m_f[n_rows - 1][1] - m_f[n_rows - 1][0]);
    m_fdy[n_rows - 1][0] = 0.5 * (m_f[n_rows - 1][0] - m_f[n_rows - 2][0]);
    m_fdxdy[n_rows - 1][0] = 0.25 * (m_f[n_rows - 1][1] - m_f[n_rows - 2][1] - m_f[n_rows - 1][0] + m_f[n_rows - 2][0]);

    // Lower right border values - backward finite difference
    m_fdx[n_rows - 1][n_cols - 1] = 0.5 * (m_f[n_rows - 1][n_cols - 1] - m_f[n_rows - 1][n_cols - 2]);
    m_fdy[n_rows - 1][n_cols - 1] = 0.5 * (m_f[n_rows - 1][n_cols - 1] - m_f[n_rows - 2][n_cols - 1]);
    m_fdxdy[n_rows - 1][n_cols - 1] = 0.25 * (m_f[n_rows - 1][n_cols - 1] - m_f[n_rows - 2][n_cols - 1] - m_f[n_rows - 1][n_cols - 2] + m_f[n_rows - 2][n_cols - 2]);

    // Left and right border lines
    for(int i=1; i<n_rows-1; i++){
        // Forward difference - left side
        m_fdx[i][0] = 0.5 * (m_f[i][1] - m_f[i][0]);
        // Centered difference
        m_fdy[i][0] = 0.5 * (m_f[i+1][0] - m_f[i-1][0]);
        m_fdxdy[i][0] = 0.25 * (m_f[i+1][1] - m_f[i-1][1] - m_f[i+1][0] + m_f[i-1][0]);
        // Backward difference - right side
        m_fdx[i][n_cols - 1]   = 0.5 * (m_f[i][n_cols - 1] - m_f[i][n_cols - 2]);
        m_fdy[i][n_cols - 1]   = 0.5 * (m_f[i+1][n_cols - 1] - m_f[i-1][n_cols - 1]);
        m_fdxdy[i][n_cols - 1] = 0.25 * (m_f[i+1][n_cols - 1] - m_f[i-1][n_cols - 1] - m_f[i+1][n_cols - 2] + m_f[i-1][n_cols - 2]);
    }

    // Upper and bottom border lines
    for(int j=1; j<n_cols-1;j++){
        // Upper side
        m_fdx[0][j] = 0.5 * (m_f[0][j+1] - m_f[0][j-1]);
        m_fdy[0][j] = 0.5 * (m_f[1][j] - m_f[0][j]);
        m_fdxdy[0][j] = 0.25 * (m_f[1][j+1] - m_f[0][j+1] - m_f[1][j-1] + m_f[0][j-1]);

        // Lower side
        m_fdx[n_rows-1][j] = 0.5 * (m_f[n_rows-1][j+1] - m_f[n_rows-1][j-1]);
        m_fdy[n_rows-1][j] = 0.5 * (m_f[n_rows-1][j] - m_f[n_rows-2][j]);
        m_fdxdy[n_rows-1][j] = 0.25 * (m_f[n_rows-1][j+1] - m_f[n_rows-2][j+1] - m_f[n_rows-1][j-1] + m_f[n_rows-2][j-1]);
    }

    // Centered finite differences without the step distance!
    // Now for the rest of the space
    // This is lazy programming, basically, this order of derivatives should
    // only be on the 2 row/col width band and not the whole area. But no
    // problems since it is overwritten by the higher order equations.
    for(int i=1; i<n_rows-1;i++){
        for(int j=1; j<n_cols-1;j++){
            // m_fdx[i][j] = 0.5 *(m_f[i][j+1] - m_f[i][j-1]);
            m_fdx[i][j] = 0.5 *(m_f[i][j+1] - m_f[i][j-1]);
            m_fdy[i][j] = 0.5 *(m_f[i+1][j] - m_f[i-1][j]);
            m_fdxdy[i][j] = 0.25 * (m_f[i+1][j+1] - m_f[i-1][j+1] - m_f[i+1][j-1] + m_f[i-1][j-1]);
        }
    }

    // 4th order centered difference approximations
    for(int i=2; i<n_rows-2;i++){
        for(int j=2; j<n_cols-2;j++){
            m_fdx[i][j] = (-m_f[i][j+2] + 8.0 * m_f[i][j+1] - 8.0 * m_f[i][j-1] + m_f[i][j-2]) / 12.0;
            m_fdy[i][j] = (-m_f[i+2][j] + 8.0 * m_f[i+1][j] - 8.0 * m_f[i-1][j] + m_f[i-2][j]) / 12.0;
            m_fdxdy[i][j] = (m_f[i+2][j+2] - 8.0 * m_f[i+1][j+2] + 8.0 * m_f[i-1][j+2] - m_f[i-2][j+2] + \
                             8.0 * (-m_f[i+2][j+1] + 8.0 * m_f[i+1][j+1] - 8.0 * m_f[i-1][j+1] + m_f[i-2][j+1]) - \
                             8.0 * (-m_f[i+2][j-1] + 8.0 * m_f[i+1][j-1] - 8.0 * m_f[i-1][j-1] + m_f[i-2][j-1]) - \
                             m_f[i+2][j-2] + 8.0 * m_f[i+1][j-2] - 8.0 * m_f[i-1][j-2] + m_f[i-2][j-2]) / 144.0;
        }
    }

    // 6th order centered difference approximations
    for(int i=3; i<n_rows-3;i++){
        for(int j=3; j<n_cols-3;j++){
            m_fdx[i][j] = (-m_f[i][j-3] + 9.0 * m_f[i][j-2] - 45.0 * m_f[i][j-1] + 45.0 * m_f[i][j+1] - 9.0 * m_f[i][j+2] + m_f[i][j+3]) / 60.0;
            m_fdy[i][j] = (-m_f[i-3][j] + 9.0 * m_f[i-2][j] - 45.0 * m_f[i-1][j] + 45.0 * m_f[i+1][j] - 9.0 * m_f[i+2][j] + m_f[i+3][j]) / 60.0;

            m_fdxdy[i][j] = (-1.0 * (-m_f[i-3][j-3] + 9.0 * m_f[i-2][j-3] - 45.0 * m_f[i-1][j-3] + 45.0 * m_f[i+1][j-3] - 9.0 * m_f[i+2][j-3] + m_f[i+3][j-3]) + \
                              9.0 * (-m_f[i-3][j-2] + 9.0 * m_f[i-2][j-2] - 45.0 * m_f[i-1][j-2] + 45.0 * m_f[i+1][j-2] - 9.0 * m_f[i+2][j-2] + m_f[i+3][j-2]) + \
                            -45.0 * (-m_f[i-3][j-1] + 9.0 * m_f[i-2][j-1] - 45.0 * m_f[i-1][j-1] + 45.0 * m_f[i+1][j-1] - 9.0 * m_f[i+2][j-1] + m_f[i+3][j-1]) + \
                             45.0 * (-m_f[i-3][j+1] + 9.0 * m_f[i-2][j+1] - 45.0 * m_f[i-1][j+1] + 45.0 * m_f[i+1][j+1] - 9.0 * m_f[i+2][j+1] + m_f[i+3][j+1]) + \
                             -9.0 * (-m_f[i-3][j+2] + 9.0 * m_f[i-2][j+2] - 45.0 * m_f[i-1][j+2] + 45.0 * m_f[i+1][j+2] - 9.0 * m_f[i+2][j+2] + m_f[i+3][j+2]) + \
                              1.0 * (-m_f[i-3][j+3] + 9.0 * m_f[i-2][j+3] - 45.0 * m_f[i-1][j+3] + 45.0 * m_f[i+1][j+3] - 9.0 * m_f[i+2][j+3] + m_f[i+3][j+3]))/ 3600.0;
        }
    }

    // 8th order centered difference approximation
    for(int i=4; i<n_rows-4;i++){
        for(int j=4; j<n_cols-4;j++){
            m_fdx[i][j] = (3.0 * m_f[i][j-4] - 32.0 * m_f[i][j-3] + 168.0 * m_f[i][j-2] - 672.0 * m_f[i][j-1] + 672.0 * m_f[i][j+1] - 168.0 * m_f[i][j+2] + 32.0 * m_f[i][j+3] - 3.0 * m_f[i][j+4]) / 840.0;
            m_fdy[i][j] = (3.0 * m_f[i-4][j] - 32.0 * m_f[i-3][j] + 168.0 * m_f[i-2][j] - 672.0 * m_f[i-1][j] + 672.0 * m_f[i+1][j] - 168.0 * m_f[i+2][j] + 32.0 * m_f[i+3][j] - 3.0 * m_f[i+4][j]) / 840.0;

            m_fdxdy[i][j] = (3.0 * (3.0 * m_f[i-4][j-4] - 32.0 * m_f[i-3][j-4] + 168.0 * m_f[i-2][j-4] - 672.0 * m_f[i-1][j-4] + 672.0 * m_f[i+1][j-4] - 168.0 * m_f[i+2][j-4] + 32.0 * m_f[i+3][j-4] - 3.0 * m_f[i+4][j-4]) - \
                            32.0 * (3.0 * m_f[i-4][j-3] - 32.0 * m_f[i-3][j-3] + 168.0 * m_f[i-2][j-3] - 672.0 * m_f[i-1][j-3] + 672.0 * m_f[i+1][j-3] - 168.0 * m_f[i+2][j-3] + 32.0 * m_f[i+3][j-3] - 3.0 * m_f[i+4][j-3]) + \
                           168.0 * (3.0 * m_f[i-4][j-2] - 32.0 * m_f[i-3][j-2] + 168.0 * m_f[i-2][j-2] - 672.0 * m_f[i-1][j-2] + 672.0 * m_f[i+1][j-2] - 168.0 * m_f[i+2][j-2] + 32.0 * m_f[i+3][j-2] - 3.0 * m_f[i+4][j-2]) - \
                           672.0 * (3.0 * m_f[i-4][j-1] - 32.0 * m_f[i-3][j-1] + 168.0 * m_f[i-2][j-1] - 672.0 * m_f[i-1][j-1] + 672.0 * m_f[i+1][j-1] - 168.0 * m_f[i+2][j-1] + 32.0 * m_f[i+3][j-1] - 3.0 * m_f[i+4][j-1]) + \
                           672.0 * (3.0 * m_f[i-4][j+1] - 32.0 * m_f[i-3][j+1] + 168.0 * m_f[i-2][j+1] - 672.0 * m_f[i-1][j+1] + 672.0 * m_f[i+1][j+1] - 168.0 * m_f[i+2][j+1] + 32.0 * m_f[i+3][j+1] - 3.0 * m_f[i+4][j+1]) - \
                           168.0 * (3.0 * m_f[i-4][j+2] - 32.0 * m_f[i-3][j+2] + 168.0 * m_f[i-2][j+2] - 672.0 * m_f[i-1][j+2] + 672.0 * m_f[i+1][j+2] - 168.0 * m_f[i+2][j+2] + 32.0 * m_f[i+3][j+2] - 3.0 * m_f[i+4][j+2]) + \
                            32.0 * (3.0 * m_f[i-4][j+3] - 32.0 * m_f[i-3][j+3] + 168.0 * m_f[i-2][j+3] - 672.0 * m_f[i-1][j+3] + 672.0 * m_f[i+1][j+3] - 168.0 * m_f[i+2][j+3] + 32.0 * m_f[i+3][j+3] - 3.0 * m_f[i+4][j+3]) - \
                             3.0 * (3.0 * m_f[i-4][j+4] - 32.0 * m_f[i-3][j+4] + 168.0 * m_f[i-2][j+4] - 672.0 * m_f[i-1][j+4] + 672.0 * m_f[i+1][j+4] - 168.0 * m_f[i+2][j+4] + 32.0 * m_f[i+3][j+4] - 3.0 * m_f[i+4][j+4])) / 705600.0;
        }
    }

    int rp1;
    int cp1;
    int offset;
    double b[16];

    // r - row, c - col
    for(int r=0; r<n_rows-1;r++){
        for(int c=0; c<n_cols-1;c++){
            rp1 = r + 1;
            cp1 = c + 1;
            // First 4 values are the values of the function.
            //std::vector<std::vector<double>> p{m_f};
            b[0] =  m_f[r][c];
            b[1] =  m_f[r][cp1];
            b[2] =  m_f[rp1][c];
            b[3] =  m_f[rp1][cp1];
            // First partial derivative by x
            b[4] =  m_fdx[r][c];
            b[5] =  m_fdx[r][cp1];
            b[6] =  m_fdx[rp1][c];
            b[7] =  m_fdx[rp1][cp1];
            /*First partial derivative by y*/
            b[8] =  m_fdy[r][c];
            b[9] =  m_fdy[r][cp1];
            b[10] = m_fdy[rp1][c];
            b[11] = m_fdy[rp1][cp1];
            /*Mixed partial derivative by xy*/
            b[12] = m_fdxdy[r][c];
            b[13] = m_fdxdy[r][cp1];
            b[14] = m_fdxdy[rp1][c];
            b[15] = m_fdxdy[rp1][cp1];
            offset = (r*n_cols + c)*16;
            m_a[offset + 0] = b[0];
            m_a[offset + 1] = b[4];
            m_a[offset + 2] = -3.0 * b[0]+3.0 * b[1]-2.0 * b[4]-1.0 * b[5];
            m_a[offset + 3] = 2.0 * b[0]-2.0 * b[1]+b[4]+b[5];
            m_a[offset + 4] = b[8];
            m_a[offset + 5] = b[12];
            m_a[offset + 6] = -3.0 * b[8]+3.0 * b[9]-2.0 * b[12]-1.0 * b[13];
            m_a[offset + 7] = 2.0 * b[8]-2.0 * b[9]+b[12]+b[13];
            m_a[offset + 8] = -3.0 * b[0]+3.0 * b[2]-2.0 * b[8]-1.0 * b[10];
            m_a[offset + 9] = -3.0 * b[4]+3.0 * b[6]-2.0 * b[12]-1.0 * b[14];
            m_a[offset + 10] = 9.0 * b[0]-9.0 * b[1]-9.0 * b[2]+9.0 * b[3]+6.0 * b[4]+3.0 * b[5]-6.0 * b[6]-3.0 * b[7]+6.0 * b[8]-6.0 * b[9]+3.0 * b[10]-3.0 * b[11]+4.0 * b[12]+2.0 * b[13]+2.0 * b[14]+b[15];
            m_a[offset + 11] = -6.0 * b[0]+6.0 * b[1]+6.0 * b[2]-6.0 * b[3]-3.0 * b[4]-3.0 * b[5]+3.0 * b[6]+3.0 * b[7]-4.0 * b[8]+4.0 * b[9]-2.0 * b[10]+2.0 * b[11]-2.0 * b[12]-2.0 * b[13]-1.0 * b[14]-1.0 * b[15];
            m_a[offset + 12] = 2.0 * b[0]-2.0 * b[2]+b[8]+b[10];
            m_a[offset + 13] = 2.0 * b[4]-2.0 * b[6]+b[12]+b[14];
            m_a[offset + 14] = -6.0 * b[0]+6.0 * b[1]+6.0 * b[2]-6.0 * b[3]-4.0 * b[4]-2.0 * b[5]+4.0 * b[6]+2.0 * b[7]-3.0 * b[8]+3.0 * b[9]-3.0 * b[10]+3.0 * b[11]-2.0 * b[12]-1.0 * b[13]-2.0 * b[14]-1.0 * b[15];
            m_a[offset + 15] = 4.0 * b[0]-4.0 * b[1]-4.0 * b[2]+4.0 * b[3]+2.0 * b[4]+2.0 * b[5]-2.0 * b[6]-2.0 * b[7]+2.0 * b[8]-2.0 * b[9]+2.0 * b[10]-2.0 * b[11]+b[12]+b[13]+b[14]+b[15];

        }
    }

};

// void BICUBIC_INTERP::getValues(BI_DATA *context){
void BICUBIC_INTERP::getValues(BI_DATA *context){
    // Function for obtaining values inside the X, Y domain.
#ifndef NDEBUG
    printf("GetValues %f %f\n", context->r, context->z);
#endif
    double x,y;

    /// Output
    context->val = 0;
    context->valdx = 0;
    context->valdy = 0;
#ifndef NDEBUG
    printf("Clipping\n");
#endif
    x = __clip(context->r, context->minx, context->maxx);
    y = __clip(context->z, context->miny, context->maxy);

    // Calculate relative cell position
    int ind_x, ind_y;
    double cell_x, cell_y;

    // Get the index
    ind_x = (int) ((x - context->minx)/ context->dx);
    ind_y = (int) ((y - context->miny)/ context->dy);

    // Calculate the position in cell
    cell_x = (x - m_x[ind_x]) / context->dx;
    // Upper range clip
    if (ind_x >= context->nx - 1){
        ind_x = context->nx - 2;
        cell_x = 1.0;
    }

    // Calculate the position in cell
    cell_y = (y - m_y[ind_y]) / context->dy;
    // Upper range clip
    if (ind_y >= context->ny - 1){
        ind_y = context->ny - 2;
        cell_y = 1.0;
    }

    // if (ind_x != context->col || ind_y != context->row){
    //     // We are in the same cell, no recomputation needed.
    //     context->col = ind_x;
    //     context->row = ind_y;
    //     interpolate(context->row, context->col, context->a);
    // }

    double cell_x2 = cell_x * cell_x;
    double cell_y2 = cell_y * cell_y;

    double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    int offset;
    // a0 = context->a[0];
    // a1 = context->a[1];
    // a2 = context->a[2];
    // a3 = context->a[3];
    // a4 = context->a[4];
    // a5 = context->a[5];
    // a6 = context->a[6];
    // a7 = context->a[7];
    // a8 = context->a[8];
    // a9 = context->a[9];
    // a10 = context->a[10];
    // a11 = context->a[11];
    // a12 = context->a[12];
    // a13 = context->a[13];
    // a14 = context->a[14];
    // a15 = context->a[15];
    offset = (ind_y * m_nx + ind_x)*16;
    a0 = m_a[offset + 0];
    a1 = m_a[offset + 1];
    a2 = m_a[offset + 2];
    a3 = m_a[offset + 3];
    a4 = m_a[offset + 4];
    a5 = m_a[offset + 5];
    a6 = m_a[offset + 6];
    a7 = m_a[offset + 7];
    a8 = m_a[offset + 8];
    a9 = m_a[offset + 9];
    a10 = m_a[offset + 10];
    a11 = m_a[offset + 11];
    a12 = m_a[offset + 12];
    a13 = m_a[offset + 13];
    a14 = m_a[offset + 14];
    a15 = m_a[offset + 15];

    double val = 0.0;
    val += a0  + /*a[0, 0]*/ \
          a1 * cell_x  + /*a[1, 0]*/ \
          a2 * cell_x2  + /*a[2, 0]*/ \
          a3 * cell_x2 * cell_x  + /*a[3, 0]*/ \
          a4 * cell_y + /*a[0, 1]*/ \
          a5 * cell_x * cell_y + /*a[1, 1]*/ \
          a6 * cell_x2 * cell_y + /*a[2, 1]*/ \
          a7 * cell_x2 * cell_x * cell_y + /*a[3, 1]*/ \
          a8 * cell_y * cell_y + /*a[0, 2]*/ \
          a9 * cell_x * cell_y2 + /*a[1, 2]*/ \
          a10 * cell_x2 * cell_y2 + /*a[2, 2]*/ \
          a11 * cell_x2 * cell_x * cell_y2 + /*a[3, 2]*/ \
          a12 * cell_y2 * cell_y + /*a[0, 3]*/ \
          a13 * cell_x * cell_y2 * cell_y + /*a[1, 3]*/ \
          a14 * cell_x2 * cell_y2 * cell_y + /*a[2, 3]*/ \
          a15 * cell_x2 * cell_x * cell_y2 * cell_y /*a[3, 3]*/;
    context->val = val;

    // Get partial derivative dx
    context->valdx = /*a[0, 0]=0*/ \
            a1  + /*a[1, 0]*/ \
            a2 * 2 * cell_x  + /*a[2, 0]*/ \
            a3 * 3 * cell_x2  + /*a[3, 0]*/ \
            /*a[0, 1]=0*/ \
            a5 * cell_y + /*a[1, 1]*/ \
            a6 * 2 * cell_x * cell_y + /*a[2, 1]*/ \
            a7 * 3 * cell_x2 * cell_y + /*a[3, 1]*/ \
            /*a[0, 2]=0*/ \
            a9 * cell_y * cell_y + /*a[1, 2]*/ \
            a10 * 2 * cell_x * cell_y * cell_y + /*a[2, 2]*/ \
            a11 * 3 * cell_x2 * cell_y2 + /*a[3, 2]*/ \
            /*a[0, 3]=0*/ \
            a13 * cell_y2 * cell_y + /*a[1, 3]*/ \
            a14 * 2 * cell_x * cell_y2 * cell_y + /*a[2, 3]*/ \
            a15 * 3 * cell_x2 * cell_y2 * cell_y /*a[3, 3]*/;

    // Get partial derivative dy
    context->valdy = /*a[0, 0]=0*/ \
            /*a[1, 0]=0*/ \
            /*a[2, 0]=0*/ \
            /*a[3, 0]=0*/ \
            a4  + /*a[0, 1]*/ \
            a5 * cell_x  + /*a[1, 1]*/ \
            a6 * cell_x2  + /*a[2, 1]*/ \
            a7 * cell_x2 * cell_x  + /*a[3, 1]*/ \
            a8 * 2 * cell_y + /*a[0, 2]*/ \
            a9 * cell_x * 2 * cell_y + /*a[1, 2]*/ \
            a10 * cell_x2 * 2 * cell_y + /*a[2, 2]*/ \
            a11 * cell_x2 * cell_x * 2 * cell_y + /*a[3, 2]*/ \
            a12 * 3 * cell_y * cell_y + /*a[0, 3]*/ \
            a13 * cell_x * 3 * cell_y2 + /*a[1, 3]*/ \
            a14 * cell_x2 * 3 * cell_y2 + /*a[2, 3]*/ \
            a15 * cell_x2 * cell_x * 3 * cell_y2 /*a[3, 3]*/;

    // Apply rectilinear factors
    context->valdx = context->valdx / context->dx;
    context->valdy = context->valdy / context->dy;

    // Equivalent to the following.
    // int index =0;
    // for (int j=0; j<4;j++){
    //     for (int i=0; i<4;i++){
    //         index = j * 4 + i;
    //         val = val + context->a[index] * pow(m_cell_x, i) * pow(m_celly, j);
    //         if (i != 0){
    //             valdx = valdx + context->a[index] * i * pow(m_cell_x, i - 1) * pow(m_celly, j);
    //         }
    //         if (j != 0){
    //             valdy = valdy + context->a[index] * j * pow(m_cell_x, i) * pow(m_celly, j - 1);
    //         }
    //         if (i != 0 && j != 0){
    //             valdy = valdy + context->a[index] * i * j * pow(m_cell_x, i - 1) * pow(m_celly, j - 1);
    //         }
    //     }
    // }

};
