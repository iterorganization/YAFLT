#include <bicubic.hpp>
#include <stdio.h>
#include <cmath>

double clip(double x, double lower, double upper){
    if (x < lower) x = lower;
    if (x > upper) x = upper;
    return x;
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
}

BICUBIC_INTERP::BICUBIC_INTERP(){

};

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

    /*Also set the coefficients to zero*/
    // m_omp_a.clear();
    // m_omp_a.resize(16);
    // m_b.clear();
    // m_b.resize(16);
    // std::fill(m_omp_a.begin(), m_omp_a.end(), 0);
    // std::fill(m_b.begin(), m_b.end(), 0);

    // m_cell_row = -1;
    // m_cell_col = -1;

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
    // o-------->> X
    // |
    // |
    // |
    // |
    // v
    // Y
    //m_fdx is derivative in R direciton
    m_fdx.resize(n_rows, std::vector<double>(n_cols));
    //m_fdy is derivative in Z direction
    m_fdy.resize(n_rows, std::vector<double>(n_cols));

    m_fdxdy.resize(n_rows, std::vector<double>(n_cols));

    // Upper left border values
    m_fdx[0][0] = 0.5 * (m_f[0][1] - m_f[0][0]);
    m_fdy[0][0] = 0.5 * (m_f[1][0] - m_f[0][0]);
    m_fdxdy[0][0] = 0.25 * (m_f[1][1] - m_f[0][1] - m_f[1][0] + m_f[0][0]);

    // Upper right border values
    m_fdx[0][n_cols - 1] = 0.5 * (m_f[0][n_cols - 1] - m_f[0][n_cols - 2]);
    m_fdy[0][n_cols - 1] = 0.5 * (m_f[1][n_cols - 1] - m_f[0][n_cols - 1]);
    m_fdxdy[0][n_cols - 1] = 0.25 * (m_f[1][n_cols - 1] - m_f[0][n_cols - 1] - m_f[1][n_cols - 2] + m_f[0][n_cols - 2]);

    // Lower left border values
    m_fdx[n_rows - 1][0] = 0.5 * (m_f[n_rows - 1][1] - m_f[n_rows - 1][0]);
    m_fdy[n_rows - 1][0] = 0.5 * (m_f[n_rows - 1][0] - m_f[n_rows - 2][0]);
    m_fdxdy[n_rows - 1][0] = 0.25 * (m_f[n_rows - 1][1] - m_f[n_rows - 2][1] - m_f[n_rows - 1][0] + m_f[n_rows - 2][0]);

    // Lower right border values
    m_fdx[n_rows - 1][n_cols - 1] = 0.5 * (m_f[n_rows - 1][n_cols - 1] - m_f[n_rows - 1][n_cols - 2]);
    m_fdy[n_rows - 1][n_cols - 1] = 0.5 * (m_f[n_rows - 1][n_cols - 1] - m_f[n_rows - 2][n_cols - 1]);
    m_fdxdy[n_rows - 1][n_cols - 1] = 0.25 * (m_f[n_rows - 1][n_cols - 1] - m_f[n_rows - 2][n_cols - 1] - m_f[n_rows - 1][n_cols - 2] + m_f[n_rows - 2][n_cols - 2]);


    // Left and right border lines
    for(int i=1; i<n_rows-1; i++){
        // Forward derivative - left side
        m_fdx[i][0] = 0.5 * (m_f[i][1] - m_f[i][0]);
        m_fdy[i][0] = 0.5 * (m_f[i+1][0] - m_f[i-1][0]);
        m_fdxdy[i][0] = 0.25 * (m_f[i+1][1] - m_f[i-1][1] - m_f[i+1][0] + m_f[i-1][0]);

        // Backward derivative - right side
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

    // Centered finite differences.
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

    // Higher order finite difference
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

};

void BICUBIC_INTERP::prepareContainers(int number_of_omp_threads){
    // Prepare the vector for the coefficients.
    m_omp_a.clear();
    m_cell_row.clear();
    m_cell_col.clear();

    m_omp_a.resize(number_of_omp_threads);
    m_cell_row.resize(number_of_omp_threads);
    m_cell_col.resize(number_of_omp_threads);
    for (int i=0; i < number_of_omp_threads; i++){
        m_omp_a[i].resize(16); // No need to set values as they are always set.
        m_cell_row[i] = -1;
        m_cell_col[i] = -1;
    }
}

void BICUBIC_INTERP::rcin(double x, double y, double &out_cell_x, double &out_cell_y, int omp_index){
    /*rcin - Recompute Coefficients If Needed*/
    int ind_x, ind_y;
    double cellx, celly;

    // Get the index
    ind_x = (int) ((x - m_minx)/ m_dx);
    ind_y = (int) ((y - m_miny)/ m_dy);

    // Calculate the position in cell
    cellx = (x - m_x[ind_x]) / m_dx;
#ifndef NDEBUG
    printf("Calculating index and position in cell\n");
    fflush(stdout);
    printf("Now for x\n");
    fflush(stdout);
#endif
    // First clipping.
    // if (cellx > 1) {
    //     cellx = 0;
    //     ind_x = (ind_x + 1) % m_nx;
    // }
    // Upper range clip
    if (ind_x >= m_nx - 1){
        ind_x = ind_x - 1;
        cellx = 1.0;
    }
#ifndef NDEBUG
    printf("cellx=%f\n", cellx);
    printf("nx=%d\n", m_nx);
    fflush(stdout);
    printf("Now for y\n");
    fflush(stdout);
#endif

    celly = (y - m_y[ind_y]) / m_dy;
#ifndef NDEBUG
    printf("celly=%f\n", celly);
    printf("ny=%d\n", m_ny);
    fflush(stdout);
#endif
    // if (celly > 1) {
    //     celly = 0;
    //     ind_y = (ind_y + 1) % m_ny;
    // }
    // Upper range clip
    if (ind_y >= m_ny - 1){
        ind_y = ind_y - 1;
        celly = 1.0;
    }
#ifndef NDEBUG
    printf("Calculated in cell relative position: %f %f\n", cellx, celly);
    fflush(stdout);
    printf("Calculated cell index: %d %d\n", ind_x, ind_y);
    fflush(stdout);
#endif

    out_cell_x = cellx;
    out_cell_y = celly;
    if (ind_x == m_cell_col[omp_index] && ind_y == m_cell_row[omp_index]){
            // We are in the same cell, no recomputation needed.
            return;
        }
    m_cell_col[omp_index] = ind_x;
    m_cell_row[omp_index] = ind_y;
    interpolate(m_cell_row[omp_index], m_cell_col[omp_index], omp_index);
};

void BICUBIC_INTERP::interpolate(int r, int c, int omp_index){
    /*Recalculate the m_omp_a coefficients.*/

    // Calculate the m_omp_a coefficients.

    // Take into account that the array indexes 0,0 starts from upper left
    // corner

    int rp1, rm1;//, rp2;
    int cp1, cm1;//, cp2;

    double b[16];

    rp1 = r + 1;
    cp1 = c + 1;

    // rp2 = r + 2;
    // cp2 = c + 2;

    /*First 4 values is the value of the function in the corner.*/
    b[0] =  m_f[r][c];
    b[1] =  m_f[r][cp1];
    b[2] =  m_f[rp1][c];
    b[3] =  m_f[rp1][cp1];
    /*First partial derivative by x*/
    // b[4] =  0.5*(m_f[rp1][c]   - m_f[rm1][c]);
    // b[5] =  0.5*(m_f[rp2][c]   - m_f[r][c]);
    // b[6] =  0.5*(m_f[rp1][cp1] - m_f[rm1][cp1]);
    // b[7] =  0.5*(m_f[rp2][cp1] - m_f[r][cp1]);
    b[4] =  m_fdx[r][c];
    b[5] =  m_fdx[r][cp1];
    b[6] =  m_fdx[rp1][c];
    b[7] =  m_fdx[rp1][cp1];

    /*First partial derivative by y*/
    // b[8] =  0.5*(m_f[r][cp1]   - m_f[r][cm1]);
    // b[9] =  0.5*(m_f[rp1][cp1] - m_f[rp1][cm1]);
    // b[10] = 0.5*(m_f[r][cp2]   - m_f[r][c]);
    // b[11] = 0.5*(m_f[rp1][cp2] - m_f[rp1][c]);
    b[8] =  m_fdy[r][c];
    b[9] =  m_fdy[r][cp1];
    b[10] = m_fdy[rp1][c];
    b[11] = m_fdy[rp1][cp1];

    /*First partial derivative by xy*/
    // b[12] = 0.25*(m_f[rp1][cp1] - m_f[rm1][cp1] - m_f[rp1][cm1] + m_f[rm1][cm1]);
    // b[13] = 0.25*(m_f[rp2][cp1] - m_f[r][cp1] - m_f[rp2][cm1] + m_f[r][cm1]);
    // b[14] = 0.25*(m_f[rp1][cp2] - m_f[rm1][cp2] - m_f[rp1][c] + m_f[rm1][c]);
    // b[15] = 0.25*(m_f[rp2][cp2] - m_f[r][cp2] - m_f[rp2][c] + m_f[r][c]);
    b[12] = m_fdxdy[r][c];
    b[13] = m_fdxdy[r][cp1];
    b[14] = m_fdxdy[rp1][c];
    b[15] = m_fdxdy[rp1][cp1];

    // for(int i=0;i<16;i++){
    //     m_omp_a[omp_index][i] = 0;
    //     for(int j=0;j<16;j++){
    //         m_omp_a[omp_index][i] = m_omp_a[omp_index][i] + bicubic_A[i][j] * b[j];
    //     }
    // }

    m_omp_a[omp_index][0] = b[0];
    m_omp_a[omp_index][1] = b[4];
    m_omp_a[omp_index][2] = -3.0 * b[0]+3.0 * b[1]-2.0 * b[4]-1.0 * b[5];
    m_omp_a[omp_index][3] = 2.0 * b[0]-2.0 * b[1]+b[4]+b[5];
    m_omp_a[omp_index][4] = b[8];
    m_omp_a[omp_index][5] = b[12];
    m_omp_a[omp_index][6] = -3.0 * b[8]+3.0 * b[9]-2.0 * b[12]-1.0 * b[13];
    m_omp_a[omp_index][7] = 2.0 * b[8]-2.0 * b[9]+b[12]+b[13];
    m_omp_a[omp_index][8] = -3.0 * b[0]+3.0 * b[2]-2.0 * b[8]-1.0 * b[10];
    m_omp_a[omp_index][9] = -3.0 * b[4]+3.0 * b[6]-2.0 * b[12]-1.0 * b[14];
    m_omp_a[omp_index][10] = 9.0 * b[0]-9.0 * b[1]-9.0 * b[2]+9.0 * b[3]+6.0 * b[4]+3.0 * b[5]-6.0 * b[6]-3.0 * b[7]+6.0 * b[8]-6.0 * b[9]+3.0 * b[10]-3.0 * b[11]+4.0 * b[12]+2.0 * b[13]+2.0 * b[14]+b[15];
    m_omp_a[omp_index][11] = -6.0 * b[0]+6.0 * b[1]+6.0 * b[2]-6.0 * b[3]-3.0 * b[4]-3.0 * b[5]+3.0 * b[6]+3.0 * b[7]-4.0 * b[8]+4.0 * b[9]-2.0 * b[10]+2.0 * b[11]-2.0 * b[12]-2.0 * b[13]-1.0 * b[14]-1.0 * b[15];
    m_omp_a[omp_index][12] = 2.0 * b[0]-2.0 * b[2]+b[8]+b[10];
    m_omp_a[omp_index][13] = 2.0 * b[4]-2.0 * b[6]+b[12]+b[14];
    m_omp_a[omp_index][14] = -6.0 * b[0]+6.0 * b[1]+6.0 * b[2]-6.0 * b[3]-4.0 * b[4]-2.0 * b[5]+4.0 * b[6]+2.0 * b[7]-3.0 * b[8]+3.0 * b[9]-3.0 * b[10]+3.0 * b[11]-2.0 * b[12]-1.0 * b[13]-2.0 * b[14]-1.0 * b[15];
    m_omp_a[omp_index][15] = 4.0 * b[0]-4.0 * b[1]-4.0 * b[2]+4.0 * b[3]+2.0 * b[4]+2.0 * b[5]-2.0 * b[6]-2.0 * b[7]+2.0 * b[8]-2.0 * b[9]+2.0 * b[10]-2.0 * b[11]+b[12]+b[13]+b[14]+b[15];
};

void BICUBIC_INTERP::getValues(double x, double y, double &val, double &valdx,
                               double &valdy, int omp_index){
    // Function for obtaining values inside the X, Y domain.
#ifndef NDEBUG
    printf("GetValues %f %f\n", x, y);
    printf("omp_index %d\n", omp_index);
#endif
    val = 0;
    valdx = 0;
    valdy = 0;
#ifndef NDEBUG
    printf("Clipping\n");
#endif
    x = clip(x, m_minx, m_maxx);
    y = clip(y, m_miny, m_maxy);
#ifndef NDEBUG
    printf("After clip: %f %f\n", x, y);
#endif
    double cell_x, cell_y;
    rcin(x, y, cell_x, cell_y, omp_index);
#ifndef NDEBUG
    printf("Rcin successfull\n");
#endif
    double cell_x2 = cell_x * cell_x;
    double cell_y2 = cell_y * cell_y;

    val = m_omp_a[omp_index][ 0]  + /*a[0, 0]*/ \
          m_omp_a[omp_index][ 1] * cell_x  + /*a[1, 0]*/ \
          m_omp_a[omp_index][ 2] * cell_x2  + /*a[2, 0]*/ \
          m_omp_a[omp_index][ 3] * cell_x2 * cell_x  + /*a[3, 0]*/ \
          m_omp_a[omp_index][ 4] * cell_y + /*a[0, 1]*/ \
          m_omp_a[omp_index][ 5] * cell_x * cell_y + /*a[1, 1]*/ \
          m_omp_a[omp_index][ 6] * cell_x2 * cell_y + /*a[2, 1]*/ \
          m_omp_a[omp_index][ 7] * cell_x2 * cell_x * cell_y + /*a[3, 1]*/ \
          m_omp_a[omp_index][ 8] * cell_y * cell_y + /*a[0, 2]*/ \
          m_omp_a[omp_index][ 9] * cell_x * cell_y2 + /*a[1, 2]*/ \
          m_omp_a[omp_index][10] * cell_x2 * cell_y2 + /*a[2, 2]*/ \
          m_omp_a[omp_index][11] * cell_x2 * cell_x * cell_y2 + /*a[3, 2]*/ \
          m_omp_a[omp_index][12] * cell_y2 * cell_y + /*a[0, 3]*/ \
          m_omp_a[omp_index][13] * cell_x * cell_y2 * cell_y + /*a[1, 3]*/ \
          m_omp_a[omp_index][14] * cell_x2 * cell_y2 * cell_y + /*a[2, 3]*/ \
          m_omp_a[omp_index][15] * cell_x2 * cell_x * cell_y2 * cell_y /*a[3, 3]*/;

    // By definition this should be valdx, but... the way data is presented
    // it is valdy
    valdx = /*a[0, 0]=0*/ \
            m_omp_a[omp_index][ 1]  + /*a[1, 0]*/ \
            m_omp_a[omp_index][ 2] * 2 * cell_x  + /*a[2, 0]*/ \
            m_omp_a[omp_index][ 3] * 3 * cell_x2  + /*a[3, 0]*/ \
            /*a[0, 1]=0*/ \
            m_omp_a[omp_index][ 5] * cell_y + /*a[1, 1]*/ \
            m_omp_a[omp_index][ 6] * 2 * cell_x * cell_y + /*a[2, 1]*/ \
            m_omp_a[omp_index][ 7] * 3 * cell_x2 * cell_y + /*a[3, 1]*/ \
            /*a[0, 2]=0*/ \
            m_omp_a[omp_index][ 9] * cell_y * cell_y + /*a[1, 2]*/ \
            m_omp_a[omp_index][10] * 2 * cell_x * cell_y * cell_y + /*a[2, 2]*/ \
            m_omp_a[omp_index][11] * 3 * cell_x2 * cell_y2 + /*a[3, 2]*/ \
            /*a[0, 3]=0*/ \
            m_omp_a[omp_index][13] * cell_y2 * cell_y + /*a[1, 3]*/ \
            m_omp_a[omp_index][14] * 2 * cell_x * cell_y2 * cell_y + /*a[2, 3]*/ \
            m_omp_a[omp_index][15] * 3 * cell_x2 * cell_y2 * cell_y /*a[3, 3]*/;

    // Same here, it should be valdy, but it is valdx
    valdy = /*a[0, 0]=0*/ \
            /*a[1, 0]=0*/ \
            /*a[2, 0]=0*/ \
            /*a[3, 0]=0*/ \
            m_omp_a[omp_index][ 4]  + /*a[0, 1]*/ \
            m_omp_a[omp_index][ 5] * cell_x  + /*a[1, 1]*/ \
            m_omp_a[omp_index][ 6] * cell_x2  + /*a[2, 1]*/ \
            m_omp_a[omp_index][ 7] * cell_x2 * cell_x  + /*a[3, 1]*/ \
            m_omp_a[omp_index][ 8] * 2 * cell_y + /*a[0, 2]*/ \
            m_omp_a[omp_index][ 9] * cell_x * 2 * cell_y + /*a[1, 2]*/ \
            m_omp_a[omp_index][10] * cell_x2 * 2 * cell_y + /*a[2, 2]*/ \
            m_omp_a[omp_index][11] * cell_x2 * cell_x * 2 * cell_y + /*a[3, 2]*/ \
            m_omp_a[omp_index][12] * 3 * cell_y * cell_y + /*a[0, 3]*/ \
            m_omp_a[omp_index][13] * cell_x * 3 * cell_y2 + /*a[1, 3]*/ \
            m_omp_a[omp_index][14] * cell_x2 * 3 * cell_y2 + /*a[2, 3]*/ \
            m_omp_a[omp_index][15] * cell_x2 * cell_x * 3 * cell_y2 /*a[3, 3]*/;

    // Apply rectilinear factors

    // Originally it should be like this until I figure out the indexing
    // valdx = valdx / m_dx;
    // valdy = valdy / m_dy;
    valdx = valdx / m_dx;
    valdy = valdy / m_dy;

    // Equivalent to the following.
    // int index =0;
    // for (int j=0; j<4;j++){
    //     for (int i=0; i<4;i++){
    //         index = j * 4 + i;
    //         val = val + m_omp_a[index] * pow(m_cell_x, i) * pow(m_celly, j);
    //         if (i != 0){
    //             valdx = valdx + m_omp_a[index] * i * pow(m_cell_x, i - 1) * pow(m_celly, j);
    //         }
    //         if (j != 0){
    //             valdy = valdy + m_omp_a[index] * j * pow(m_cell_x, i) * pow(m_celly, j - 1);
    //         }
    //     }
    // }

};
