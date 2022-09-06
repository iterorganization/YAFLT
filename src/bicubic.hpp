#ifndef BICUBIC_H
#define BICUBIC_H

#include <vector>

static const double bicubic_A[16][16]={
{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,},
{0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,},
{-3.0,3.0,0.0,0.0,-2.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,},
{2.0,-2.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-3.0,3.0,0.0,0.0,-2.0,-1.0,0.0,0.0,},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,-2.0,0.0,0.0,1.0,1.0,0.0,0.0,},
{-3.0,0.0,3.0,0.0,0.0,0.0,0.0,0.0,-2.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,},
{0.0,0.0,0.0,0.0,-3.0,0.0,3.0,0.0,0.0,0.0,0.0,0.0,-2.0,0.0,-1.0,0.0,},
{9.0,-9.0,-9.0,9.0,6.0,3.0,-6.0,-3.0,6.0,-6.0,3.0,-3.0,4.0,2.0,2.0,1.0,},
{-6.0,6.0,6.0,-6.0,-3.0,-3.0,3.0,3.0,-4.0,4.0,-2.0,2.0,-2.0,-2.0,-1.0,-1.0,},
{2.0,0.0,-2.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,},
{0.0,0.0,0.0,0.0,2.0,0.0,-2.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,},
{-6.0,6.0,6.0,-6.0,-4.0,-2.0,4.0,2.0,-3.0,3.0,-3.0,3.0,-2.0,-1.0,-2.0,-1.0,},
{4.0,-4.0,-4.0,4.0,2.0,2.0,-2.0,-2.0,2.0,-2.0,2.0,-2.0,1.0,1.0,1.0,1.0,},
};

void process_input_arrays(std::vector<double> x, std::vector<double> y,
                          std::vector<std::vector<double>> f, double *out_x,
                          double *out_y, double *out_f, double *out_fdx,
                          double *out_fdy, double *out_fdxdy);
class BICUBIC_INTERP
{
private:
    // Constants, set with setArray
    double m_dx=1.0, m_dy=1.0; /*Rectilinear weights*/
    double m_minx, m_maxx, m_miny, m_maxy;
    std::vector<double> m_x, m_y; /*Input arrays for interpolations*/

    int m_nx, m_ny;

    // Container for the coefficients.
    std::vector<double> m_a;
    int m_cell_row, m_cell_col;

    void rcin(double x, double y, double &out_cell_x, double &out_cell_y);
    void interpolate(int r, int c);

public:
    std::vector<std::vector<double>> m_f;
    std::vector<std::vector<double>> m_fdx;
    std::vector<std::vector<double>> m_fdy;
    std::vector<std::vector<double>> m_fdxdy;
    BICUBIC_INTERP();
    ~BICUBIC_INTERP(){};

    void prepareContainers(int number_of_omp_threads=1);
    void setArrays(std::vector<double> x, std::vector<double> y,
                   std::vector<std::vector<double>> f);
    void getValues(double x, double y, double &val, double &valdx,
                   double &valdy);
    void debugGetValues(double x, double y, double &val, double &valdx,
                        double &valdy, double &valdxdy);

};
#endif /*BICUBIC_H*/
