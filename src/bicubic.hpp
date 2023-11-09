#ifndef BICUBIC_H
#define BICUBIC_H

#include <vector>

void process_input_arrays(std::vector<double> x, std::vector<double> y,
                          std::vector<std::vector<double>> f, double *out_x,
                          double *out_y, double *out_f, double *out_fdx,
                          double *out_fdy, double *out_fdxdy);
class BICUBIC_INTERP
{
private:
    /// Container for the coefficients. These coefficients are calculated for
    /// each cell.
    std::vector<double> m_a;

    /// To avoid recalculating the m_a coefficients, we remember the last cell
    /// we were in case the new query point lies inside the same cell.
    int m_cell_row, m_cell_col;


public:
    /// The Recompute Coefficients If Needed function. When calling the
    /// getValues or getAllValues function, the index of the cell where the
    /// queried point lies has to be computed. Then this function checks if for
    /// that cell we already have the interpolation coefficients. If not then
    /// the interpolate function is called.
    void rcin(double x, double y, double &out_cell_x, double &out_cell_y);

    /// Actually interpolates the function by generating the interpolant
    /// coefficients stored in m_a. These coefficients are then used when
    /// calling getValues function.
    void interpolate(int r, int c);
    /// Rectilinear weights. Obtained from the input X- and Y-axis arrays and
    /// are needed when obtaining the partial derivatives.
    double m_dx=1.0, m_dy=1.0;

    /// Domain borders.
    double m_minx, m_maxx, m_miny, m_maxy;

    /// Input 1D arrays of the X- and Y-axis. These two arrays form the cell
    /// grid or a rectilinear mesh grid.
    std::vector<double> m_x, m_y;

    /// Size of the X- and Y-axis arrays.
    int m_nx, m_ny;

    /// 2D vectors for storing the discrete values of the function we wish to
    /// interpolate, the first partial derivatives and the mixed partial
    /// derivative.
    std::vector<std::vector<double>> m_f;
    std::vector<std::vector<double>> m_fdx;
    std::vector<std::vector<double>> m_fdy;
    std::vector<std::vector<double>> m_fdxdy;

    BICUBIC_INTERP(){};
    ~BICUBIC_INTERP(){};

    /// Resizes the vectors so that OpenMP threads will access their respective
    /// thread related storages.
    void prepareContainers(int number_of_omp_threads=1);

    /// Function that takes the X-axis and Y-axis 1D arrays and the 2D array
    /// which we wish to interpolate. For the bicubic interpolation the first
    /// partial derivatives are required, as well as the second partial mixed
    /// derivative. This function obtains the derivative values by using higher
    /// finite difference table (combination of central, forward and backward,
    /// depending where are we.).
    void setArrays(std::vector<double> x, std::vector<double> y,
                   std::vector<std::vector<double>> f);

    /// Obtain the function value and its first partial derivatives at point
    /// (x, y). **DOES NOT EXPTRAPOLATE**. If the Query point points outside
    /// the numerical domain defined by the function domain, it is
    /// automatically clipped o the border.
    void getValues(double x, double y, double &val, double &valdx,
                   double &valdy);

    /// Same as getValues, except it also provides the second mixed partial
    /// derivative.
    void getAllValues(double x, double y, double &val, double &valdx,
                      double &valdy, double &valdxdy);

    /// Get the second order derivatives.
    void getSecondDerivatives(double x, double y, double &valdxdx, double &valdydy);

};
#endif /*BICUBIC_H*/
