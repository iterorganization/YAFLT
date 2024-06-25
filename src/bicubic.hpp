#ifndef BICUBIC_H
#define BICUBIC_H

#include <vector>

// BI_DATA struct is used for providing context when obtaining interpolation
// values. Also more OpenMP friendly than plain shared class attributes.
struct alignas(64) BI_DATA
{
    double r=0;
    double z=0;
    double val=0;
    double valdx=0;
    double valdy=0;
    double dx=0.0;
    double dy=0.0;
    double minx=0.0;
    double miny=0.0;
    double maxx=0.0;
    double maxy=0.0;
    int nx=0;
    int ny=0;
};

class BICUBIC_INTERP
{

public:
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
    /// 1D vector that holds the constants. The Size is (Nx*Ny*16);
    std::vector<double> m_a;

    BICUBIC_INTERP(){};
    ~BICUBIC_INTERP(){};

    /// Function that takes the X-axis and Y-axis 1D arrays and the 2D array
    /// which we wish to interpolate. For the bicubic interpolation the first
    /// partial derivatives are required, as well as the second partial mixed
    /// derivative. This function obtains the derivative values by using higher
    /// finite difference table (combination of central, forward and backward,
    /// depending where are we.).
    void setArrays(std::vector<double> x, std::vector<double> y,
                   std::vector<std::vector<double>> f);

    void populateContext(BI_DATA *context){
        context->dx = m_dx;
        context->dy = m_dy;
        context->minx = m_minx;
        context->miny = m_miny;
        context->nx = m_nx;
        context->ny = m_ny;
        context->maxx = m_maxx;
        context->maxy = m_maxy;
    }

    /// Obtain the function value and its first partial derivatives at point
    /// (x, y). **DOES NOT EXPTRAPOLATE**. If the Query point points outside
    /// the numerical domain defined by the function domain, it is
    /// automatically clipped to the border.
    void getValues(BI_DATA *context);
    void getSecondDerivativeValues(BI_DATA *context);

};
#endif /*BICUBIC_H*/
