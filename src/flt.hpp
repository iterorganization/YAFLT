#ifndef FLT_H
#define FLT_H

#include <rkf45.hpp>

#include <interpolation.h>
#include <cmath>
#include <accell_embree.hpp>

class RKF45; // Forward declaration

/// Class that solves the FL equation to obtain FL segments and checking if
/// there is any intersection with the shadowing geometry.
class FLT
{
private:


    /// Vector that holds the RKF45 solvers. Reason for that is that the solver
    /// is not thread-safe, therefore for parallel reasons a vector is created
    /// to hold the RKF45 solvers.
    std::vector<RKF45*> m_rkf45_solvers;

    /// ALLIB 2D spline interpolation object m_interp_psi is the interpolation
    /// object on PSI (or flux, units Webb/rad)
    alglib::spline2dinterpolant m_interp_psi;
    /// ALGLIB 1D spline interpolation m_interp_fpol is the interpolation
    /// object for FPOL (toroidal current function, units m T).
    alglib::spline1dinterpolant m_interp_fpol;

    /// Flag for stating if the interpolating objects are prepared
    bool m_prepared=false;

    /// Dimensions of the R, Z space.
    size_t m_NDIMR=-1;
    size_t m_NDIMZ=-1;

    /// Array of R points describing the radial coordinates in unit of meter.\n
    /// Together the R and Z points create the (R, Z) plane in units (m, m) if
    /// one make pairs for every value in one array with the values of the other
    /// array. Or a meshgrid.
    std::vector<double> m_r_points;
    /// Array of Z points, describing the vertical coordinates in unit of
    /// meter.\n
    /// Together the R and Z points create the (R, Z) plane in units (m, m) if
    /// one make pairs for every value in one array with the values of the other
    /// array. Or a meshgrid.
    std::vector<double> m_z_points;

    /// 1D array of poloidal flux values in unit of Webb/rad. The size of this
    /// array must be, m_NDIMR * m_NDIMZ. In other words the rows of the PSI
    /// matrix are concatenated into one array.
    std::vector<double> m_psi_values;

    /// 1D array of poloidal flux values in units of Webb/rad. This defines the
    /// definition space of the FPOL function, which gives us the toroidal
    /// component of the magnetic field times the major radius. NOTE: Not sure
    /// if this is the correct way of obtaining the toroidal component of the
    /// magnetic field outside of the plasma.\n
    /// Note PSI == FLUX in terms of units. Different names are used to
    /// distinguish the role, as one acts as function values and other as
    /// function points.
    std::vector<double> m_flux_points; //FPOL domain

    /// 1D array of poloidal current function values in unit of m T for the
    /// toroidal component interpolator
    std::vector<double> m_fpol_values;

    /// The poloidal current function F = B_t R in the vacuum.
    double m_vacuum_fpol;

    /// FLT option index. Determines the termination condition when following
    /// FLs. By default 0.\n
    /// 0 - BY_TIME - terminate until a set time (toroidal angle)\n
    /// 1 - BY_LENGTH - terminate until a set connection length.
    int m_flt_option = BY_TIME;

    /// Initial values in (R, Z, Phi) in units (m, m, rad).\n
    /// A vector is used regardless of parallel or sequential.
    std::vector<double> m_initial_y;

    // Solver parameters
    // The following parameters are passed to the RKF45 solver.

    /// Absolute error
    double m_abserr=1e-5;
    /// Relative error
    double m_relerr=1e-5;
    /// Time resolution. This basically tells the solver, we want points at
    /// these steps. In radians.
    double m_t_step=1e-4;
    /// End time in radians.
    double m_t_end=0.0;


    /// Self intersection avoidance length. Length at which FLT does not check
    /// for intersections.
    double m_self_intersection_avoidance_length=5e-3;
    /// m_max_conlen states the maximum length to follow a FL. In meters.
    double m_max_conlen=1.0;

    /// Starting direction of a fieldline. By default to determine the
    /// direction of a fieldline, we first take the sign of the FPol or the
    /// toroidal component of the magnetic field. Then, depending on the normal
    /// of our point, by checking the sign of the dot product between the
    /// magnetic field vector and the target normal, we flip the direction so
    /// that we are following the FL on the right side of our element.
    std::vector<int> m_directions;

    /// Output flag of the RKF45 solver. Currently if the data is good enough,
    /// there is no need to process this value.
    int m_flag;


    /// The radial displacement moves in unit of meter. This shift corresponds
    /// to the shift of the plasma. Instead of applying this shift before
    /// interpolation, it is added negatively to the location points when
    /// solving a fieldline.
    double m_r_move = 0.0;
    double m_z_move = 0.0;

    /// Pointer to the EmbreeAccell objects which does the FieldLineTrace. Main
    /// difference between RayTrace is that we use the method RayCast, which
    /// traces finite rays.

    EmbreeAccell* m_embree_obj;

    /// This states if a embree object is set or loaded to the FLT class. Since
    /// the FLT class does not handle the creation of Embree objects then this
    /// is the way to check if m_embree_obj is a valid pointer
    bool m_embree_obj_loaded=false;

public:
    FLT();
    ~FLT(){};

    /// FLT options enumerator. FL can be followed to a certain length or to a
    /// certain toroidal angle
    enum FLT_OPTION {BY_TIME, BY_LENGTH};

    /// Sets the FLT option. If outside the range, default is BY_TIME
    void setFltOption(int option);

    /// Sets dimension of the (R,Z) grid. This has to be set before adding any
    /// radial or vertical points, psi values or flux points and fpol values as
    /// if raises an error if the dimensions do not match.
    void setNDIM(size_t NDIMR, size_t NDIMZ);

    /// Sets the radial points array of size NDIMR.
    void setRARR(std::vector<double> r);
    /// Sets the vertical points array of size NDIMZ
    void setZARR(std::vector<double> z);
    /// Sets the psi values array of size NDIMR*NDIMZ.
    void setPSI(std::vector<double> psi);
    /// Sets the flux points of size NDIMR. Size subject to change
    void setFARR(std::vector<double> flux);
    /// Sets the value of the poloidal current function (Bt R) of size NDIMR.
    /// Size subject to change
    void setFPOL(std::vector<double> fpol);
    /// Sets the poloidal current value (Bt R) in vacuum. Units (m T)
    void setVacuumFPOL(double vacuum_fpol){m_vacuum_fpol = vacuum_fpol;};

    /// Sets the plasma shift in (m, m)
    void setShift(double r_move, double z_move);

    /// Sets the desired absolute error for the RKF45 solver
    void setAbsError(double abserr){m_abserr = abserr;};
    /// Sets the desired relative error for the RKF45 solver
    void setRelError(double relerr){m_relerr = relerr;};
    /// Sets the initial value for running in single thread. When running in
    /// serial only for one thread the containers are prepared and the user does
    /// not need to give the thread id.
    void setIV(double r, double z, double phi, int omp_thread=0);
    /// Sets the time span where t_end denotes the end for the toroidal angle and
    /// t_step denotes the resolution at which we wish to gather FL points.
    void setTimeSpan(double t_end, double t_step);
    /// Sets the maximum connection length in meters to follow a FL.
    void setMaximumConnectionLength(double max_conlen){m_max_conlen = max_conlen;};
    /// Sets the direction of a FL (1 or -1) for running in single thread. When
    /// running in serial only for one thread the containers are prepared and the
    /// user does not need to give the thread id.
    void setDirection(int direction, int omp_thread=0);
    /// Calculates the barycenter of a set of triangle points. Orientation is
    /// determined from the order of the ID list of a triangle. :
    ///                      p3 - - p2
    ///                       \     /
    ///                        \   /
    ///                          p1
    void getBaryCenter(double p1[3], double p2[3], double p3[3]);
    /// Sets the self-intersection avoidance length in order to avoid getting
    /// false positive self intersection results.
    void setSelfIntersectionAvoidanceLength(double value) {m_self_intersection_avoidance_length = value;};
    ///Prepares the interpolation objects.
    bool prepareInterpolation();

    /// Returns a 1D array of points for a fieldline that is tracked until end
    /// of time (toroidal angle) or the maximum connection length from the
    /// initial starting point set by setIV. The points saved in the vector are
    /// in cylindrical coordinates (R, Z, phi) with units (m, m, rad).
    ///
    /// FLT can be activated with the argument with_flt.
    void getFL(std::vector<double>& storage, bool with_flt=false, int omp_thread=0);

    /// Runs a FLT from a starting initial value set by setIV until the end
    /// of time (toroidal angle) or the maximum connection length. The final
    /// result is a value stored in m_hits and m_conlens.
    /// By default if only one thread is run, the default thread ID is 0. If
    /// openMP is activated, then the used must provide thread ID as argument,
    /// otherwise
    void runFLT(int omp_thread=0);

    /// m_hit tells us if the FL, observed in FLT functions reported a hit with
    /// the shadowing geometry.
    /// By default if only one thread is run, aka sequential, this is a vector
    /// of size 1.
    std::vector<bool> m_hits;

    /// m_conlens tells us the length of the FL, evaluated in FLT functions.
    /// By default if only one thread is run, aka sequential, this is a vector
    /// of size 1.
    std::vector<double> m_conlens;

    /// Current position of FLs. Used in case we wish to continue FLT from
    /// current state.
    std::vector<double> m_current_y;

    /// Derivative function for the RKF45 object. Since the FL are traced in the
    /// parametric time or toroidal angle, the values stored in the derivative
    /// variable yp yield the ratio between the poloidal components of the
    /// magnetic field and the toroidal components of the magnetic field
    /// (tan of the angle between the components).
    /// The derivative values for the direction in the (R, Z) space
    /// respectfully.
    /// \f[\frac{d s}{d R} = - \frac{d \Psi}{d Z} * \frac{R}{Fpol}\f]
    /// \f[\frac{d s}{d Z} =   \frac{d \Psi}{d R} * \frac{R}{Fpol}\f]
    void r8_flt(double t, double y[2], double yp[2]);

    /// For a given point in the (R, Z, Phi) space in units of (m, m, rad) return
    /// the magnetic field vector in Cartesian coordinate system. The values
    /// are written in the out variable.
    void getBCart(double r, double z, double phi, std::vector<double> &out);
    /// For a given point in the (R, Z) space in units of (m, m) return the
    /// poloidal current function in units of m T (Bt R).
    double getFPol(double r, double z);
    /// For a given point in the flux space in units of (Webb/rad) return the
    /// poloidal current function in units of m T (Bt R).
    double getFPol(double flux);
    /// Get's the Poloidal current function in vacuum
    double getVacuumFPOL(){return m_vacuum_fpol;};
    /// For a given point in the (R, Z) space in units of (m, m) return the
    /// poloidal flux value in Webb/rad
    double getPoloidalFlux(double r, double z);
    /// Sets the pointer to the Embree object responsible for the tracing.
    void setEmbreeObj(EmbreeAccell* accellObj);
    /// The following function allocates n of each local thread
    /// values. To be called before running the OpenMP functions in parallel
    /// OpenMP block. When running sequentially user must just call this
    /// function so that the containers are prepared.
    void prepareThreadContainers(int n=1);
};

#endif //FLT_H