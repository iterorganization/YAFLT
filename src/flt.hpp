#ifndef FLT_H
#define FLT_H

// RKF45 solver method
#include <rkf45.hpp>

// Bicubic interpolation method
#include <bicubic.hpp>


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
    BICUBIC_INTERP *m_interp_psi;

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

    /// Left border of major radius R value that defines the computation
    /// domain. Along with m_r_max, m_z_min and m_z_max this is used to catch
    /// FLs that might escape the tokamak area. The values are set when the
    /// R and Z points are delivered via setRARR and setZARR methods. In
    /// meters.
    double m_r_min=3.0;
    /// Right border of major radius R value that defines the computation
    /// domain. In meters.
    double m_r_max=9.0;

    /// Bottom border of the Z value that defines the computation domain. In
    /// meters.
    double m_z_min=-6.0;
    /// Top border of the Z value that defines the computation domain. In
    /// meters.
    double m_z_max=6.0;

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
    /// @param[in] option is whether 0, BY_TIME, or 1, BY_LENGTH.
    void setFltOption(int option);

    /// Sets dimension of the (R,Z) grid. This has to be set before adding any
    /// radial or vertical points, psi values or flux points and fpol values as
    /// if raises an error if the dimensions do not match.
    /// @param[in] NDIMR number of points in radial direction
    /// @param[in] NDIMZ number of points in vertical direction.
    void setNDIM(size_t NDIMR, size_t NDIMZ);

    /// Sets the radial points array of size NDIMR.
    /// @param[in] r is the array, which contains NDIMR number of radial points
    void setRARR(std::vector<double> r);
    /// Sets the vertical points array of size NDIMZ
    /// @param[in] r is the array, which contains NDIMZ number of vertical
    ///              points
    void setZARR(std::vector<double> z);
    /// Sets the psi values array of size NDIMR*NDIMZ.
    /// @param[in] psi is a 1D array, long NDIMR*NDIMZ which contains values of
    ///                the magnetic poloidal flux. Row oriented. Meaning, first
    ///                NDIMR points is the first row, etc...
    void setPSI(std::vector<double> psi);
    /// Sets the flux points of size NDIMR. Size subject to change
    /// @param[in] flux is a 1D array of flux values going from plasma center
    ///                 to plasma boundary.
    void setFARR(std::vector<double> flux);
    /// Sets the value of the poloidal current function (Bt R) of size NDIMR.
    /// Size subject to change
    /// @param[in] fpol is a 1D array which contains NDIMR number of values of
    ///                 the poloidal current function.
    void setFPOL(std::vector<double> fpol);
    /// Sets the poloidal current value (Bt R) in vacuum. Units (m T)
    /// @param[in] vacuum_fpol is the F=Bt*R value in vacuum.
    void setVacuumFPOL(double vacuum_fpol){m_vacuum_fpol = vacuum_fpol;};
    double getVacuumFPOL(double vacuum_fpol){return m_vacuum_fpol;};

    /// Sets the plasma shift in (m, m)
    /// @param[in] r_move is the radial shift.
    /// @param[in] z_move is the vertical shift.
    void setShift(double r_move, double z_move);

    /// Sets the desired absolute error for the RKF45 solver
    /// @param[in] abserr is the RKF45 parameter for absolute accuracy.
    void setAbsError(double abserr){m_abserr = abserr;};
    /// Sets the desired relative error for the RKF45 solver
    /// @param[in] relerr is the RKF45 parameter for relative accuracy.
    void setRelError(double relerr){m_relerr = relerr;};
    /// Sets the initial value for running in single thread. When running in
    /// serial only for one thread the containers are prepared and the user does
    /// not need to give the thread id.
    /// The initial value is defined in the Cylindrical coordinate system.
    /// @param[in] r is the radial value of the starting point
    /// @param[in] z is the vertical value of the starting point
    /// @param[in] phi is the toroidal (cylindrical) angle of the starting
    ///                point.
    /// @param[in] omp_thread is the index of the OpenMP, if activated, thread.
    ///                       Defaults to 0 as "serial" execution.
    void setIV(double r, double z, double phi, int omp_thread=0);
    /// Sets the time span where t_end denotes the end for the toroidal angle and
    /// t_step denotes the resolution at which we wish to gather FL points.
    /// @param[in] t_end is the end toroidal angle or the end "time", to which
    ///                  a FL is followed if option=BY_TIME. For each starting
    ///                  point the "time" or angle starts at 0, which means if
    ///                  the t_end is set to 10 degrees then FLs are followed
    ///                  from their starting points for 10 degrees toroidially.
    /// @param[in] t_step is the toroidal angle step at which we wish to gather
    ///                   FLT points, used for FLT or basically tracking.
    void setTimeSpan(double t_end, double t_step);
    /// Sets the maximum connection length in meters to follow a FL.
    /// @param[in] max_conlen is the maximum connection length to which we
    ///                       follow a FL, if the termination is determined by
    ///                       its length. Otherwise ignored if following till a
    ///                       specified toroidal angle.
    void setMaximumConnectionLength(double max_conlen){m_max_conlen = max_conlen;};
    /// Sets the direction of a FL (1 or -1) for running in single thread. When
    /// running in serial only for one thread the containers are prepared and the
    /// user does not need to give the thread id.
    /// @param[in] direction specifies the direction we traverse the FL. Mainly
    ///                      this is used to correctly travel away from the
    ///                      geometry we are tracing.
    /// @param[in] omp_thread is the index of the OpenMP, if activated, thread.
    ///                       Defaults to 0 as "serial" execution.
    void setDirection(int direction, int omp_thread=0);
    /// Calculates the barycenter of a set of triangle points. Orientation is
    /// determined from the order of the ID list of a triangle. :
    ///                      p3 - - p2
    ///                       \     /
    ///                        \   /
    ///                          p1
    /// @param[in] p1 is the first point of the triangle
    /// @param[in] p2 is the second point of the triangle
    /// @param[in] p3 is the third point of the triangle
    void getBaryCenter(double p1[3], double p2[3], double p3[3]);
    /// Sets the self-intersection avoidance length in order to avoid getting
    /// false positive self intersection results. In meters.
    ///
    /// Note that when this is set to a positive value, the RKF45 solver steps
    /// is shortened so that we do not jump out of this length in the first
    /// solver step.
    ///
    /// @param[in] value is the value at which we would like to avoid
    ///                  intersection test when following a FL. This value
    ///                  should be small, in the range of millimeters.
    void setSelfIntersectionAvoidanceLength(double value) {m_self_intersection_avoidance_length = value;};
    ///Prepares the interpolation objects.
    bool prepareInterpolation();

    /// Returns a 1D array of points for a fieldline that is tracked until end
    /// of time (toroidal angle) or the maximum connection length from the
    /// initial starting point set by setIV. The points saved in the vector are
    /// in cylindrical coordinates (R, Z, phi) with units (m, m, rad).
    ///
    /// FLT can be activated with the argument with_flt.
    /// @param[in, out] storage is a 1D array containing points of the FL.
    /// @param[in] with_flt activates or deactivates intersection tests when
    ///                     following FL.
    /// @param[in] omp_thread is the index of the OpenMP, if activated, thread.
    ///                       Defaults to 0 as "serial" execution.
    void getFL(std::vector<double>& storage, bool with_flt=false, int omp_thread=0);

    /// Runs a FLT from a starting initial value set by setIV until the end
    /// of time (toroidal angle) or the maximum connection length. The final
    /// result is a value stored and m_conlens and m_geom_hit_ids.
    /// By default if only one thread is run, the default thread ID is 0. If
    /// openMP is activated, then the used must provide thread ID as argument,
    /// otherwise
    /// @param[in] omp_thread is the index of the OpenMP, if activated, thread.
    ///                       Defaults to 0 as "serial" execution.
    void runFLT(int omp_thread=0);

    /// m_conlens tells us the length of the FL, evaluated in FLT functions.
    /// By default if only one thread is run, aka sequential, this is a vector
    /// of size 1.
    std::vector<double> m_conlens;

    /// m_geom_hit_ids tells us which geometry, an integer registered by embree
    /// when commiting a geometry to it's shadowing structure. Since the
    /// geometry ID's in Embree are non-negative, the negative values tells us
    /// the following:
    /// -1 : Connection length has reached it's end toroidal time or the
    ///      maximum connection length
    /// -2 : The FL has left the computational RZ domain, defined in the
    ///      equilibrium. Since interpolation is used, there is no use to
    ///      follow FLs outside the interpolation area, as extrapolation would
    ///      be required
    /// -3 : The FL wasn't followed at all, i.e., the initialized value
    std::vector<int> m_geom_hit_ids;

    /// Derivative function for the RKF45 object. Since the FL are traced in the
    /// parametric time or toroidal angle, the values stored in the derivative
    /// variable yp yield the ratio between the poloidal components of the
    /// magnetic field and the toroidal components of the magnetic field
    /// (tan of the angle between the components).
    /// The derivative values for the direction in the (R, Z) space
    /// respectfully.
    /// \f[\frac{d s}{d R} = - \frac{d \Psi}{d Z} * \frac{R}{Fpol}\f]
    /// \f[\frac{d s}{d Z} =   \frac{d \Psi}{d R} * \frac{R}{Fpol}\f]
    /// @param[in] t is the parametric time, in this case toroidal angle
    /// @param[in] y contains the current location of the FL on the R,Z map
    /// @param[in,out] yp holds the derivative values for the next step in
    ///                   solving
    /// @param[in] omp_thread is the OpenMP thread.
    void r8_flt(double t, double y[2], double yp[2]);
    /// For a given point in the (R, Z, Phi) space in units of (m, m, rad)
    /// return the poloidal and toroidal component of the magnetic field.
    /// @param[in] r is the radial position of a point
    /// @param[in] z is the vertical position of a point
    /// @param[in] phi is the toroidal position of a point
    /// @param[in, out] out holds the value of the poloidal and toroidal
    ///                     component of the magnetic field
    void getBCyln(double r, double z, std::vector<double> &out);
    /// For a given point in the (R, Z, Phi) space in units of (m, m, rad) return
    /// the magnetic field vector in Cartesian coordinate system. The values
    /// are written in the out variable.
    /// @param[in] r is the radial position of a point
    /// @param[in] z is the vertical position of a point
    /// @param[in, out] out holds the value of the magnetic field vector in the
    ///                     Cartesian coordinate system.
    void getBCart(double r, double z, double phi, std::vector<double> &out);
    /// For a given point in the (R, Z) space in units of (m, m) return the
    /// poloidal current function in units of m T (Bt R).
    /// @returns fpol which is the value of the poloidal current function
    /// @param[in] r is the radial position
    /// @param[in] z is the vertical position
    double getFPol(double r, double z);
    /// For a given point in the flux space in units of (Webb/rad) return the
    /// poloidal current function in units of m T (Bt R).
    /// @returns fpol which is the value of the poloidal current function
    /// @param[in] flux is the poloidal magnetic flux value, or the magnetic
    ///                 surface at which we are interested to obtain the fpol.
    double getFPol(double flux);
    /// Get's the Poloidal current function in vacuum
    /// @returns m_vacuum_fpol is the value of the poloidal current function
    ///                        in the vacuum. It is constant for whole vacuum.
    double getVacuumFPOL(){return m_vacuum_fpol;};
    /// For a given point in the (R, Z) space in units of (m, m) return the
    /// poloidal flux value in Webb/rad
    /// @returns flux which is the poloidal magnetic flux value
    /// @param[in] r is the radial position
    /// @param[in] z is the vertical position
    double getPoloidalFlux(double r, double z);
    /// Sets the pointer to the Embree object responsible for the tracing.
    /// @param[in] accellObj is the Embree accelerated structure used for
    ///                      performing the intersection tests.
    void setEmbreeObj(EmbreeAccell* accellObj);
    /// The following function allocates n of each local thread
    /// values. To be called before running the OpenMP functions in parallel
    /// OpenMP block. When running sequentially user must just call this
    /// function so that the containers are prepared.
    /// @param[in] n tells us how many threads we have. By default, when no
    ///              OpenMP blocks are used this value should be 1. Otherwise
    ///              this value should equal to the number of threads we spawn
    ///              when using OpenMP
    void prepareThreadContainers(int n=1);

    void debug_getValues(double r, double z, double &val, double &valdx, double &valdy, int omp_index=0);
};

#endif //FLT_H