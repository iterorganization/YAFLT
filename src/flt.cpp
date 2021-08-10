#include <flt.hpp>
#include <accell_embree.hpp>

// ALGLIB
#include <ap.h>
#include <interpolation.h>

// RKF45 includes
#include <rkf45.hpp>

// STD includes
#include <cassert>     // For assert
#include <vector>
#include <cmath>

FLT::FLT(){
    return;
}

void FLT::prepareThreadContainers(int num_threads){
    // Prepares thread variables. Needs to be called before trying to run
    // parallel, otherwise a segfault is excepted.
    //
    // In essence create a copy of a local variable for each thread.
    //
    // Allocates the following variables:
    //     RKF45 solver - This class has some state variables needed for
    //         successful solving, therefore each thread needs to have one
    //     bool intersection - This holds the information on whether the FLT
    //         detected a hit of the FL.
    //     double conlen - This holds the information on the connection length of
    //         the FLT.
    //     vector<double> IV - initial value vector of 3 doubles.
    //     double timeStop - End of toroidal (angle) time
    //     double timeStep - The time resolution at which we wish to obtain
    //         results
    //
    // After clearing the vectors, populate with a default value or object.
    //
    // Subsequently call the Embree allocation to do the same as structs are
    // populated during RTC and cannot be shared between threads.

    m_rkf45_solvers.clear();
    m_hits.clear();
    m_conlens.clear();
    m_initial_y.clear();
    m_current_y.clear();
    m_directions.clear();

    for(int i=0; i < num_threads; i++){
        m_rkf45_solvers.push_back(new RKF45());
        m_hits.push_back(false);
        m_conlens.push_back(0.0);
        m_initial_y.push_back(0.0);
        m_initial_y.push_back(0.0);
        m_initial_y.push_back(0.0);
        m_current_y.push_back(0.0);
        m_current_y.push_back(0.0);
        m_current_y.push_back(0.0);
        m_directions.push_back(1);
    }

    if (m_embree_obj_loaded){
        // Prepare embree RayHits
        m_embree_obj->prepareThreadContainers(num_threads);
    }
}

void FLT::setNDIM(size_t NDIMR, size_t NDIMZ){
    m_NDIMR = NDIMR; // Number of rows
    m_NDIMZ = NDIMZ; // Number of columns

    // /* Resize vectors */
    // m_r_points.resize(NDIMR); /*R points*/
    // m_z_points.resize(NDIMZ); /*Z points*/
    // m_psi_values.resize(NDIMR, std::vector<double>(NDIMZ)); /*Psi values*/

    // m_flux_points.resize(NDIMR); /*Psi domain from center to boundary*/
    // m_fpol_values.resize(NDIMR); /*FPOL value*/
}

void FLT::setRARR(std::vector<double> r){
    // Sets the radial points in meters with size NDIMR.
    m_r_points = r;

    m_r_min = r[0];
    m_r_max = r[r.size() - 1];
    assert(("Number of elements in m_r_points is not the same as m_NDIMR", \
            m_r_points.size() == m_NDIMR));
}
void FLT::setZARR(std::vector<double> z){
    // Sets the vertical points in meters with size NDIMZ.
    m_z_points = z;
    m_z_min = z[0];
    m_z_max = z[z.size() - 1];
    assert(("Number of elements in m_z_points is not the same as m_NDIMZ", \
            m_z_points.size() == m_NDIMZ));
}

void FLT::setPSI(std::vector<double>psi){
    // Sets the poloidal flux values in Web/rad with size NDIMR*NDIMZ.
    m_psi_values = psi;
    assert(("Number of elements in m_psi_values is not the same as m_NDIMR * m_NDIMZ", \
            m_psi_values.size() == m_NDIMR * m_NDIMZ));
}

void FLT::setFARR(std::vector<double> flux){
    // Sets the poloidal flux points in Web/rad with size NDIMR.
    m_flux_points = flux;
    assert(("Number of elements in m_flux_points is not the same as m_NDIMR", \
            m_flux_points.size() == m_NDIMR));
}

void FLT::setFPOL(std::vector<double> fpol){
    // Sets the poloidal current values in m T with size NDIMR.
    m_fpol_values = fpol;
    assert(("Number of elements in m_fpol_values is not the same as m_NDIMR", \
            m_fpol_values.size() == m_NDIMR));
}


void FLT::setFltOption(int option){
    if (option >= BY_TIME && option <= BY_LENGTH){
        m_flt_option = option;
    }
    else {
        m_flt_option = BY_TIME;
    }
}

void FLT::setIV(double r, double z, double phi, int omp_thread){
    // Sets the initial value in (R, Z, Phi) space in units (m,m,rad) for
    // a given OpenMP thread.
    //
    // When running sequentially, the default value for omp_thread is 0.
    m_initial_y[3*omp_thread] = r;
    m_initial_y[3*omp_thread + 1] = z;
    m_initial_y[3*omp_thread + 2] = phi;
}

void FLT::setShift(double rmove, double zmove){
    // Sets the plasma shifts in units of m, m. The shift is applied when
    // calling the interpolation method and is applied negatively as this is
    // intended for the plasma shift.
    m_r_move = rmove;
    m_z_move = zmove;
}

void FLT::setTimeSpan(double t_end, double t_step){
    // Sets the time end and time step in units rad, rad (actually toroidal
    // angle). The time end corresponds to how far toroidally we wish to follow
    // the FL and time step is the resolution (i.e., at what steps) we wish to
    // gather points. The RKF45 solver while solving returns values at these
    // steps but in between it changes step sizes according to RKF45 scheme.
    m_t_step=t_step;
    m_t_end=t_end;
    // std::cout << m_t_end << " " << m_t_step << std::endl;
}

void FLT::setDirection(int direction, int omp_thread){
    // Direction of a FL to follow from the starting point set by setIV. In
    // principle we always wish to start away from the front surface of a PFC
    // component but for certain scans, such as the connection length at
    // midplane, we wish to follow in both directions.
    //
    // When running sequentially, the default value for omp_thread is 0.
    //
    // VALUE SHOULD BE ONLY 1 or -1*/
    m_directions[omp_thread] = direction;
}

bool FLT::prepareInterpolation(){
    // Allocates and creates the interpolation objects. As well as creates the
    // RKF45 object.
    alglib::real_1d_array R, Z, Psi, Flux, FPol;

    R.setcontent(m_r_points.size(), &(m_r_points[0]));
    Z.setcontent(m_z_points.size(), &(m_z_points[0]));
    Psi.setcontent(m_psi_values.size(), &(m_psi_values[0]));
    Flux.setcontent(m_flux_points.size(), &(m_flux_points[0]));
    FPol.setcontent(m_fpol_values.size(), &(m_fpol_values[0]));

    spline1dbuildcubic(Flux, FPol, m_interp_fpol);
    spline2dbuildbicubicv(R, m_NDIMR, Z, m_NDIMZ, Psi, 1, m_interp_psi);

    m_prepared = true;
    return true;
}

void FLT::r8_flt(double t, double y[2], double yp[2]){
    // Derivative function for the RKF45. It evaluates the ratio between the
    // (R, Z) magnetic components and the toroidal magnetic component in the
    // cylindrical coordinate system.
    double derivFluxdX, derivFluxdY, factor;

    spline2ddiff(m_interp_psi, y[0] - m_r_move, y[1] - m_z_move, factor,
                 derivFluxdX, derivFluxdY, factor);
    factor = y[0] / m_vacuum_fpol;
    yp[0] = - derivFluxdY * factor;
    yp[1] =   derivFluxdX * factor;
}

void FLT::getFL(std::vector<double>& storage, bool with_flt, int omp_thread){
    // Follows and stores FL points into the vector container.
    //
    // In this case the FL is followed with FLT activated.
    //
    // The state_index tells us what is the termination condition.
    //  0 - by time (toroidal angle) - default
    //  1 - by length (connection length)

    int option_index;


    // state_data contains information on current time and length
    double state_data[2];
    // stop_data contains information on termination options, such as target
    //  toroidal angle and target connection length
    double stop_data[2];

    // relerr - relative error
    // abserr - absolute error
    // t - current time (toroidal angle). Can be negative or positive
    // t_step - desired toroidal angle resolution
    // t_offset - Position of the starting point.
    // new_time - next new time for RKF45 to solve
    // y[2] - R, Z values in Cylindrical coordinate system
    // yp[2] - Derivative of R, Z values in Cylindrical coordinate system
    double relerr, abserr, t, t_step, t_offset, new_time, y[2], yp[2];

    /*
        flag - RKF45 flag value
        direction - toroidal direciton, 1 Counter-CW, -1 CW
    */
    int flag, direction;

    // pre_t_step - Initial following of FL for avoiding self-intersection
    double pre_t_step;

    option_index = m_flt_option;


    relerr = m_relerr;
    abserr = m_abserr;

    direction = m_directions[omp_thread];

    t_step = m_t_step * direction;

    stop_data[BY_TIME] = m_t_end;
    stop_data[BY_LENGTH] = m_max_conlen;
    // Initial values for state_data
    state_data[BY_TIME] = 0.0;
    state_data[BY_LENGTH] = 0.0;

    // Initial condition;
    y[0] = m_initial_y[3*omp_thread];
    y[1] = m_initial_y[3*omp_thread + 1];
    t_offset = m_initial_y[3*omp_thread + 2];

    t = 0;

    storage.push_back(y[0]);
    storage.push_back(y[1]);
    storage.push_back(t_offset);

    // FLT values
    double ox, oy, oz, x1, y1, z1, norm, dx, dy, dz, tnear, tfar;
    bool intersect=false;

    // Flag=1 means solver will run until end time. If Flag=-1 then solver
    // solves in singular steps where the user can check values.
    flag = 1; // First time set the flag to 1

    if (with_flt){
        // For avoiding self-intersection test. Follow the FL for a small length
        // and then start from there. The value should be around 1 cm but not too
        // much.

        // In order to determine what should be the first step we simply take the
        // major radius of the element and make sure that there are at least 2
        // steps between the origin point and the self avoiding intersection length
        pre_t_step = std::asin(0.1 * m_self_intersection_avoidance_length / y[0]) * direction;

        while (state_data[BY_LENGTH] < m_self_intersection_avoidance_length){
            ox = y[0] * cos(t_offset + t);
            oy = y[0] * sin(t_offset + t);
            oz = y[1];
            new_time = t + pre_t_step;
            while (0 < (new_time - t) * direction){
                flag = m_rkf45_solvers[omp_thread]->r8_rkf45(this, y, yp, &t,
                                                            new_time, &relerr,
                                                            abserr, flag);
                flag = 2; // Until we find an actual problem, ignore messages from
                          // rkf45 as they do not indicate a problem*/
            }
            x1 = y[0] * cos(t_offset + t);
            y1 = y[0] * sin(t_offset + t);
            z1 = y[1];

            // Calculate length of field-line segment for the current step
            // calculate direction of field-line segment for the current step
            dx = (x1-ox);
            dy = (y1-oy);
            dz = (z1-oz);
            norm = sqrt(dx * dx + dy * dy + dz * dz);
            dx = dx / norm;
            dy = dy / norm;
            dz = dz / norm;
            state_data[BY_LENGTH] = state_data[BY_LENGTH] + norm;
            state_data[BY_TIME] = t * direction;
        }
    }

    while (state_data[option_index] < stop_data[option_index] && !intersect){
        ox = y[0] * cos(t_offset + t);
        oy = y[0] * sin(t_offset + t);
        oz = y[1];
        new_time = t + t_step;

        while (0 < (new_time - t) * direction){
            flag = m_rkf45_solvers[omp_thread]->r8_rkf45(this, y, yp, &t,
                                                        new_time, &relerr,
                                                        abserr, flag);
            flag = 2; // Until we find an actual problem, ignore messages from
                      // rkf45 as they do not indicate a problem
        }

        // Check for intersection
        // calculate end of field-line segment for the current step
        x1 = y[0] * cos(t_offset + t);
        y1 = y[0] * sin(t_offset + t);
        z1 = y[1];

        // Calculate length of field-line segment for the current step
        // calculate direction of field-line segment for the current step
        dx = (x1-ox);
        dy = (y1-oy);
        dz = (z1-oz);
        norm = sqrt(dx * dx + dy * dy + dz * dz);
        dx = dx / norm;
        dy = dy / norm;
        dz = dz / norm;

        tnear = 0;
        tfar = norm;

        state_data[BY_LENGTH] = state_data[BY_LENGTH] + norm;
        state_data[BY_TIME] = t * direction;
        // std::cout << y[0] << " " << y[1] << " " << t<< std::endl;
        storage.push_back(y[0]);
        storage.push_back(y[1]);
        storage.push_back(t_offset + t);

        if (with_flt){
           // Check intersection of field-line segment for the current step
            m_embree_obj->castRay(ox, oy, oz, dx, dy, dz, tnear, tfar,
                                      omp_thread);
            intersect = m_embree_obj->checkIfHit(omp_thread);
            // Check if point is outside the computation domain
            if (y[0] < m_r_min || y[0] > m_r_max){
                // For now treat this as an intersect.
                intersect = true;
                continue;
            }
            if (y[1] < m_z_min || y[1] > m_z_max){
                // For now treat this as an intersect.
                intersect = true;
                continue;
            }
        }

    }

}

void FLT::runFLT(int omp_thread){

    // This functions follows a fieldline from the set initial start and uses
    // the provided m_embree_obj embree structure which checks if the FL
    // hits anytging in the shadowing geometry.
    //
    //
    // The omp_thread is the ID of the thread. When running sequentially, the
    // default value is 0 but nevertheless to avoid duplicating too much code the
    // same function is used when calling sequentially and running in a OMP
    // parallel block.
    //
    //
    // The state_index tells us what is the termination condition.
    //   0 - by time (toroidal angle) - default
    //   1 - by length (connection length)
    int option_index;


    // state_data contains information on current time and length
    double state_data[2];
    // stop_data contains information on termination options, such as target
    // toroidal angle and target connection length
    double stop_data[2];

    //  relerr - relative error
    //  abserr - absolute error
    //  t - current time (toroidal angle). Can be negative or positive.
    //      Starts at 0
    //  t_step - desired toroidal angle resolution
    //  t_offset - Position of the starting point.
    //  new_time - next new time for RKF45 to solve
    //  y[2] - R, Z values in Cylindrical coordinate system
    //  yp[2] - Derivative of R, Z values in Cylindrical coordinate system
    double relerr, abserr, t, t_step, t_offset, new_time, y[2], yp[2];

    //    flag - RKF45 flag value
    //    direction - toroidal direciton, 1 Counter-CW, -1 CW
    int flag, direction;

    // pre_t_step - Initial following of FL for avoiding self-intersection
    double pre_t_step;

    option_index = m_flt_option;


    relerr = m_relerr;
    abserr = m_abserr;

    direction = m_directions[omp_thread];

    t_step = m_t_step * direction;

    stop_data[BY_TIME] = m_t_end;
    stop_data[BY_LENGTH] = m_max_conlen;
    // Initial values for state_data
    state_data[BY_TIME] = 0.0;
    state_data[BY_LENGTH] = 0.0;

    // Initial condition;
    y[0] = m_initial_y[3*omp_thread];
    y[1] = m_initial_y[3*omp_thread + 1];
    t_offset = m_initial_y[3*omp_thread + 2];

    t = 0;

    // FLT values
    double ox, oy, oz, x1, y1, z1, norm, dx, dy, dz, tnear, tfar;
    bool intersect=false;

    // Flag=1 means solver will run until end time. If Flag=-1 then solver
    // solves in singular steps where the user can check values.
    flag = 1; // First time set the flag to 1

    // For avoiding self-intersection test. Follow the FL for a small length
    // and then start from there. The value should be around 1 cm but not too
    // much.

    // In order to determine what should be the first step we simply take the
    // major radius of the element and make sure that there are at least 2
    // steps between the origin point and the self avoiding intersection length
    pre_t_step = std::asin(0.1 * m_self_intersection_avoidance_length / y[0]) * direction;
    while (state_data[BY_LENGTH] < m_self_intersection_avoidance_length){
        ox = y[0] * cos(t_offset + t);
        oy = y[0] * sin(t_offset + t);
        oz = y[1];
        new_time = t + pre_t_step;
        while (0 < (new_time - t) * direction){
            flag = m_rkf45_solvers[omp_thread]->r8_rkf45(this, y, yp, &t,
                                                        new_time, &relerr,
                                                        abserr, flag);
            flag = 2; // Until we find an actual problem, ignore messages from
                      // rkf45 as they do not indicate a problem*/
        }

        x1 = y[0] * cos(t_offset + t);
        y1 = y[0] * sin(t_offset + t);
        z1 = y[1];

        // Calculate length of field-line segment for the current step
        // calculate direction of field-line segment for the current step
        dx = (x1-ox);
        dy = (y1-oy);
        dz = (z1-oz);
        norm = sqrt(dx * dx + dy * dy + dz * dz);

        state_data[BY_LENGTH] = state_data[BY_LENGTH] + norm;
        state_data[BY_TIME] = t * direction;
    }

    while (state_data[option_index] < stop_data[option_index] && !intersect){
        ox = y[0] * cos(t_offset + t);
        oy = y[0] * sin(t_offset + t);
        oz = y[1];
        new_time = t + t_step;
        while (0 < (new_time - t) * direction){
            flag = m_rkf45_solvers[omp_thread]->r8_rkf45(this, y, yp, &t,
                                                        new_time, &relerr,
                                                        abserr, flag);
            flag = 2; // Until we find an actual problem, ignore messages from
                      // rkf45 as they do not indicate a problem*/
        }
        // Check for intersection
        // calculate end of field-line segment for the current step
        x1 = y[0] * cos(t_offset + t);
        y1 = y[0] * sin(t_offset + t);
        z1 = y[1];

        // Calculate length of field-line segment for the current step
        // calculate direction of field-line segment for the current step
        dx = (x1-ox);
        dy = (y1-oy);
        dz = (z1-oz);
        norm = sqrt(dx * dx + dy * dy + dz * dz);
        dx = dx / norm;
        dy = dy / norm;
        dz = dz / norm;

        tnear = 0;
        tfar = norm;

        state_data[BY_LENGTH] = state_data[BY_LENGTH] + norm;
        state_data[BY_TIME] = t * direction;
        // Check intersection of field-line segment for the current step
        m_embree_obj->castRay(ox, oy, oz, dx, dy, dz, tnear, tfar,
                                      omp_thread);
        intersect = m_embree_obj->checkIfHit(omp_thread);

        // Check if point is outside the computation domain
        if (y[0] < m_r_min || y[0] > m_r_max){
            // For now treat this as an intersect.
            intersect = true;
            continue;
        }
        if (y[1] < m_z_min || y[1] > m_z_max){
            // For now treat this as an intersect.
            intersect = true;
            continue;
        }
    }
    m_hits[omp_thread] = intersect;
    m_conlens[omp_thread] = state_data[BY_LENGTH];
    // Save current position of fieldline if we wish to continue running FLT
    m_current_y[3 * omp_thread] = y[0];
    m_current_y[3 * omp_thread + 1] = y[1];
    m_current_y[3 * omp_thread + 2] = t_offset + t;
}

void FLT::getBCart(double r, double z, double phi, std::vector<double> &out){
    // Calculate the magnetic vector at point P(r, z, phi).
    //
    // This function transforms the magnetic vector from cylindrical to
    // Cartesian system.
    //
    // Transformation matrix from Cylindrical to Cartesian:
    //
    //  cos Phi | -sin Phi | 0
    // ---------+----------+---
    //  sin Phi |  cos Phi | 0
    // ---------+----------+---
    //    0     |    0     | 1
    //
    // Phi is taken from the physical point where we wish to calculate the BCart.
    std::vector<double> magVec;
    double derivFluxdR, derivFluxdZ, bphi, br, bz;

    spline2ddiff(m_interp_psi, r - m_r_move, z - m_z_move, bphi, derivFluxdR,
                 derivFluxdZ, bphi); // Bphi here is a dummy value


    br = -derivFluxdZ / r;
    bphi = m_vacuum_fpol / r;
    bz = derivFluxdR / r;

    out[0] = br * cos(phi) - sin(phi) * bphi;
    out[1] = br * sin(phi) + cos(phi) * bphi;
    out[2] =bz;
    // magVec.push_back(br * cos(phi) - sin(phi) * bphi);
    // magVec.push_back(br * sin(phi) + cos(phi) * bphi);
    // magVec.push_back(bz);
    // return magVec;
}

double FLT::getFPol(double r, double z){
    // Returns value of the poloidal flux current.
    //
    // In order to get the toroidal component one has to divide by major radius,
    // or essentially in this case with R.
    //
    // Bt = FPol / r
    //
    // Used more or less to determine the starting direction of FLT.
    double fpol;
    fpol = spline2dcalc(m_interp_psi, r - m_r_move, z - m_z_move);
    fpol = spline1dcalc(m_interp_fpol, fpol);
    return fpol;
}

double FLT::getFPol(double flux){
    // Returns value of the poloidal flux current.
    //
    // In order to get the toroidal component one has to divide by major radius.
    //
    // Bt = FPol / r
    //
    // Used more or less to determine the starting direction of FLT.
    double fpol;
    fpol = spline1dcalc(m_interp_fpol, flux);
    return fpol ;
}


double FLT::getPoloidalFlux(double r, double z){
    // Returns the value of the poloidal flux at point (r,z).
    double flux=spline2dcalc(m_interp_psi, r - m_r_move, z - m_z_move);
    return flux;
}

void FLT::setEmbreeObj(EmbreeAccell* accellObj){
    // Set the RayTrace pointer to the class
    m_embree_obj = accellObj;
    m_embree_obj_loaded = true;
}