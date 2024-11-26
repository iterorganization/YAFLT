#include <flt.hpp>

// STD includes
#include <vector>

#include <stdio.h>

#include <omp.h>

FLT::FLT(){
    m_interp_psi = new BICUBIC_INTERP();
    return;
}

FLT::~FLT(){
    delete m_interp_psi;
    return;
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
}
void FLT::setZARR(std::vector<double> z){
    // Sets the vertical points in meters with size NDIMZ.
    m_z_points = z;
    m_z_min = z[0];
    m_z_max = z[z.size() - 1];
}

void FLT::setPSI(std::vector<double>psi){
    // Sets the poloidal flux values in Web/rad with size NDIMR*NDIMZ.
    m_psi_values = psi;
}

void FLT::setFARR(std::vector<double> flux){
    // Sets the poloidal flux points in Web/rad with size NDIMR.
    m_flux_points = flux;
}

void FLT::setFPOL(std::vector<double> fpol){
    // Sets the poloidal current values in m T with size NDIMR.
    m_fpol_values = fpol;
}


void FLT::setShift(double rmove, double zmove){
    // Sets the plasma shifts in units of m, m. The shift is applied when
    // calling the interpolation method and is applied negatively as this is
    // intended for the plasma shift.
    m_r_move = rmove;
    m_z_move = zmove;
}

bool FLT::prepareInterpolation(){
    // Allocates and creates the interpolation objects. As well as creates the
    // RKF45 object.

    int n_rows, n_cols;
    std::vector<double> buffer;
    n_rows = m_z_points.size();
    n_cols = m_r_points.size();
    m_reshaped_psi.clear();
    for (int i=0; i<n_rows; i++){
        buffer.clear();
        for (int j=0; j<n_cols; j++){
            buffer.push_back(m_psi_values[i * n_cols + j]);
        }
        // m_reshaped_psi.insert(m_reshaped_psi.begin(), buffer);
        m_reshaped_psi.push_back(buffer);
    }

    m_interp_psi->setArrays(m_r_points, m_z_points, m_reshaped_psi);
    return true;
}

void FLT::getBCyln(double r, double z, std::vector<double> &out){
    // Get the poloidal and toroidal component of the magnetic vector in
    // Cylindrical coordinate system at point P(r, z, phi).

    // Poloidal component is calculated as the sqrt of sum of squares of the
    // derivatives of the Flux in (R, Z) space.

    BI_DATA *ctx = new BI_DATA();
    m_interp_psi->populateContext(ctx);

    ctx->r = r - m_r_move;
    ctx->z = z - m_z_move;
    m_interp_psi->getValues(ctx);

    out[0] = sqrt(ctx->valdx * ctx->valdx + ctx->valdy * ctx->valdy) / r;
    out[1] = m_vacuum_fpol / r; // Bphi
    delete ctx;
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
    double bphi, br, bz;

    BI_DATA *ctx = new BI_DATA();
    m_interp_psi->populateContext(ctx);

    ctx->r = r - m_r_move;
    ctx->z = z - m_z_move;
    m_interp_psi->getValues(ctx);

    br = -ctx->valdy / r;
    bphi = m_vacuum_fpol / r;
    bz = ctx->valdx / r;

    out[0] = br * cos(phi) - sin(phi) * bphi;
    out[1] = br * sin(phi) + cos(phi) * bphi;
    out[2] = bz;
    delete ctx;
    // magVec.push_back(br * cos(phi) - sin(phi) * bphi);
    // magVec.push_back(br * sin(phi) + cos(phi) * bphi);
    // magVec.push_back(bz);
    // return magVec;
}

double FLT::getPoloidalFlux(double r, double z){
    // Returns the value of the poloidal flux at point (r,z).
    double out;
    BI_DATA *ctx = new BI_DATA();
    m_interp_psi->populateContext(ctx);

    ctx->r = r - m_r_move;
    ctx->z = z - m_z_move;
    m_interp_psi->getValues(ctx);
    out = ctx->val;
    delete ctx;
    return out;
}

void FLT::setEmbreeObj(EmbreeAccell* accellObj){
    // Set the RayTrace pointer to the class
    m_embree_obj = accellObj;
    m_embree_obj_loaded = true;
}


void FLT::setPoints(std::vector<double> origin_points){
    /// Set's the points and resizes the output arrays.
    m_origin_points=origin_points;
    const int n_points = (int) m_origin_points.size() / 3.0;
    // Clear is not necessarily needed since all elements get overwritten
    //
    // Clear the output arrays
    m_out_fieldline_lengths.clear();
    m_out_geom_hit_ids.clear();
    m_out_prim_hit_ids.clear();

    // Resize them
    m_out_fieldline_lengths.resize(n_points);
    m_out_geom_hit_ids.resize(n_points);
    m_out_prim_hit_ids.resize(n_points);
}

void FLT::runFLT(){
    // This function follows a fieldline from the set initial start and uses
    // the provided m_embree_obj embree structure which checks if the FL
    // hits anything in the shadowing geometry.

    int n_points = m_origin_points.size();
    if (n_points == 0){
        return;
    }
    n_points = n_points / 3;

    // double n_points_check = n_points / 3.0;
    // if (n_points != n_points_check){
    //     // Raise error
    //     return;
    // }

#ifndef NDEBUG
    printf("Starting up on %d cores\n", m_number_of_threads);
    printf("In total %d points\n", n_points);
#endif

    #pragma omp parallel num_threads(m_number_of_threads)
    {
        // Control options and values set by the user!!
        // Termination length of the field line. Meaning, tracing the field line
        // to this maximum length, unless other conditions terminate the trace
        // earlier.
        const double terminate_length = m_max_fieldline_length;
        // const int direction = m_directions[0];
        // Length where we don't check for RT intersections.
        const double self_intersection_avoidance_length = m_self_intersection_avoidance_length;
        // pre_t_step - Initial following of FL for avoiding self-intersection
        // const double pre_t_step = std::asin(0.1 * m_self_intersection_avoidance_length / y[0]) * direction;
        // t_step is the desired step we wish to make
        const double t_step = m_t_step;
        // Relative and absolute error
        const double relerr = m_relerr;
        const double abserr = m_abserr;
        // Direction - toroidal direction, 1 Counter-CW, -1 CW


        // Length of the followed field line so far.
        double total_length;
        bool intersect, left_computational_domain;
        int direction;
        double pre_t_step;

        // *** RKF4(5) variables.
        // h is the actual step performed during the RKF4(5) steps.
        double h;

        // Step increment arrays
        double k2[2];
        double k3[2];
        double k4[2];
        double k5[2];
        double k6[2];
        // Solution array
        double s[2];
        // Error values
        const double error_scale = 2.0/relerr;
        const double error_absolute = error_scale*abserr;
        double error_et;
        double error_ete;
        double error_eot;
        double error_tolerance;

        // Embree constructs
        RTCRayHit rayHit = RTCRayHit();
        RTCIntersectContext rayContext = RTCIntersectContext();
        rtcInitIntersectContext(&rayContext);
        // Ray tracing variables that are put into the RTCRayHit struct.
        double ox, oy, oz, x1, y1, z1, norm, dx, dy, dz;
        // y[2] and yp[2] correspond to (R, Z) and value of the derivative of the
        // function we wish to solve.
        double y[2], yp[2];
        // t equals to the parameteric time or toroidal angle
        // dt equals to the
        double t;
        double dt, h_scale;

        // Interpolation struct
        BI_DATA *interp_context = new BI_DATA();
        m_interp_psi->populateContext(interp_context);

        #pragma omp for
        for(int i=0; i < n_points; i++){
            // printf("%d\n", i);
            // Initial values;
            y[0] = m_origin_points[3*i];
            y[1] = m_origin_points[3*i+1];
            flt_pde(y, yp, interp_context);

            t = m_origin_points[3*i+2];
            direction = m_directions[i];

            if (self_intersection_avoidance_length > 0.0){
                pre_t_step = std::asin(0.1 * self_intersection_avoidance_length / y[0]) * direction;
                dt = pre_t_step;
            }
            else {
                dt = t_step * direction;
            }

            // Reset field line length.
            total_length = 0.0;

            // Termination flags for the main while loop
            intersect=false;
            left_computational_domain=false;

            // Calculate the (X, Y, Z) starting point of a FL segment
            ox = y[0] * std::cos(t);
            oy = y[0] * std::sin(t);
            oz = y[1];

            // printf("r=%f z=%f t0=%f d=%d dt=%f pre_t_step=%f self_intersection_avoidance_length=%f\n", y[0], y[1], t, direction, dt, pre_t_step, m_self_intersection_avoidance_length);

            while (total_length <= terminate_length && !intersect){
                // Obtain the derivatives
                // flt_pde(t, y, yp);

                // We do not allow limiting the step size. If the minimum step
                // size is set to 1e-6, or if this is our limiting lower value
                // we would have to leave the method running until time
                // 173215370, which corresponds to more than 27 million
                // revolutions around the tokamak. At least if we see how
                // the lower limit dt is calculated.
                //m_h = max(h, 26.0 * eps * max(abs(time), abs(dt)));


                // Single RKF45 step
                // if (std::abs(pre_t_step) <= 26.0 * EPSILON * std::abs(t)){
                //     t = t + pre_t_step;
                //     y[0] = y[0] + dt * yp[0];
                //     y[1] = y[1] + dt * yp[1];
                // }

                h = dt;

                // Start the integration
                while (true){
                    // Advance an approximate solution over one step of length
                    // h. If there is no need to lower the step size then the
                    // solution in variable s is accepted.
                    fehl_step(y, t, h, yp, k2, k3, k4, k5, k6, s, interp_context);
                    // Compute and test allowable tolerances versus local error
                    // estimates and remove scaling of tolerances. The relative
                    // error is measured with respect to the average of the
                    // magnitudes of the solution at the beginning and end of
                    // the step.

                    // FORMULA 2 Table III in Fehlberg for determining
                    // truncation error.

                    // Get average of the function magnitudes at t and t + h
                    error_et = std::abs(y[0]) + std::abs(k2[0]) + error_absolute;
                    // We assume that the errors set are positive!

                    // Truncation error
                    error_ete = std::abs(-2090.0*yp[0] + 21970.0 * k4[0] -
                                         15048.0 * k5[0] + 22528.0 * k3[0] -
                                         27360.0 * k6[0]);

                    // Get the relative error
                    error_eot = error_ete/error_et;

                    // Do the same at second dimension
                    error_et = std::abs(y[1]) + std::abs(k2[1]) + error_absolute;
                    // We assume that the errors set are positive!

                    // Truncation error
                    error_ete = std::abs(-2090.0*yp[1] + 21970.0 * k4[1] -
                                         15048.0 * k5[1] + 22528.0 * k3[1] -
                                         27360.0 * k6[1]);

                    // Take the max relative error of the two
                    error_eot = std::max(error_eot, error_ete/error_et);

                    // Final error tolerance evaluation.
                    error_tolerance = std::abs(h) * error_eot * error_scale / 752400.0;

                    if (error_tolerance <= 1.0){
                        // Successful step
                        break;
                    }

                    // Unsuccessful step
                    if (error_tolerance < 59049.0){
                        h_scale = 0.9 / std::pow(error_tolerance, 0.2);
                    }
                    else {
                        h_scale = 0.1;
                    }
                    h = h_scale * h;
                } // END SUBSTEPPING. Successful integration!
                // Increase the time
                t = t + h;

                // Copy the result.
                y[0] = s[0];
                y[1] = s[1];
                // Calculate the new derivative values;
                flt_pde(y, yp, interp_context);

                // Check immediately if we left the domain
                if (y[0] < m_r_min || y[0] > m_r_max){
                    left_computational_domain = true;
                    intersect=true;
                    continue;
                }
                if (y[1] < m_z_min || y[1] > m_z_max){
                    left_computational_domain = true;
                    intersect=true;
                    continue;
                }

                x1 = y[0] * cos(t);
                y1 = y[0] * sin(t);
                z1 = y[1];

                // Calculate length of field-line segment for the current step
                dx = (x1-ox);
                dy = (y1-oy);
                dz = (z1-oz);
                norm = std::sqrt(dx * dx + dy * dy + dz * dz);
                dx = dx / norm;
                dy = dy / norm;
                dz = dz / norm;
                total_length = total_length + norm;

                // No FLT check self intersection avoidance length
                if (total_length < self_intersection_avoidance_length){
                    // Instead of re-calculating the starting point of a FL,
                    // assign the end point from the current FL segment as the
                    // start point of the next FL segment. (To avoid
                    // re-computing it).
                    ox = x1;
                    oy = y1;
                    oz = z1;
                    continue;
                }
                // Set the dt to t_step (before we tried to avoid the
                // self_intersection_avoidance_length)
                dt = t_step * direction;

                // Prepare the struct for RayTracing
                rayHit.ray.org_x = ox;
                rayHit.ray.org_y = oy;
                rayHit.ray.org_z = oz;
                rayHit.ray.dir_x = dx;
                rayHit.ray.dir_y = dy;
                rayHit.ray.dir_z = dz;
                rayHit.ray.tnear = 0.0; // Segment going from origin point
                rayHit.ray.tfar = norm; // and is of length norm in (dx, dy,
                rayHit.ray.mask = 0;    // dz) direction
                rayHit.ray.flags = 0;
                rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
                rayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
                m_embree_obj->castRay(&rayHit, &rayContext);
                if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                    intersect=true;
                }

                // Instead of re-calculating the starting point of a FL, assign
                // the end point from the current FL segment as the start point
                //of the next FL segment. (To avoid re-computing it).
                ox = x1;
                oy = y1;
                oz = z1;

            } // END WHILE LOOP. EITHER INRESECTION OR MAX LENGTH ACCHIEVED

            // Register which geometry was hit.
            if (intersect){
                if (!left_computational_domain){
                    m_out_geom_hit_ids[i] = rayHit.hit.geomID;
                    m_out_prim_hit_ids[i] = rayHit.hit.primID;
                }
                else {
                    m_out_geom_hit_ids[i] = -2;
                    m_out_prim_hit_ids[i] = -2;
                }
            }
            else {
                // Since there was no hit, then we set the geomId to -1,
                // meaning the FL has been followed to the end of it's toroidal
                // time or maximum connection length. In case the FL leaves the
                // computational domain this is already taken care of.
                m_out_geom_hit_ids[i] = -1;
                m_out_prim_hit_ids[i] = -1;

            }
            m_out_fieldline_lengths[i] = total_length;
            } // END MAIN FOR LOOP
        delete interp_context;
    } // END OF PRAGMA OMP
}

void FLT::fehl_step(double y[2], double t, double h, double yp[2],
                    double k2[2], double k3[2], double k4[2], double k5[2],
                    double k6[2], double s[2], BI_DATA *interp_context){

    // Calculate the k1, k2, k3, k4, k5, k6 factors of the RKF4(5) method and
    // the approximate solution. The k* factors are used for error estimation
    // while the s holds the approximate solution. If the error is acceptable
    // then the approximate solution is the new solution.

    // The factors are calculated based on Table II in Fehlberg, NASA Technical
    // Report 287.

    // In this case we have a non-time dependant PDE

    // Use k6 as container for  calculating other k factors (used as argument
    // for the flt_pde method).
    double ch;

    // k1 is h flt_pde(t, y)

    // k2
    ch = h * 0.25;
    k6[0] = y[0] + ch * yp[0];
    k6[1] = y[1] + ch * yp[1];

    // Evaluate the function for k2
    // flt_pde(t + ch, k6, k2);
    flt_pde(k6, k2, interp_context);

    // k3
    ch = 3.0 * h / 32.0;
    k6[0] = y[0] + ch * (yp[0] + 3.0 * k2[0]);
    k6[1] = y[1] + ch * (yp[1] + 3.0 * k2[1]);

    // Evaluate the function for k3
    // flt_pde(t + 3.0*h/8.0, k6, k3, interp_context);
    flt_pde(k6, k3, interp_context);

    // k4
    ch = h / 2197.0;
    k6[0] = y[0] + ch * (1932.0 * yp[0] + (7296.0 * k3[0] - 7200.0 * k2[0]));
    k6[1] = y[1] + ch * (1932.0 * yp[1] + (7296.0 * k3[1] - 7200.0 * k2[1]));

    // Evaluate the function for k5
    // flt_pde(t + 12.0 * h / 13.0, k6, k4, interp_context);
    flt_pde(k6, k4, interp_context);

    // k5
    ch = h / 4104.0;
    k6[0] = y[0] + ch * ((8341.0 * yp[0] - 845.0 * k4[0]) + (29440.0 * k3[0] - 32832.0 * k2[0]));
    k6[1] = y[1] + ch * ((8341.0 * yp[1] - 845.0 * k4[1]) + (29440.0 * k3[1] - 32832.0 * k2[1]));

    // Evaluate the function at the new position
    // flt_pde(t + h, k6, k5, interp_context);
    flt_pde(k6, k5, interp_context);

    // Finally k6
    // The k2 factor is not used in the last evaluations, so we will use k2
    // here for storing result for k6.
    ch = h / 20520.0;
    k2[0] = y[0] + ch * (-6080.0 * yp[0] + (9295.0 * k4[0] - 5643.0 * k5[0]) + (41040.0 * k2[0] - 28352.0 * k3[0]));
    k2[1] = y[1] + ch * (-6080.0 * yp[1] + (9295.0 * k4[1] - 5643.0 * k5[1]) + (41040.0 * k2[1] - 28352.0 * k3[1]));

    // Evaluate the function for k6
    // flt_pde(t + 0.5 * h, k2, k6, interp_context);
    flt_pde(k2, k6, interp_context);

    // Approximate solution using up to order 6
    ch = h / 7618050.0;
    s[0] = y[0] + ch * (902880.0 * yp[0] + (3855735.0 * k4[0] - 1371249.0 * k5[0]) + (3953664.0 * k3[0] + 277020.0 * k6[0]));
    s[1] = y[1] + ch * (902880.0 * yp[1] + (3855735.0 * k4[1] - 1371249.0 * k5[1]) + (3953664.0 * k3[1] + 277020.0 * k6[1]));
}

void FLT::getFL(const double r, const double z, const double phi,
                 const int direction, std::vector<double>& storage,
                 const bool with_flt){
    // This function is used to obtain the fieldlines point.

    // Control options and values set by the user!!
    // Termination length of the field line. Meaning, tracing the field line
    // to this maximum length, unless other conditions terminate the trace
    // earlier.
    const double terminate_length = m_max_fieldline_length;
    // const int direction = m_directions[0];
    // Length where we don't check for RT intersections.
    const double self_intersection_avoidance_length = m_self_intersection_avoidance_length;
    // pre_t_step - Initial following of FL for avoiding self-intersection
    // const double pre_t_step = std::asin(0.1 * m_self_intersection_avoidance_length / y[0]) * direction;
    // t_step is the desired step we wish to make
    const double t_step = m_t_step;
    // Relative and absolute error
    const double relerr = m_relerr;
    const double abserr = m_abserr;
    // Direction - toroidal direction, 1 Counter-CW, -1 CW


    // Length of the followed field line so far.
    double total_length;
    bool intersect;
    double pre_t_step;

    // RKF4(5) variables.
    // h is the actual step performed during the RKF4(5) steps.
    double h;

    // Step increment arrays
    double k2[2];
    double k3[2];
    double k4[2];
    double k5[2];
    double k6[2];
    // Solution array
    double s[2];
    // Error values
    const double error_scale = 2.0/relerr;
    const double error_absolute = error_scale*abserr;
    double error_et;
    double error_ete;
    double error_eot;
    double error_tolerance;

    // Embree constructs
    RTCRayHit rayHit = RTCRayHit();
    RTCIntersectContext rayContext = RTCIntersectContext();
    rtcInitIntersectContext(&rayContext);

    // Ray tracing variables that are put into the RTCRayHit struct.
    double ox, oy, oz, x1, y1, z1, norm, dx, dy, dz;
    // y[2] and yp[2] correspond to (R, Z) and value of the derivative of the
    // function we wish to solve.
    double y[2], yp[2];
    // t equals to the parameteric time or toroidal angle
    // dt equals to the
    double t;
    double dt, h_scale;

    BI_DATA *interp_context = new BI_DATA();
    m_interp_psi->populateContext(interp_context);;
    // Initial values;
    y[0] = r;
    y[1] = z;
    flt_pde(y, yp, interp_context); // Calculate the derivative values

    if (with_flt && self_intersection_avoidance_length > 0.0){
        pre_t_step = std::asin(0.1 * self_intersection_avoidance_length / y[0]) * direction;
        dt = pre_t_step;
    }
    else {
        dt = t_step * direction;
    }
    t = phi;

    total_length = 0.0;

    // Termination flags for the main while loop
    intersect=false;

    // Calculate the (X, Y, Z) starting point of a FL segment
    ox = y[0] * std::cos(t);
    oy = y[0] * std::sin(t);
    oz = y[1];

    // Put first point into storage
    storage.push_back(r);
    storage.push_back(z);
    storage.push_back(phi);

    while (total_length <= terminate_length && !intersect){
        // new_time = t + pre_t_step;
        // while (0 < (new_time - t) * direction){
        //     flag = m_rkf45_solvers[omp_thread]->r8_rkf45(y, yp, &t,
        //                                                 new_time, &relerr,
        //                                                 abserr, flag);
        //     flag = 2; // Until we find an actual problem, ignore messages from
        //               // rkf45 as they do not indicate a problem*/
        // }

        // Obtain the derivatives
        // flt_pde(t, y, yp);

        // We do not allow limiting the step size. Even if the step size is
        // 1e-6, the actual time would have to be greater than 173215370,
        // which corresponds to more than 27 million revolutions around
        // the tokamak.

        // Single RKF45 step
        // if (std::abs(pre_t_step) <= 26.0 * EPSILON * std::abs(t)){
        //     t = t + pre_t_step;
        //     y[0] = y[0] + dt * yp[0];
        //     y[1] = y[1] + dt * yp[1];
        // }

        h = dt;

        // Start the integration
        while (true){
            // Advance an approximate solution over one step of length h
            // If there is no need to lower the step size then the solution
            // in variable s is accepted.
            fehl_step(y, t, h, yp, k2, k3, k4, k5, k6, s, interp_context);
            // Compute and test allowable tolerances versus local error
            // estimates and remove scaling of tolerances. The relative
            // error is measured with respect to the average of the
            // magnitudes of the solution at the beginning and end of the
            // step.

            // FORMULA 2 Table III in Fehlberg for determining truncation
            // error.

            // Get average of the function magnitudes at t and t + h
            error_et = std::abs(y[0]) + std::abs(k2[0]) + error_absolute;
            // We assume that the errors set are positive!

            // Truncation error
            error_ete = std::abs(-2090.0*yp[0] + 21970.0 * k4[0] - 15048.0 * k5[0] + 22528.0 * k3[0] - 27360.0 * k6[0]);

            // Get the relative error
            error_eot = error_ete/error_et;

            // Do the same at second dimension
            error_et = std::abs(y[1]) + std::abs(k2[1]) + error_absolute;
            // We assume that the errors set are positive!

            // Truncation error
            error_ete = std::abs(-2090.0*yp[1] + 21970.0 * k4[1] - 15048.0 * k5[1] + 22528.0 * k3[1] - 27360.0 * k6[1]);

            // Take the max relative error of the two
            error_eot = std::max(error_eot, error_ete/error_et);

            // Final error tolerance evaluation.
            error_tolerance = std::abs(h) * error_eot * error_scale / 752400.0;

            if (error_tolerance <= 1.0){
                // Successful step
                break;
            }

            // Unsuccessful step
            if (error_tolerance < 59049.0){
                h_scale = 0.9 / std::pow(error_tolerance, 0.2);
            }
            else {
                h_scale = 0.1;
            }
            h = h_scale * h;
        } // END SUBSTEPPING. Successful integration!
        // Increase the time
        t = t + h;

        // Copy the result.
        y[0] = s[0];
        y[1] = s[1];
        // Calculate the new derivative values;
        flt_pde(y, yp, interp_context);

        if (y[0] < m_r_min || y[0] > m_r_max){
            intersect=true;
            continue;
        }
        if (y[1] < m_z_min || y[1] > m_z_max){
            intersect=true;
            continue;
        }

        x1 = y[0] * cos(t);
        y1 = y[0] * sin(t);
        z1 = y[1];
        storage.push_back(y[0]);
        storage.push_back(y[1]);
        storage.push_back(t);
        // Calculate length of field-line segment for the current step
        dx = (x1-ox);
        dy = (y1-oy);
        dz = (z1-oz);
        norm = std::sqrt(dx * dx + dy * dy + dz * dz);
        dx = dx / norm;
        dy = dy / norm;
        dz = dz / norm;
        // Set the original point of the ray into the struct here even if do
        // not need FLT.
        rayHit.ray.org_x = ox;
        rayHit.ray.org_y = oy;
        rayHit.ray.org_z = oz;

        // Instead of re-calculating the starting point of a FL, assign the end
        // point from the current FL segment as the start point of the next
        // FL segment. (To avoid re-computing it).
        ox = x1;
        oy = y1;
        oz = z1;

        total_length = total_length + norm;

        // No FLT check self intersection avoidance length
        if (with_flt && total_length < self_intersection_avoidance_length){
            continue;
        }

        // Set the dt to t_step (before we tried to avoid the
        // self_intersection_avoidance_length)
        dt = t_step * direction;

        // No RayTrace check for collision with geometry.
        if (!with_flt){
            continue;
        }

        // Prepare the struct for RayTracing
        rayHit.ray.dir_x = dx;
        rayHit.ray.dir_y = dy;
        rayHit.ray.dir_z = dz;
        rayHit.ray.tnear = 0.0;
        rayHit.ray.tfar = norm;
        rayHit.ray.mask = 0;
        rayHit.ray.flags = 0;
        rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        rayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

        m_embree_obj->castRay(&rayHit, &rayContext);
        if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
            intersect=true;
        }


    } // END WHILE LOOP. EITHER INRESECTION OR MAX LENGTH ACCHIEVED
    delete interp_context;
}

void FLT::flt_pde(double y[2], double yp[2], BI_DATA *context){
    double factor;

    context->r = y[0] - m_r_move;
    context->z = y[1] - m_z_move;
    m_interp_psi->getValues(context);
    factor = y[0] / m_vacuum_fpol;

    yp[0] = -context->valdy * factor;
    yp[1] = context->valdx * factor;
}