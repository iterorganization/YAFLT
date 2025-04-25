// An example limiter equilibrium test where we obtain the lengths of the
// fieldlines and the mask showing which fieldlines are effectively "shadowed",
// or which fieldlines intersect early with a geometry we loaded and we call as
// shadow geometry.

#include <cstdio>
#include <flt.hpp>

#include <accell_embree.hpp>
#include <data_for_spline.h>
#include <data_for_shadow.h>
#include <omp.h>


int main(){

    setbuf(stdout, NULL);

    // Create the Embree object
    EmbreeAccell *embreeObj = new EmbreeAccell();

    float* shadow_vertices = svec;
    unsigned int* shadow_triangles = stri;
    embreeObj->commitMesh(shadow_vertices, n_svec, shadow_triangles, n_stri);
    // Create the fieldline object
    FLT *obj = new FLT();

    // Set the Embree object - It's empty but it's okay
    obj->setEmbreeObj(embreeObj);

    // Load the EQDSK data - limiter equilibrium.

    // Set input data
    obj->setPoloidalMagneticFlux(Rs, Zs, PSIs);
    obj->setVacuumFPOL(FPOLs[FPOLs.size() - 1]);

    // Set ODE parameters
    obj->setAbsError(1e-4);
    obj->setRelError(1e-4);

    // Set the integration toroidal angle step
    obj->setDesiredStep(0.05);
    obj->setMaximumFieldlineLength(100.0); // meters
    obj->setSelfIntersectionAvoidanceLength(0.001);

    // Create the points
    std::vector<double> points;
    std::vector<int> directions;

    // Creating points inside plasma core so that every fieldline is followed
    // Until it's length.
    double Rl, Rr, Zu, Zl, rmin, zmin, rdiff, zdiff;
    Rl = 5.0,
    Zl = -0.5;
    Rr = 7.0;
    Zu = 1.5;

    rmin = Rl;
    rdiff = Rr - Rl;
    zmin = Zl;
    zdiff = Zu - Zl;

    int n_points = 1000000;
    printf("Creating points...");
    for (int i = 0; i < n_points; ++i)
    {
        points.push_back(rmin + rdiff * rand() / RAND_MAX);
        points.push_back(zmin + zdiff * rand() / RAND_MAX);
        points.push_back(0.0);
        directions.push_back(1.0);
    }
    printf("Done.\n");

    obj->setStartingFLDirection(directions);
    obj->setPoints(points);

    std::vector<double> times;

    int MAX_CPUS = omp_get_num_procs();

    double begin, end;

    double one_core_time;

    obj->setNumberOfThreads(1);

    begin = omp_get_wtime();
    obj->runFLT();
    end = omp_get_wtime();

    one_core_time = end - begin;
    printf("One core time: %f\n", one_core_time);

    double max_core_time;

    obj->setNumberOfThreads(MAX_CPUS);
    begin = omp_get_wtime();
    obj->runFLT();
    end = omp_get_wtime();

    max_core_time = end - begin;

    printf("Max (%d) core time: %f\n", 16, max_core_time);
    // printf("Speedup at %d cores: %f\n", MAX_CPUS, one_core_time / max_core_time);

    delete embreeObj;
    delete obj;

    return 0;
}