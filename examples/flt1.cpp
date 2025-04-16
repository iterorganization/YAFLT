// This example shows how to obtain an arbitrary fieldline. No intersection
// tests.

#include <flt.hpp>
#include <polylineWriter.cpp>

#include <iostream>    // std::cout

#include <data_for_spline.h>

int main(){
    FLT *obj = new FLT();

    // Set input data
    obj->setPoloidalMagneticFlux(Rs, Zs, PSIs);
    obj->setVacuumFPOL(FPOLs[FPOLs.size() - 1]);

    // Set ODE parameters
    obj->setAbsError(1e-4);
    obj->setRelError(1e-4);

    obj->setDesiredStep(0.01);

    std::vector<double> result;
    int direction=1;
    bool with_flt=false; // No shadow geometry used anyway;
    obj->setMaximumFieldlineLength(2565.4);
    obj->getFL(4.043846846215178, 0.284476, 0.09093153458023334, direction, result, with_flt);
    if (result.size() == 0){
        std::cout << "No results stored. Stopping" << std::endl;
    }

    write2VTK(result, "fl_cpp2_1.vtk");

    return 0;
}