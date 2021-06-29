#include <flt.hpp>
#include <polylineWriter.cpp>

#include <iostream>    // std::cout

#include <data_for_spline.h>

int main(){
    FLT *obj = new FLT();

    int Nr = Rs.size();
    int Nz = Zs.size();


    // Set input data
    obj->setNDIM(Nr, Nz);

    obj->setRARR(Rs);
    obj->setZARR(Zs);
    obj->setPSI(PSIs);
    obj->setFARR(SIMAG_BDRYs);

    obj->setFPOL(FPOLs);
    obj->setVacuumFPOL(FPOLs[FPOLs.size() - 1]);

    // Init and run
    obj->prepareInterpolation();
    /*Even if running sequentially the following function must be called.*/
    obj->prepareThreadContainers();

    // Set ODE parameters
    obj->setAbsError(1e-4);
    obj->setRelError(1e-4);
    obj->setIV(4.043846846215178, 0.284476, 0.09093153458023334);
    obj->setTimeSpan(471.23, 1e-2);

    std::vector<double> result;
    obj->getFL(result);
    if (result.size() == 0){
        std::cout << "No results stored. Stopping" << std::endl;
    }

    write2VTK(result, "fl_cpp1.vtk");

    return 0;
}