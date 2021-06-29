#include <flt.hpp>
#include <polylineWriter.cpp>

#include <iostream>    // std::cout

#include <data_for_spline.h>
#include <data_for_target.h>
#include <vector>

std::vector<double> getBaryCenter(double p1[3], double p2[3], double p3[3]){
    std::vector<double> bc;
    bc.push_back((p1[0] + p2[0] + p3[0]) / 3);
    bc.push_back((p1[1] + p2[1] + p3[1]) / 3);
    bc.push_back((p1[2] + p2[2] + p3[2]) / 3);
    return bc;
}

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
    obj->setTimeSpan(1.0, 1e-1);

    double p1[3], p2[3], p3[3];
    std::vector<double> barycenter;
    double br, bz, bphi;
    //std::vector<double> p1, p2, p3;
    int nTriangle = triv0.size();
    for(int i = 1; i < nTriangle + 1; i++){
        /*Run only for 10% of triangles*/
        if(i % (triv0.size() / 10) == 0){
            printf("Triangle No. %d has vertices:\n", i);
            printf("Vertex No. %d: {%f, %f, %f}\n",  triv0[i-1], vx[triv0[i-1]-1], vy[triv0[i-1]-1], vz[triv0[i-1]-1]);
            printf("Vertex No. %d: {%f, %f, %f}\n",  triv1[i-1], vx[triv1[i-1]-1], vy[triv1[i-1]-1], vz[triv1[i-1]-1]);
            printf("Vertex No. %d: {%f, %f, %f}\n",  triv2[i-1], vx[triv2[i-1]-1], vy[triv2[i-1]-1], vz[triv2[i-1]-1]);

            p1[0] = vx[triv0[i-1]-1]; p1[1] = vy[triv0[i-1]-1]; p1[2] = vz[triv0[i-1]-1];
            p2[0] = vx[triv1[i-1]-1]; p2[1] = vy[triv1[i-1]-1]; p2[2] = vz[triv1[i-1]-1];
            p3[0] = vx[triv2[i-1]-1]; p3[1] = vy[triv2[i-1]-1]; p3[2] = vz[triv2[i-1]-1];
            barycenter = getBaryCenter(p1, p2, p3);
            br = 0.001 * sqrt(pow(barycenter[0], 2.0) +  pow(barycenter[1], 2.0));
            bz = 0.001 * barycenter[2];
            bphi = atan2(barycenter[1], barycenter[0]);

            obj->setIV(br, bz, bphi);

            std::vector<double> result;
            obj->getFL(result);
            if (result.size() == 0){
                std::cout << "No results stored. Stopping" << std::endl;
            }

            write2VTK(result, "fl_cpp2_" + std::to_string(i) + ".vtk");
            writeTri2VTK(p1, p2, p3, "target_tri_" + std::to_string(i) + ".vtk");
        }
    }

    return 0;
}