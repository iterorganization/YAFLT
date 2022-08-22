#define _USE_MATH_DEFINES
#include <flt.hpp>
#include <polylineWriter.cpp>

#include <iostream>    // std::cout

#include <data_for_spline.h>
#include <data_for_target.h>
#include <vector>
#include <numeric> // std::inner_product
#include <algorithm> // std::transform
#include <functional> // std::negate

std::vector<double> getBaryCenter(double p1[3], double p2[3], double p3[3]){
    std::vector<double> bc;
    bc.push_back((p1[0] + p2[0] + p3[0]) / 3);
    bc.push_back((p1[1] + p2[1] + p3[1]) / 3);
    bc.push_back((p1[2] + p2[2] + p3[2]) / 3);
    return bc;
}

std::vector<double> getTriangleNormal(double p1[3], double p2[3], double p3[3]){
    std::vector<double> normal;

    double a1,a2,a3,b1,b2,b3,norm2;

    a1 = p2[0] - p1[0];
    a2 = p2[1] - p1[1];
    a3 = p2[2] - p1[2];

    b1 = p3[0] - p1[0];
    b2 = p3[1] - p1[1];
    b3 = p3[2] - p1[2];

    normal.push_back(a2*b3 - a3*b2);
    normal.push_back(a3*b1 - a1*b3);
    normal.push_back(a1*b2 - a2*b1);
    norm2 = sqrt(std::inner_product(normal.begin(), normal.end(), normal.begin(), 0.0));
    normal[0] = normal[0] / norm2;
    normal[1] = normal[1] / norm2;
    normal[2] = normal[2] / norm2;
    return normal;
}

double angleBetweenVectors(std::vector<double> x1, std::vector<double> x2,
                           int direction=1){

    if (direction == -1){
        /*When the dot product between vectors is negative flip one of the
          vectors.*/
        std::transform(x2.cbegin(),x2.cend(),x2.begin(),std::negate<double>());
    }
    double dot = std::inner_product(x1.begin(), x1.end(), x2.begin(), 0.0);
    double sq1 = std::inner_product(x1.begin(), x1.end(), x1.begin(), 0.0);
    double sq2 = std::inner_product(x2.begin(), x2.end(), x2.begin(), 0.0);

    double angle = acos(dot / sqrt(sq1 * sq2));
    return angle;
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

    // Inres -6 mm radial shift and 0 vertical shift
    obj->setShift(-0.006, 0.0);

    double p1[3], p2[3], p3[3];
    std::vector<double> barycenter;

    std::vector<double> result;

    /*Variables for direction and incident angle*/
    std::vector<double> normal, magneticVector;
    magneticVector.resize(3);
    double br, bz, bphi;
    double dot_product;
    double angle_rad, angle_deg, angle_inc_deg;
    double fpolVacuum=FPOLs[FPOLs.size()-1];

    int direction;
    // Sign of FPol tell us the initial starting direction.
    direction = (fpolVacuum < 0.0) ? -1 : 1;

    //std::vector<double> p1, p2, p3;
    int nTriangle = triv0.size();
    for(int i = 0; i < nTriangle; i++){
        /*Run only for 10% of triangles*/
        std::cout << "Triangle No. " << i << " has vertices:" << std::endl;
        std::cout << "Vertex No. " << triv0[i] << ": {" << vx[triv0[i]-1] << " " << vy[triv0[i]-1] << " " << vz[triv0[i]-1] << "}" << std::endl;
        std::cout << "Vertex No. " << triv1[i] << ": {" << vx[triv1[i]-1] << " " << vy[triv1[i]-1] << " " << vz[triv1[i]-1] << "}" << std::endl;
        std::cout << "Vertex No. " << triv2[i] << ": {" << vx[triv2[i]-1] << " " << vy[triv2[i]-1] << " " << vz[triv2[i]-1] << "}" << std::endl;

        p1[0] = vx[triv0[i]-1]; p1[1] = vy[triv0[i]-1]; p1[2] = vz[triv0[i]-1];
        p2[0] = vx[triv1[i]-1]; p2[1] = vy[triv1[i]-1]; p2[2] = vz[triv1[i]-1];
        p3[0] = vx[triv2[i]-1]; p3[1] = vy[triv2[i]-1]; p3[2] = vz[triv2[i]-1];
        barycenter = getBaryCenter(p1, p2, p3);
        br = 0.001 * sqrt(pow(barycenter[0], 2.0) +  pow(barycenter[1], 2.0));
        bz = 0.001 * barycenter[2];
        bphi = atan2(barycenter[1], barycenter[0]);

        obj->setIV(br, bz, bphi);

        /*Figure out direction*/

        /*Get triangle normal*/
        normal = getTriangleNormal(p1, p2, p3);
        /*Get direction and angle from the magnetic */
        obj->getBCart(br, bz, bphi, magneticVector);

        std::cout << "Normal vector: " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
        std::cout << "Magnetic vector: " << magneticVector[0] << " " << magneticVector[1] << " " << magneticVector[2] << std::endl;

        dot_product = std::inner_product(normal.begin(), normal.end(), magneticVector.begin(), 0.0);
        std::cout << "Dot product: " << dot_product << std::endl;

        // But what if the normal of the triangle is facing the other direction?
        direction = (dot_product < 0.0) ? -direction : direction;

        std::cout << "Starting direction - CW: 1, CCW: 0 -> " << direction << std::endl;

        /*Get incident angle*/
        angle_rad = angleBetweenVectors(normal, magneticVector, direction);
        angle_deg = angle_rad * 180.0 / M_PI;

        angle_inc_deg = (angle_deg > 90.0) ? angle_deg - 90.0 : 90.0 - angle_deg;

        std::cout << "Angle in radians: " << angle_rad << std::endl;
        std::cout << "Angle in degrees: " << angle_deg << std::endl;
        std::cout << "Angle in degrees (incident): " << angle_inc_deg << std::endl;


        /*Apply direction */

        obj->setDirection(direction);
        obj->setTimeSpan(1.0, 1e-2);
        result.clear();
        result.push_back(br);
        result.push_back(bz);
        result.push_back(bphi);
        obj->getFL(result);
        if (result.size() == 3){
            std::cout << "No results stored. Stopping" << std::endl;
        }
        else {
            if (i == 1683){
                // Triangle max Q for inres
                write2VTK(result, "fl_cpp2_" + std::to_string(i) + ".vtk");
                writeTri2VTK(p1, p2, p3, "target_tri_" + std::to_string(i) + ".vtk");
            }
        }

    }

    return 0;
}