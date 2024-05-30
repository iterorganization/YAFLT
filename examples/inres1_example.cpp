#define _USE_MATH_DEFINES
#include <flt.hpp>
#include <polylineWriter.cpp>

#include <iostream> // std::cout

#include <accell_embree.hpp>
#include <data_for_shadow.h>
#include <data_for_spline.h>
#include <data_for_target.h>

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

    // Create and load the Embree object
    EmbreeAccell *embreeObj = new EmbreeAccell();
    // Create pointers

    float* shadow_vertices = svec;
    unsigned int* shadow_triangles = stri;

    std::cout << "Loading shadow mesh... ";
    embreeObj->commitMesh(shadow_vertices, n_svec, shadow_triangles, n_stri);
    std::cout << "Done!" << std::endl;

    // Create the fieldline object
    FLT *obj = new FLT();

    // Set the Embree object
    obj->setEmbreeObj(embreeObj);

    // Load the EQDSK data

    // Necessary to input the number of points of R, Z grid.
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

    // Prepare interpolation objects
    obj->prepareInterpolation();


    // Set ODE parameters
    obj->setAbsError(1e-4);
    obj->setRelError(1e-4);

    // Inres -6 mm radial shift and 0 vertical shift
    obj->setShift(-0.006, 0.0);

    // Set the integration toroidal angle step
    obj->setDesiredStep(0.01);
    obj->setMaximumFieldlineLength(4.3); // meters
    obj->setSelfIntersectionAvoidanceLength(0.001);


    /*Variables*/
    double p1[3], p2[3], p3[3];
    std::vector<double> barycenter;

    std::vector<double> normal, magneticVector;
    magneticVector.resize(3);
    double br, bz, bphi;
    double dot_product;
    bool is_hit;
    int direction;

    double FPol_vaccuum = FPOLs[FPOLs.size() - 1];
    int default_direction = (FPol_vaccuum < 0.0) ? -1 : 1;

    int nTriangle = triv0.size();

    std::vector<double> conlen;
    std::vector<int> mask;

    conlen.resize(nTriangle);
    mask.resize(nTriangle);

    int nShadowed = 0;

    std::vector<double> points;
    std::vector<int> directions;
    for(int i = 0; i < nTriangle; i++){
        std::cout << std::endl;
        std::cout << "Triangle No. " << i << " has vertices:" << std::endl;
        std::cout << "Vertex No. " << triv0[i] << ": {" << vx[triv0[i]-1] << " " << vy[triv0[i]-1] << " " << vz[triv0[i]-1] << "}" << std::endl;
        std::cout << "Vertex No. " << triv1[i] << ": {" << vx[triv1[i]-1] << " " << vy[triv1[i]-1] << " " << vz[triv1[i]-1] << "}" << std::endl;
        std::cout << "Vertex No. " << triv2[i] << ": {" << vx[triv2[i]-1] << " " << vy[triv2[i]-1] << " " << vz[triv2[i]-1] << "}" << std::endl;

        p1[0] = vx[triv0[i]-1]; p1[1] = vy[triv0[i]-1]; p1[2] = vz[triv0[i]-1];
        p2[0] = vx[triv1[i]-1]; p2[1] = vy[triv1[i]-1]; p2[2] = vz[triv1[i]-1];
        p3[0] = vx[triv2[i]-1]; p3[1] = vy[triv2[i]-1]; p3[2] = vz[triv2[i]-1];
        barycenter = getBaryCenter(p1, p2, p3);

        // Get the barycenter of the triangles and convert from mm to m.
        br = 0.001 * sqrt(pow(barycenter[0], 2.0) +  pow(barycenter[1], 2.0));
        bz = 0.001 * barycenter[2];
        bphi = atan2(barycenter[1], barycenter[0]);
        // Store the points
        points.push_back(br);
        points.push_back(bz);
        points.push_back(bphi);

        /*Figure out direction*/

        /*Get triangle normal*/
        normal = getTriangleNormal(p1, p2, p3);
        /*Get direction and angle from the magnetic */
        obj->getBCart(br, bz, bphi, magneticVector);

        dot_product = std::inner_product(normal.begin(), normal.end(), magneticVector.begin(), 0.0);

        // Sign of FPol tell us the initial starting direction.
        direction = default_direction;
        // But what if the normal of the triangle is facing the other direction?
        direction = (dot_product < 0.0) ? -direction : direction;
        directions.push_back(direction);
        printf("i=%d br=%f bz=%f bphi=%f direction=%d\n", i, br, bz, bphi, direction);
    }


    obj->setStartingFLDirection(directions);
    obj->setPoints(points);

    // Set number of threads
    obj->setNumberOfThreads(16);

    // Run the FLT
    obj->runFLT();

    double fieldline_length=0.0;
    for(int i = 0; i < obj->m_out_fieldline_lengths.size(); i++){
        fieldline_length = obj->m_out_fieldline_lengths[i];
        printf("%d length=%f geom=%d prim=%d\n", i, fieldline_length, obj->m_out_geom_hit_ids[i], obj->m_out_prim_hit_ids[i]);
        conlen[i] = fieldline_length;

        is_hit=false;
        if (fieldline_length < 4.3){
            is_hit = true;
            nShadowed += is_hit;
        }
        mask[i] = is_hit;

    }


    std::cout << "Number of shadowed triangles: " << nShadowed << std::endl;
    float percentage = 100.0 * nShadowed / nTriangle;
    std::cout << "Percentage: " << percentage << "%." << std::endl;

    int nPoints = vx.size();
    writePolyData2VTK(nPoints, vx, vy, vz, nTriangle, triv0, triv1, triv2,
                      conlen, mask, "inres1.vtk");
    return 0;
}