#include <tlas.hpp>

#include <numeric> // std::inner_product
#include <limits> // std::numeric_limits

int main(){

    TLAS *tlas = new TLAS();

    // Simple triangle with 6 points

    float vertices[18];
    vertices[0] = -0.5f; vertices[1] = -0.5f; vertices[2] = 0.f;
    vertices[3] = -1.f; vertices[4] = 0.f; vertices[5] = 0.f;
    vertices[6] = 0.f; vertices[7] = -1.f; vertices[8] = 0.f;
    vertices[9] = 0.f; vertices[10] = 0.f; vertices[11] = 0.f;
    vertices[12] = 1.f; vertices[13] = 0.f; vertices[14] = 0.f;
    vertices[15] = 0.f; vertices[16] = 1.f; vertices[17] = 0.f;

    // Two triangles
    unsigned indices[6];
    indices[0] = 0; indices[1] = 1; indices[2] = 2;
    indices[3] = 3; indices[4] = 4; indices[5] = 5;

    int vertices_size = 18;
    int triangles_size = 2;
    double tfar = 1.0E30;
    bool intersect;

    tlas->commitMesh(vertices, vertices_size, indices, triangles_size);

    tinybvh::Ray rayHit = tinybvh::Ray();
    // This one will hit
    rayHit.O.x = 0;
    rayHit.O.y = 0;
    rayHit.O.z = -1;
    rayHit.D.x = 0;
    rayHit.D.y = 0;
    rayHit.D.z = 1;
    rayHit.hit.t = tfar;

    tlas->castRay(&rayHit);
    intersect = rayHit.hit.t < tfar;

    if (intersect != true){
        printf("castRay failed detecting ");
        return 1;
    }

    // Now try it the other way, using the castRay
    tlas->castRay(0, 0, -1, 0, 0, 1, tfar);
    if (tlas->checkIfHit() != true){
        printf("castRay failed detecting intersection");
        return 1;
    }

    // This one will miss!
    rayHit.O.x = 1;
    rayHit.O.y = 1;
    rayHit.O.z = -1;
    rayHit.D.x = 0;
    rayHit.D.y = 0;
    rayHit.D.z = 1;
    rayHit.hit.t = tfar;

    tlas->castRay(&rayHit);
    intersect = rayHit.hit.t < tfar;
    if (intersect == true){
        printf("castRay falsly detects intersection");
        return 1;
    }
    // Now try it the other way
    // Second one NOT!
    tlas->castRay(1, 1, -1, 0, 0, 1, tfar);
    if (tlas->checkIfHit() == true){
        printf("castRay falsly detects intersection");
        return 1;
    }

    delete tlas;
    return 0;
}
