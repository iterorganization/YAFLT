#include <accell_embree.hpp>

#include <iostream>    // std::cout
#include <numeric> // std::inner_product
#include <limits> // std::numeric_limits

int main(){

    EmbreeAccell *obj2 = new EmbreeAccell();
    // obj2->createDevice();

    float vertices[18];
    vertices[0] = -0.5f; vertices[1] = -0.5f; vertices[2] = 0.f;
    vertices[3] = -1.f; vertices[4] = 0.f; vertices[5] = 0.f;
    vertices[6] = 0.f; vertices[7] = -1.f; vertices[8] = 0.f;
    vertices[9] = 0.f; vertices[10] = 0.f; vertices[11] = 0.f;
    vertices[12] = 1.f; vertices[13] = 0.f; vertices[14] = 0.f;
    vertices[15] = 0.f; vertices[16] = 1.f; vertices[17] = 0.f;

    unsigned indices[6];
    indices[0] = 0; indices[1] = 1; indices[2] = 2;
    indices[3] = 3; indices[4] = 4; indices[5] = 5;

    int vertices_size = 18;
    int triangles_size = 2;
    double tnear = 0;
    double tfar = std::numeric_limits<float>::infinity();
    bool intersect;

    obj2->commitMesh(vertices, vertices_size, indices, triangles_size);
    /*Even if running sequentially the following function must be called.*/
    obj2->prepareThreadContainers();
    obj2->castRay(0, 0, -1, 0, 0, 1, tnear, tfar);

    intersect = obj2->checkIfHit();
      if (intersect == true)
         std::cout << "intersect: " << intersect << std::endl;
      else
        std::cout << "intersect: 0" << std::endl;

    obj2->castRay( 1, 1, -1, 0, 0, 1, tnear, tfar);

    intersect = obj2->checkIfHit();
      if (intersect == true)
         std::cout << "intersect: " << intersect << std::endl;
      else
        std::cout << "intersect: 0" << std::endl;

    return 0;
}
