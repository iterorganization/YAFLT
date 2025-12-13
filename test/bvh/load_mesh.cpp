#include <iostream>

#include <tlas.hpp>
#include <data_for_shadow.h>

int main(){
    TLAS *tlas = new TLAS();
    // Create pointers

    float* shadow_vertices = svec;
    unsigned int* shadow_triangles = stri;

    std::cout << "Loading shadow mesh... " << std::flush;
    tlas->commitMesh(shadow_vertices, n_svec, shadow_triangles, n_stri);
    std::cout << "Done!" << std::endl;

    delete tlas;
    return 0;
}
