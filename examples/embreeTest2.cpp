#include <iostream>

#include <accell_embree.hpp>
#include <data_for_shadow.h>

int main(){

    EmbreeAccell *embreeObj = new EmbreeAccell();
    // Create pointers

    float* shadow_vertices = svec;
    unsigned int* shadow_triangles = stri;

    std::cout << "Loading shadow mesh... ";
    embreeObj->commitMesh(shadow_vertices, n_svec, shadow_triangles, n_stri);
    std::cout << "Done!" << std::endl;
    return 0;
}
