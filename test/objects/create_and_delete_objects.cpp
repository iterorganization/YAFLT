#include <bicubic.hpp>
#include <accell_embree.hpp>
#include <cstdio>
#include <flt.hpp>
#include <rkf45.hpp>

// run with valgrind --leak-check=full ./create_and_delete_objects

int main(){

    EmbreeAccell *embree_obj;
    FLT *flt_obj;
    BICUBIC_INTERP *bicubic_obj;

    int N = 10;
    for (int i = 0; i < N; ++i)
    {
        embree_obj = new EmbreeAccell();
        flt_obj = new FLT();
        bicubic_obj = new BICUBIC_INTERP();

        delete embree_obj;
        delete flt_obj;
        delete bicubic_obj;
    }
    printf("If you see this, then it's okay\n");
    return 0;
}