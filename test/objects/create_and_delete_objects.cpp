#include <bicubic.hpp>
#include <tlas.hpp>
#include <cstdio>
#include <flt.hpp>
#include <rkf45.hpp>

// run with valgrind --leak-check=full ./create_and_delete_objects

int main(){

    TLAS *tlas_obj;
    FLT *flt_obj;
    BICUBIC_INTERP *bicubic_obj;

    int N = 10;
    for (int i = 0; i < N; ++i)
    {
        tlas_obj = new TLAS();
        flt_obj = new FLT();
        bicubic_obj = new BICUBIC_INTERP();

        delete tlas_obj;
        delete flt_obj;
        delete bicubic_obj;
    }
    return 0;
}