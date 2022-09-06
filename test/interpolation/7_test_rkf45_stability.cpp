//TODO
#include <eqdsk_drsep10_data.hpp>
#include <rkf45.hpp>
#include <bicubic.hpp>
int main(){

    RKF45 *solver = new RKF45();
    BICUBIC_INTERP *naive_interp = new BICUBIC_INTERP();

    solver->set_interpolator(naive_interp);

    int n_cols = drsep10_Rs.size();
    int n_rows = drsep10_Zs.size();

    std::vector<std::vector<double>> reshaped_Psi;
    std::vector<double> buffer;
    for (int i=0; i<n_rows; i++){
        buffer.clear();
        for (int j=0; j<n_cols; j++){
            buffer.push_back(drsep10_Psi[i * n_cols + j]);
        }
        // reshaped_Psi.insert(reshaped_Psi.begin(), buffer);
        reshaped_Psi.push_back(buffer);
    }
    naive_interp->prepareContainers();
    naive_interp->setArrays(drsep10_Rs, drsep10_Zs, reshaped_Psi);

    return 0;
}