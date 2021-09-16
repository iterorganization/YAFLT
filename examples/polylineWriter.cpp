#include <vector>
#include <string>

#include <cmath>       // Sin and cos
#include <iostream>    // std::cout
#include <fstream>     // std::ofstream
#include <iomanip>     // std::setprecision

void write2VTK(std::vector<double> data, std::string filePath, int mul=1000){
    std::ofstream outFile;

    // Container contains 3D vectors.
    double x,y,z;
    int size = (int) data.size() / 3;

    outFile.open(filePath, std::ios::trunc);
    outFile << "# vtk DataFile Version 2.0\n";
    outFile << "comment\n";
    outFile << "ASCII\n\n";
    outFile << "DATASET UNSTRUCTURED_GRID\n";

    outFile << "POINTS"
          << std::setw(9)  << size << "  "
          << " double\n";

    for (int i=0; i<size;i++){
        x = mul * data[i * 3] * cos(data[i*3 + 2]);
        y = mul * data[i * 3] * sin(data[i*3 + 2]);
        z = mul * data[i * 3 + 1];
        outFile << std::setw(12) << std::setprecision(10) << x  << "  "
                << std::setw(12) << std::setprecision(10) << y  << "  "
                << std::setw(12) << std::setprecision(10) << z  << "\n";
    }

    outFile << "CELLS"
          << std::setw(12)  << size - 1  << "  "
          << std::setw(12)  << (size - 1) * 3  << "  " << "\n";
    for ( int i = 0; i < size - 1; i++ )
    {
    outFile << "2"
            << std::setw(12) << i  << "  "
            << std::setw(12) << i+1  << "\n";
    }

    outFile << "CELL_TYPES"
          << std::setw(9)  << size  << "  " << "\n";
    for ( int i = 0; i < size; i++ )
    {
        outFile << "3" << "\n";
    }
    outFile.close();


}


void writeTri2VTK(double p1[3], double p2[3], double p3[3], std::string filePath){
    std::ofstream outFile;

    outFile.open(filePath, std::ios::trunc);
    outFile << "# vtk DataFile Version 2.0\n";
    outFile << "comment\n";
    outFile << "ASCII\n\n";
    outFile << "DATASET UNSTRUCTURED_GRID\n";

    outFile << "POINTS"
            << std::setw(9)  << 3 << "  "
            << " double\n";

    outFile << std::setw(12) << std::setprecision(10) << p1[0]  << "  "
            << std::setw(12) << std::setprecision(10) << p1[1]   << "  "
            << std::setw(12) << std::setprecision(10) << p1[2]   << "\n";

    outFile << std::setw(12) << std::setprecision(10) << p2[0]  << "  "
            << std::setw(12) << std::setprecision(10) << p2[1]   << "  "
            << std::setw(12) << std::setprecision(10) << p2[2]   << "\n";

    outFile << std::setw(12) << std::setprecision(10) << p3[0]  << "  "
            << std::setw(12) << std::setprecision(10) << p3[1]   << "  "
            << std::setw(12) << std::setprecision(10) << p3[2]   << "\n";

    outFile << "CELLS"
            << std::setw(12)  << 3  << "  "
            << std::setw(12)  << 9  << "  " << "\n";

    for ( int i = 0; i < 3; i++ )
    {
        outFile << "2"
                << std::setw(12) << i  << "  "
                << std::setw(12) << i+1  << "\n";
    }

    outFile << "CELL_TYPES"
            << std::setw(9)  << 3  << "  " << "\n";

    for ( int i = 0; i < 3; i++ )
    {
        outFile << "3" << "\n";
    }
    outFile.close();
}

// Used by inres1_example
void writePolyData2VTK(int n_p, std::vector<double> vx, std::vector<double> vy,
    std::vector<double> vz, int n_t, std::vector<int> t0, std::vector<int> t1,
    std::vector<int> t2, std::vector<double> conlen, std::vector<int> mask,
    std::string filePath){

    std::ofstream outFile;
    outFile.open(filePath, std::ios::trunc);
    outFile << "# vtk DataFile Version 2.0\n";
    outFile << "comment\n";
    outFile << "ASCII\n\n";
    outFile << "DATASET UNSTRUCTURED_GRID\n";

    outFile << "POINTS"
          << std::setw(9)  << n_p << "  "
          << " double\n";

    for (int i=0; i<n_p;i++){
        outFile << std::setw(12) << std::setprecision(10) << vx[i]  << "  "
                << std::setw(12) << std::setprecision(10) << vy[i]  << "  "
                << std::setw(12) << std::setprecision(10) << vz[i]  << "\n";
    }

    outFile << "\n";

    outFile << "CELLS"
          << std::setw(12)  << n_t  << "  "
          << std::setw(12)  << n_t * 4  << "  " << "\n";
    for ( int i = 0; i < n_t; i++ )
    {
    outFile << std::setw(12) << "3"
            << std::setw(12) << t0[i] - 1  << " "
            << std::setw(12) << t1[i] - 1 << " "
            << std::setw(12) << t2[i] - 1 << "\n";
    }

    outFile << "\n";

    outFile << "CELL_TYPES"
          << std::setw(9)  << n_t  << "  " << "\n";
    for ( int i = 0; i < n_t; i++ )
    {
        outFile << std::setw(12) << "5" << "\n";
    }

    outFile << "CELL_DATA" << std::setw(9) << n_t;
    outFile << "SCALARS conlen float 1\n";
    outFile << "LOOKUP_TABLE default\n";

    for (int i =0; i < n_t; i++ )
    {
        outFile << std::setw(12) << conlen[i];
    }

    outFile << "SCALARS mask int 1\n";
    outFile << "LOOKUP_TABLE default\n";
    for (int i =0; i < n_t; i++ )
    {
        outFile << std::setw(12) << mask[i];
    }

    outFile.close();

}