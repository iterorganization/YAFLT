
Description {#mainpage}
=======================

## Introduction

L2G_cpp is a C++ library used for running FLT runs in a tokamak-type geometry.
FLT is performed by following FLs from a desired mesh, in an axisymmetric
magnetic field provided in the form of the poloidal magnetic flux and the
poloidal current function and performing intersecting tests on neighboring or
shadowing geometry.

The equilibrium data is provided to the library via function calls and is not
limited by I/O, meaning that an application can get it's equilibrium data
from an arbitrary source and feed it to the library via function calls.

For interpolating the input equilibrium data the ALGLIB
(https://www.alglib.net/) library is used (2D BSpline interpolation method).

FLs are followed with the Runge-Kutta-Fehlberg Method (RKF45) method to ensure
good precision and stability when traveling along a FL on a magnetic surface.

Intersection tests on the FLs are performed by using the RayCast method from
the high-performance ray tracing kernel Embree (https://embree.org).

## Build

In order to build L2G_cpp the following packages are required:

 - Embree >= 3.12.2
 - CMake >= 3.12.1

If embree is not installed via a system, but it is a custom/environment
installation then you have to be sure that CMake can find the FindX scripts
of Embree.

```console
cd path/to/L2G_cpp
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=path/to/install
make -j8
make install
```