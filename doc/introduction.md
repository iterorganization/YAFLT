
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

For interpolating the input equilibrium data the implementation of 2D Bicubic
Spline interpolation method is used.

FLs are followed with the Runge-Kutta-Fehlberg Method (RKF45) method to ensure
good precision and stability when traveling along a FL on a magnetic surface.

Intersection tests on the FLs are performed by using the RayCast method from
the high-performance ray tracing kernel Embree (https://embree.org).

## Code

### src directory

Whole source code is situated in the **src/** directory.

The code that contains that calls external Embree for intersection tests

 - accell_embree.cpp
 - accell_embree.hpp

The code that implements the bicubic spline interpolation method:

 - bicubic.hpp
 - bicubic.cpp

The code that contains the RKF45 solver method (not actually used in the flt
code but used for testing the solver):

 - rkf45.cpp
 - rkf45.hpp

And finally the code that glues everything together via objects:

 - flt.cpp
 - flt.hpp

The code itself should be self-descriptive, therefore the best way is to dive
into the code to understand the inner workings.

### examples

The **examples/** directory contains small snippets of code for testing. There
is no order here except the most commonly used test to see if the code works.

### test

The **test/** directory contains the synthetic tests to see if the
implementations of the code behaves as expected.

#### test/interpolation

To test the implementation of the 2D Bicubic Spline interpolation method, this
directory contains code snippets where the implementation is pitted against the
ALGLIB++ code in terms of accuracy and performance.

In the directory there are also  tests to see if the code is written parallel
friendly. That is, no thread locks or races when using std::vector's as
containers for thread local data.

For instance what was observed that even if functions for the Bicubic Spline
interpolation method were written to be thread conscious, that is, that every
local thread variable is actually a std::vector and we use the thread id's as
an additional parameter, a slowdown was observed. The location of the code
where the slowdown occurs was not exactly located, but apparently reading is
not completely thread safe when using std::vector's or the type was not
correctly used. This resulted in supplying an instance of the 2D Bicubic Spline
object to each thread which finally yielded a speedup compared to ALGLIB++
implementation.


## The theory

To see which equations and logic are implemented look at the documentation of
the python interface of the code. In essence the magnetic fieldlines are
followed in an axysymmetric magnetic field in the cylindrical coordinate system
and its segments are tested if they intersect any geometry with Embree.

## Parallelism

The code is run in parallel mode in the main functions runFLT in the C++ code.
Users can control the number of CPUs to be used in parallel.

## Build

### Linux

In order to build L2G_cpp the following packages are required (sub-dependencies
are not written):

 - Embree >= 3
 - CMake >= 3.10

If embree is not installed via a system, but it is a custom/environment
installation then you have to be sure that CMake can find the FindX scripts
of Embree (by modifying the CMAKE_PREFIX_PATH).

```console
cd path/to/L2G_cpp
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=path/to/install
make -j8
make install
```

### Windows

In order to build on windows you are required to install the Microsoft MSVC
compiler, either part of visual studio or other ways.

As with Linux you require the following libraries:

 - Embree >= 3
 - CMake >= 3.10

Since this will be a static compilation, you can download the Embree library
and simply tell CMake where the CMake configuration files are located.

```console
cmake -A x64 -DCMAKE_PREFIX_PATH=c:\path\to\embree\cmake\files -DCMAKE_INSTALL_PREFIX=c:\path\to\installation -B"Path\To\Build\directory" -S"Path\to\source\directry"
```

Then run the following lines for building, cleaning and installing


```console
# To build
cmake --build c:\path\to\build --target ALL_BUILD --config Release # --parallel $CPUS
# To clean
cmake --build c:\path\to\build --target clean --config Release
# To install
cmake --build c:\path\to\build --target install --config Release
```
