
# About

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


## Building the software

In order to build the software the following dependencies are required:

 - Embre >= 3
 - CMake >= 3.12
 - OpenMP support
 - Doxygen (optional)
 - Ninja (optional)

Linux:

 - GCC
 - build tools

Windows:

 - Visual Studio 17 2022
 - Visual Studio CMake support


Embree (https://github.com/RenderKit/embree) is a high-performance CPU ray 
tracing kernel used for intersection testing. Both Embree 3 and 4 are 
supported.

Doxygen is optional. It is used for generating html version of the 
documentation.

Ninja is optional. It is used for building the software.

### Building with CMake

Ensure that Embree and other dependencies can be found with CMake. Otherwise
edit the CMAKE_PREFIX_PATH in order to point to the CMake directories of the
dependencies.

```bash
# Clone or extract the code
mkdir build

# Generate the build system
cmake -S <path-to-source> -B <path-to-build> -GNinja -DCMAKE_INSTALL_PREFIX=<path-to-install> -DCMAKE_BUILD_TYPE=Release

# Build the project
cmake --build <path-to-build>

# Test the build
cmake -E chdir <path-to-build> ctest

# Install the project
cmake --install <path-to-build>
```

## Debug build CMake for IDE editor

In case you are using an editor with LSP plugins or functionality and you have
installed the LSP-clangd plugin, the following steps has to be performed in
order to have the LSP-clangd enhanced experience.

### Generate the compile build JSON file

When compiling pass the the flag ``-DCMAKE_EXPORT_COMPILE_COMMANDS=ON``. This
will generate a JSON file that gives the LSP-clangd compilation information

```bash
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug <path-to-source>
```

Now from the build directory copy the compile_commands.json to the root of the
source directory.

## Usage of the code in other CMake projects

The software can be imported in other CMake projects. The CMake directory has
to be appended to the CMAKE_PREFIX_PATH:

```bash
export CMAKE_PREFIX_PATH=<path-to-installed-cmake-files>:${CMAKE_PREFIX_PATH}
```

and then can be used in CMake as following:

```cmake
find_package(FLT REQUIRED)
find_package(OpenMP)

add_executable(<TARGET> <TARGET_CONFIGURATION>)
target_link_libraries(<TARGET> flt)
target_include_directories(<TARGET> PUBLIC ${FLT_INCLUDE_DIRS})
target_compile_definitions(<TARGET> PUBLIC ${FLT_COMPILE_DEFINITIONS})
```