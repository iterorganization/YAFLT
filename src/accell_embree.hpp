#ifndef ACCELL_EMBREE_H
#define ACCELL_EMBREE_H

// Make an #error only here, since the rest of the code requires this header
// anyway.
#if EMBREE_VERSION == 3
#include <embree3/rtcore.h>
#elif EMBREE_VERSION ==4
#include <embree4/rtcore.h>
#else
#error "EMBREE_VERSION not defined (properly) ! Should be either 3 or 4"
#endif

/// Class for loading shadowing geometries to Embree, to be used for RayCasting
/// when testing if a FL intersects the shadowing geometry. In essence the
/// setup for ray tracing is simple, for the purposes of Field Line Tracing and
/// other simple ray tracing examples. All that is required is the
/// loading/unloading of meshes (2D surface triangular meshes) and tests
/// whether a line or ray intersects with the geometry, which loaded geometry
/// and finally which triangle cell.
class EmbreeAccell
{
private:


    /// Embree device that is the default factory for Embree objects.
    RTCDevice m_device;
    /// Scene into which we commit geometries.
    RTCScene m_scene;

    /// The following booleans are used for tracking the state of RTCDevice and
    /// rtcScene.

    /// Boolean if device was created or not.
    bool m_device_created = false;
    /// Boolean if scene was created or not.
    bool m_scene_created = false;
    /// Boolean that says if Embree object is empty or filled with at least one
    /// shadowing geometry.
    bool m_isEmpty = true;


    /// Private variables for when we want to use Embree in Cython or other
    /// applications and do not want to have direct calls to Embree. Basically
    /// variables for wrapping Embree calls.
    RTCRayHit m_rayHit;
#if EMBREE_VERSION == 3
    RTCIntersectContext m_rayContext;
#endif

public:
    EmbreeAccell(bool initDevice=true);
    ~EmbreeAccell();

    /// Function that setups the RTCDevice.
    void createDevice();
    /// Function that deletes the RTCDevice, if it exists
    void deleteDevice();

    /// Function that creates the RTCScene if RTCDevice exists
    void createScene();
    /// Function that deletes the scene if RTCDevice and RTCScene exists.
    void deleteScene();

    /// Function for loading a triangular mesh into Embree. Float precision.
    /// If m_scene does not exist, the function returns the maximum value
    /// of unsigned int.
    /// @returns geom_id ID of the loaded geometry. Used for identification
    ///             or removal of the geometry from m_scene.
    /// @param[in] vertices 1D array of points (p1.x, p1.y, p1.z, p2.x, p2.y,...).
    /// @param[in] n_vertices Number of vertices (len(vertices)/3).
    /// @param[in] triangles 1D array containing all triangles (t1.p1, t1.p2,
    ///                      t1.p3, t2.p1, ...).
    /// @param[in] n_triangles Number of triangles (len(triangles)/3).
    /// @return The id of the loaded geometry. This is used to on order to 
    ///         understand which mesh or geometry is intersected, when 
    ///         intersection occur.
    unsigned int commitMesh(float* vertices, long int n_vertices,
                            unsigned* triangles, long int n_triangles);

    /// Deletes the loaded geometry from m_scene.
    /// @returns successful Returns True if the geometry was successfully
    ///                     removed from m_scene.
    /// @param[in] geom_id Id of the geometry that we wish to remove.
    bool deleteMesh(unsigned int geom_id);

#if EMBREE_VERSION == 3
    /// OpenMP friendly castRay function as the ray tracing structs are
    /// provided by the calling function
    ///
    /// @param[in, out] rayHit, contains the information of the ray and the hit
    /// @param[in] rayContext, contains other information relevant to ray tracing
    void castRay(RTCRayHit *rayHit, RTCIntersectContext *rayContext);
#elif EMBREE_VERSION == 4
    /// OpenMP friendly castRay function as the ray tracing structs are
    /// provided by the calling function
    ///
    /// @param[in, out] rayHit, contains the information of the ray and the hit
    void castRay(RTCRayHit *rayHit);
#endif

    /// Function for performing ray tracing. After calling this function you
    /// call checkIfHit to get confirmation if hit happened.
    /// @param[in] ox X component of origin point
    /// @param[in] oy Y component of origin point
    /// @param[in] oz Z component of origin point
    /// @param[in] dx X component of direction vector
    /// @param[in] dy Y component of direction vector
    /// @param[in] dz Z component of direction vector
    /// @param[in] tnear Starting length along the direction vector to
    ///                  check for intersection. 0 means that the ray starts at
    ///                  the origin point
    /// @param[in] tfar  Ending length along the direction vector to check
    ///                  for intersection
    void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                 double tnear, double tfar);

    /// Function that tells you if there is an intersection hit. Call this
    /// after caling castRay.
    ///  @param[out] hit A boolean, true if a geometry was hit, false
    ///                  otherwise.
    bool checkIfHit();
    /// This function returns the geometry ID if an intersection hit occured.
    /// @param[out] geomId An integer of the geometry that was hit. This
    ///                    integer corresponds to the integer of the geometry
    ///                    loaded with commitMesh
    int returnGeomId();
    /// This function returns the geometry ID if an intersection hit occured.
    /// @param[out] primId An integer of the geometry that was hit. This
    ///                    integer corresponds to the integer of the loaded
    ///                    geometry
    int returnPrimId();
};

#endif /*ACCELL_EMBREE_H*/