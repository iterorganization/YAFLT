#ifndef ACCELL_EMBREE_H
#define ACCELL_EMBREE_H

#include <embree3/rtcore.h>
#include <vector>

/// Class for loading shadowing geometries to Embree, to be used for RayCasting
/// when testing if a FL intersects the shadowing geometry.
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


    // Private variables for when we do not want to use only Embree in Cython
    // or other applications and do not want to have direct calls to Embree.
    RTCRayHit m_rayHit;
    RTCIntersectContext m_rayContext;

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
    unsigned int commitMesh(float* vertices, long int n_vertices,
                            unsigned* triangles, long int n_triangles);

    /// Deletes the loaded geometry from m_scene.
    /// @returns successful Returns True if the geometry was successfully
    ///                     removed from m_scene.
    /// @param[in] geom_id Id of the geometry that we wish to remove.
    bool deleteMesh(unsigned int geom_id);

    /// OpenMP friendly castRay function as the ray tracing structs are
    /// provided by the calling function
    ///
    /// @param[in, out] rayHit, contains the information of the ray and the hit
    /// @param[in] rayContext, contains other information relevant to ray tracing
    void castRay(RTCRayHit *rayHit, RTCIntersectContext *rayContext);


    /// Functions for Cython exposure!
    void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                 double tnear, double tfar);
    bool checkIfHit();
    int returnGeomId();
    int returnPrimId();
};

#endif /*ACCELL_EMBREE_H*/