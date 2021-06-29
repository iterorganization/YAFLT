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

    /// List of RTCIntersectContext, a required struct when using RT functions.
    std::vector<RTCIntersectContext> m_listContext;

    /// The following booleans are used for tracking the state of RTCDevice and
    /// rtcScene.

    /// Boolean if device was created or not.
    bool m_device_created = false;
    /// Boolean if scene was created or not.
    bool m_scene_created = false;
    /// Boolean that says if Embree object is empty or filled with at least one
    /// shadowing geometry.
    bool m_isEmpty = true;

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

    /// Cast ray function for determining if there is an intersection with the
    /// shadowing geometry on a finite ray (segment of FL).
    /// By providing the origin point, direction of the FL and the length of the
    /// segment Embree ray casts it to the loaded shadowing geometries.
    ///
    /// Note that the segment is straight, i.e., not a curve, therefore when
    /// following a FL one must be sure that the segment is short enough so that
    /// there is no precision loss.
    ///
    /// @param[in] ox X coordinate of the origin point.
    /// @param[in] oy Y coordinate of the origin point.
    /// @param[in] oz Z coordinate of the origin point.
    /// @param[in] dx X component of the direction vector (normalized).
    /// @param[in] dy Y component of the direction vector (normalized).
    /// @param[in] dz Z component of the direction vector (normalized).
    /// @param[in] tnear Basically it tells embree from where on the segment
    ///                  it should start checking. Should always be 0.
    /// @param[in] tfar Length of the FL segment. With this you specify the
    ///                 length of your segment. Usually it is the length
    ///                 between the current point and next point that your
    ///                 solver outputs.
    /// @param[in] omp_thread ID of the OpenMP thread. Specifies the location of
    ///                       thread related variables. Users should only take
    ///                       into consideration when parallel OpenMP blocks are
    ///                       used. Otherwise the default values for sequential
    ///                       run are adequate.
    void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                 double tnear, double tfar, int omp_thread=0);

    /// After calling castRay, checks if there is a hit with the shadowing
    /// geometry that is loaded to Embree.

    /// @param[in] omp_thread ID of the OpenMP thread. Specifies the location of
    ///                       thread related variables. Users should only take
    ///                       into consideration when parallel OpenMP blocks are
    ///                       used. Otherwise the default values for sequential
    ///                       runs are adequate.
    bool checkIfHit(int omp_thread=0);

    /// Function that prepares local thread related variables. In order to avoid
    /// writing parallel code inside the low-level kernel, one must prepare the
    /// variables in such a way when multiple threads are called and are
    /// calling this class functions, that when they write the data, they write
    /// to it's own local variables. Hence for example if in your application
    /// you are calling a omp pragma with 6 threads, then before the parallel
    /// block this function should be called with the argument equaled to at
    /// least 6 so that each thread has their own thread variables. Calling it
    /// with a higher value than 6 will still achieve that each thread has their
    /// own thread variables, but at the same time you have leftovers and thus
    /// a bit higher memory usage.
    ///
    /// Also in order to avoid code duplication, even when running sequentially
    /// or on "one thread", this function must be called, without argument, so
    /// that the default value is 1 or with value 1. Value should be nonnegative
    /// and higher than 0!!
    ///
    /// @param[in] num_threads Number of threads for which local thread
    ///                        variables should be prepared.
    void prepareThreadContainers(int num_threads=1);
    /// Struct containing RayHit information
    /// To be used as getting information from the intersection test.
    /// For future maybe we need std::list<RTCRayHit> for parallel runs.

    /// A vector of RTCRayHit structures, used to store information about the
    /// FL segment when calling Embree ray cast function.
    std::vector<RTCRayHit> m_listRayHit;

};

#endif /*ACCELL_EMBREE_H*/