#include <accell_embree.hpp>

#ifndef NDEBUG
#include <iostream>
#endif

// We will register this error handler with the device in createDevice(),
// so that we are automatically informed on errors.
// This is extremely helpful for finding bugs in your code, prevents you
// from having to add explicit error checking to each Embree API call.
void error_handler(void *userPtr, const RTCError code, const char *str){
#ifndef NDEBUG
    std::cout << "error "<< code << ": " << str << std::endl;
#endif
}

EmbreeAccell::EmbreeAccell(bool initialize){
    if (initialize) {
        createDevice();
        createScene();
    }

    /// The following variables are used when we want to run in serial embree
    /// functions and also used via Python for other purposes.
    m_rayHit = RTCRayHit();
    m_rayContext = RTCIntersectContext();
    rtcInitIntersectContext(&m_rayContext);
}

EmbreeAccell::~EmbreeAccell(){
    deleteDevice();
}

// Initializes Embree factory and essential objects
void EmbreeAccell::createDevice(){
    if (m_device_created){
#ifndef NDEBUG
        std::cout << "Error: device already created.";
#endif
        return;
    }

    RTCDevice device = rtcNewDevice(NULL);

    if (!device) {
#ifndef NDEBUG
        std::cout << "Error " << rtcGetDeviceError(NULL);
        std::cout << ": cannot create device" << std::endl;
#endif
        return;
    }

    // Register the error function.
    rtcSetDeviceErrorFunction(device, error_handler, NULL);
    m_device = device;
    m_device_created = true;
}

void EmbreeAccell::deleteDevice(){
    // Deletes rtcDevice, and subsequently deletes the scene.

    if (!m_device_created){
#ifndef NDEBUG
        std::cout << "Error: no device to delete!";
#endif
        return;
    }

    // First check if scene is created
    if (m_scene_created){
        // Delete the scene
        deleteScene();
    }

    rtcReleaseDevice(m_device);
    m_device_created = false;
}

void EmbreeAccell::createScene(){
    if (!m_device_created){
        // No device created, ergo we cannot create a scene.
#ifndef NDEBUG
        std::cout << "Error: cannot create a scene without a device.";
#endif
        return;
    }

    if (m_scene_created){
        // Scene already exists
#ifndef NDEBUG
        std::cout << "Error: Scene already exists!";
#endif
        return;
    }

    m_scene = rtcNewScene(m_device);

    // Set Scene flags to best quality settings.
    rtcSetSceneFlags(m_scene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(m_scene, RTC_BUILD_QUALITY_HIGH);
    m_scene_created = true;
}

void EmbreeAccell::deleteScene(){
    if (!m_scene_created){
#ifndef NDEBUG
        std::cout << "Error: No scene to delete.";
#endif
        return;
    }

    rtcReleaseScene(m_scene);
    m_scene_created = false;
}

unsigned int EmbreeAccell::commitMesh(float* vertices, long int n_vertices,
                              unsigned* triangles, long int n_triangles){

    if (!m_scene_created){
#ifndef NDEBUG
        std::cout << "Error: No scene created to which to commit the geometry!";
        return -1;
#endif
    }

    // Commits a geometry, in this case a triangular mesh to the Embree scene.
    unsigned int geom_id;

    RTCGeometry geom = rtcNewGeometry(m_device, RTC_GEOMETRY_TYPE_TRIANGLE);

    // Construct structs for adding vertices (points) and triangles
    // (contaning vertex ids)
    float *rtc_vertices = (float*) rtcSetNewGeometryBuffer(geom,
                RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
                3*sizeof(float), n_vertices);
    unsigned* rtc_indices = (unsigned*) rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3*sizeof(unsigned),
        n_triangles);

    if (rtc_vertices && rtc_indices)
    {
#ifndef NDEBUG
        std::cout << "Populating embree object" << std::endl;
#endif
        for (long int i=0; i < n_vertices; i++){
            rtc_vertices[i] = vertices[i];
        }
        for (long int i=0; i < n_triangles; i++){
            rtc_indices[3 * i] = triangles[3 * i];
            rtc_indices[3 * i + 1] = triangles[3 * i + 1];
            rtc_indices[3 * i + 2] = triangles[3 * i + 2];
        }
    }
#ifndef NDEBUG
    std::cout << "Added " << n_vertices << " vertices." << std::endl;
    std::cout << "Added " << n_triangles << " triangles." << std::endl;
#endif
    // You must commit geometry objects when you are done setting them up,
    // or you will not get any intersections.
    rtcCommitGeometry(geom);

    // In rtcAttachGeometry(...), the scene takes ownership of the geom
    // by increasing its reference count. This means that we don't have
    // to hold on to the geom handle, and may release it. The geom object
    // will be released automatically when the scene is destroyed.
    //
    // rtcAttachGeometry() returns a geometry ID. We could use this to
    // identify intersected objects later on.

    geom_id = rtcAttachGeometry(m_scene, geom);
    rtcReleaseGeometry(geom);

    // Like geometry objects, scenes must be committed. This lets
    // Embree know that it may start building an acceleration structure.
    rtcCommitScene(m_scene);
    m_isEmpty = false;
    return geom_id;
}

bool EmbreeAccell::deleteMesh(unsigned int geom_id){
    // Deletes a geometry with the assigned geometry id geom_id. If it fails
    // to delete it, simply because either the geometry with such id does not
    // exist or of some internal error, the function retuns false.
    if (!m_scene_created){
#ifndef NDEBUG
        std::cout << "Scene doesn't exist, from which to delete the mesh";
#endif
        return false;
    }

    rtcDetachGeometry(m_scene, geom_id);
    if (rtcGetDeviceError(m_device) != RTC_ERROR_NONE) {
        /*Failed to remove mesh*/
        return false;
    }

    rtcCommitScene(m_scene);
    return true;
}

void EmbreeAccell::castRay(RTCRayHit *ray, RTCIntersectContext *context){
    rtcIntersect1(m_scene, context, ray);
}

// Functions for Cython (so no Embree is needed to be included)

void EmbreeAccell::castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                           double tnear, double tfar){
    m_rayHit.ray.org_x = ox;
    m_rayHit.ray.org_y = oy;
    m_rayHit.ray.org_z = oz;
    m_rayHit.ray.dir_x = dx;
    m_rayHit.ray.dir_y = dy;
    m_rayHit.ray.dir_z = dz;
    m_rayHit.ray.tnear = tnear;
    m_rayHit.ray.tfar = tfar;
    m_rayHit.ray.mask = 0;
    m_rayHit.ray.flags = 0;
    m_rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    m_rayHit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rtcIntersect1(m_scene, &m_rayContext, &m_rayHit);
}

bool EmbreeAccell::checkIfHit(){
    if (m_rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID){
        return true;
    }
    return false;
}

int EmbreeAccell::returnGeomId(){
    return (int) m_rayHit.hit.geomID;
}

int EmbreeAccell::returnPrimId(){
    /// It makes no sense why the geomID in Embree documentation is called the
    /// primitive ID and the primID is called the geometry ID... Maybe typo.
    return (int) m_rayHit.hit.primID;
}
