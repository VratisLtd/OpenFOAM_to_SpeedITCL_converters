project(converters_common)

set(libname converters_common)

include_directories(include)
include_directories(${OPEN_FOAM_INCLUDE_DIR})

set(lib_src
    src/OF2SITFCLFieldWriter.cpp
    src/OF2SITFCLMeshWriter.cpp
    src/OF2SITFCLBoundary.cpp
    src/OF2SITFCLFieldInternals.cpp
    src/OF2SITFCLUtils.cpp
    src/OF2SITFCLErrorHandling.cpp
    src/OF2SITFCLSettings.cpp
    src/OF2SITFCLSolverLoader.cpp
)

set(lib_headers
    include/OF2SITFCLVectorWriter.h
    include/OF2SITFCLFieldWriter.h
    include/OF2SITFCLMeshWriter.h
    include/OF2SITFCLBoundary.h
    include/OF2SITFCLFieldInternals.h
    include/OF2SITFCLUtils.h
    include/OF2SITFCLConstants.h
    include/OF2SITFCLErrorHandling.h
    include/OF2SITFCLSettings.h
    include/OF2SITFCLSolverLoader.h
)

message ( STATUS "open foam library: ${OPEN_FOAM_LIBRARY}")

add_library("${libname}" STATIC ${lib_src} ${lib_headers})
target_link_libraries("${libname}" ${OPEN_FOAM_LIBRARY}
                                    incompressibleRASModels
                                    incompressibleTransportModels
                                    incompressibleTurbulenceModel
                                    OpenFOAM
                                    Pstream)
