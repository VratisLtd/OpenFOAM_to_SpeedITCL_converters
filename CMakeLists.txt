project(SpeedITFlowConverters)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)

add_definitions(
    -std=c++11 # support for new cpp standard
)

## OpenFoam
find_package(OpenFoam REQUIRED)

include_directories(common/include)
include_directories(${OPEN_FOAM_INCLUDE_DIR})

## Boost
#set(Boost_NO_SYSTEM_PATHS ON)
set(Boost_USE_STATIC_LIBS ON)

find_package(Boost COMPONENTS program_options REQUIRED)

include_directories(SYSTEM ${Boost_INCLUDE_DIR})


add_definitions(
    -DWM_DP -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor
    -DNoRepository -ftemplate-depth-100
)

#string(SUBSTRING $ENV{WM_PROJECT_VERSION} 0 1 OF_VERSION_TRIMMED)

string(SUBSTRING "2.0.1" 0 1 OF_VERSION_TRIMMED)
message(STATUS "OPEN FOAM VERSION: ${OF_VERSION_TRIMMED}")

add_definitions(-DOF_VERSION=${OF_VERSION_TRIMMED})

set(CONVERTERS
    OpenFoam2SpeedITCL
)

foreach(converter ${CONVERTERS})
        add_executable(${converter} ${converter})
        target_link_libraries(${converter} converters_common)
        target_link_libraries(${converter} incompressibleRASModels
                                    incompressibleTransportModels
                                    incompressibleTurbulenceModel)

        target_link_libraries(${converter} ${Boost_PROGRAM_OPTIONS_LIBRARY})
endforeach(converter)

message(STATUS "boost program options lib: ${Boost_PROGRAM_OPTIONS_LIBRARY}")

add_subdirectory(common)
