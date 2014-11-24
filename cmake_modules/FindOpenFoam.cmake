SET(OPEN_FOAM_INCLUDE_DIR
        $ENV{FOAM_SRC}
        $ENV{FOAM_SRC}/finiteVolume/lnInclude
        $ENV{FOAM_SRC}/turbulenceModels
        $ENV{FOAM_SRC}/turbulenceModels/incompressible/RAS/lnInclude
        $ENV{FOAM_SRC}/turbulenceModels/incompressible/turbulenceModel/lnInclude
        $ENV{FOAM_SRC}/transportModels/incompressible/lnInclude
        $ENV{FOAM_SRC}/transportModels
        $ENV{FOAM_SRC}/OpenFOAM/lnInclude
        $ENV{FOAM_SRC}/meshTools/lnInclude
        $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
        $ENV{FOAM_SRC}/OSspecific/POSIX/cpuTime)

SET(OPEN_FOAM_LIB_SEARCH_PATH
        ${OPEN_FOAM_LIB_SEARCH_PATH}
        $ENV{FOAM_LIBBIN}
)

message( "Search Path: ${OPEN_FOAM_LIB_SEARCH_PATH}" )

FIND_LIBRARY(
    OPEN_FOAM_LIBRARY_OF
    NAMES OpenFOAM
    PATHS ${OPEN_FOAM_LIB_SEARCH_PATH}
)

FIND_LIBRARY(
    OPEN_FOAM_LIBRARY_FINITE_VOLUME
    NAMES finiteVolume
    PATHS ${OPEN_FOAM_LIB_SEARCH_PATH}
)

FIND_LIBRARY(
    OPEN_FOAM_LIBRARY_PSTREAM
    NAMES Pstream
    PATHS ${OPEN_FOAM_LIB_SEARCH_PATH}/openmpi-system ${OPEN_FOAM_LIB_SEARCH_PATH}/openmpi-1.5.3
)

set(OPEN_FOAM_LIBRARY)

if(OPEN_FOAM_LIBRARY_OF)
    set(OPEN_FOAM_LIBRARY ${OPEN_FOAM_LIBRARY} ${OPEN_FOAM_LIBRARY_OF})
endif(OPEN_FOAM_LIBRARY_OF)

if(OPEN_FOAM_LIBRARY_FINITE_VOLUME)
    set(OPEN_FOAM_LIBRARY ${OPEN_FOAM_LIBRARY} ${OPEN_FOAM_LIBRARY_FINITE_VOLUME})
endif(OPEN_FOAM_LIBRARY_FINITE_VOLUME)

if(OPEN_FOAM_LIBRARY_PSTREAM)
    set(OPEN_FOAM_LIBRARY ${OPEN_FOAM_LIBRARY} ${OPEN_FOAM_LIBRARY_PSTREAM})
endif(OPEN_FOAM_LIBRARY_PSTREAM)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  OPEN_FOAM
  DEFAULT_MSG
  OPEN_FOAM_LIBRARY
  OPEN_FOAM_INCLUDE_DIR
)

if(OPEN_FOAM_FOUND)
  set(OPEN_FOAM_LIBRARIES ${OPEN_FOAM_LIBRARY})
else(OPEN_FOAM_FOUND)
  set(OPEN_FOAM_LIBRARIES)
endif(OPEN_FOAM_FOUND)

mark_as_advanced(
  OPEN_FOAM_INCLUDE_DIR
  OPEN_FOAM_LIBRARY
)