#ifndef OF_2_SITFCL_BOUNDARY_H
#define OF_2_SITFCL_BOUNDARY_H

#include "fvCFD.H"
#include <vector>

namespace speeditcl
{
namespace io
{

namespace internal
{

struct BFaceValues
{
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> types;
};

// WARNING - The values defined here MUST be identical
// with values from arael.h file
enum BoundaryFaceType
{
    FIXED_VALUE                = 0,
    CALCULATED                 = 1,
    FIXED_VALUE_USER_BEGIN     = 2,
    SLIP_BOUNDARY              = 3,

    NUT_WALL_FUNCTION          = 65,

    INLET_OUTLET               = 96,

    ZERO_GRADIENT              = 128,
    FIXED_GRADIENT             = 129,
    FIXED_GRADIENT_USER_BEGIN  = 130,
    OMEGA_WALL_FUNCTION        = 131,

    UNKNOWN                    = 255
};

std::vector<int> bfacesReorderTable(const fvMesh & mesh);

BoundaryFaceType bfaceType(const std::string& type_name);

BFaceValues boundaryValues(const volScalarField& field,
                           bool values_only = false);
BFaceValues boundaryValues(const volVectorField& field,
                           bool values_only = false);
BFaceValues boundaryValues(const surfaceScalarField& field,
                           bool values_only = false );
BFaceValues boundaryValues(const surfaceVectorField& field,
                           bool values_only = false);

} // namespace internal
} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_BOUNDARY_H
