#ifndef OF_2_SITFCL_FIELD_INTERNALS_H
#define OF_2_SITFCL_FIELD_INTERNALS_H

#include "fvCFD.H"
#include <vector>

namespace speeditcl
{
namespace io
{
namespace internal
{

std::vector<std::vector<double>> internalValues(const volScalarField& field);
std::vector<std::vector<double>> internalValues(const volVectorField& field);
std::vector<std::vector<double>> internalValues(const surfaceVectorField& field);
std::vector<std::vector<double>> internalValues(const surfaceScalarField& field);

} // namespace internal
} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_FIELD_INTERNALS_H
