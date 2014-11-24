#include "OF2SITFCLFieldInternals.h"

#include "OF2SITFCLUtils.h"

namespace speeditcl
{
namespace io
{
namespace internal
{

template<class T>
int numberOfElements(const Foam::Field<T>& field)
{
    return field.size();
}

template<template<class> class Field>
void initVec(std::vector<std::vector<double>>& v, const Field<scalar>& field)
{
    v.resize(1);
    v[0].resize(numberOfElements(field));
}

template<template<class> class Field>
void initVec(std::vector<std::vector<double>>& v, const Field<vector>& field)
{
    int nelements = numberOfElements(field);

    v.resize(3);
    v[0].resize(nelements);
    v[1].resize(nelements);
    v[2].resize(nelements);
}

template<class T>
std::vector<std::vector<double>> fieldValues (const Field<T>& field)
{
    std::vector< std::vector<double>> vres;
    initVec(vres, field);

    int nelements = numberOfElements(field);

    for(int i = 0; i < nelements; i++)
    {
        setVecVal(vres, field[i], i);
    }

    return vres;
}

template<class T, template<class> class TBoundary, class TMesh>
std::vector<std::vector<double>> internalValuesTempl(
                               const GeometricField<T, TBoundary, TMesh>& field)
{
    if(typeid(T) != typeid(scalar) && typeid(T) != typeid(vector))
    {
        throw std::string("Unsupported field data type");
    }

    return fieldValues(field.internalField());
}


std::vector<std::vector<double>> internalValues(const volScalarField& field)
{
    return internalValuesTempl(field);
}

std::vector<std::vector<double>> internalValues(const volVectorField& field)
{
    return internalValuesTempl(field);
}

std::vector<std::vector<double>> internalValues(const surfaceVectorField& field)
{
    return internalValuesTempl(field);
}

std::vector<std::vector<double>> internalValues(const surfaceScalarField& field)
{
    return internalValuesTempl(field);
}

} // namespace internal
} // namespace io
} // namespace speeditcl
