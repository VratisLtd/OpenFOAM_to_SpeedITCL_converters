#include "OF2SITFCLBoundary.h"

#include "OF2SITFCLUtils.h"

namespace speeditcl
{
namespace io
{
namespace internal
{

std::vector<int> bfacesReorderTable(const fvMesh & mesh)
{
    struct BFacePos
    {
        int pos, owner;
    };

    std::vector< BFacePos > vbp ;
    int pos = 0;
    auto& boundary = mesh.boundary();

    for(int z = 0; z < boundary.size(); z++)
    {
        for (int i = 0; i < boundary[z].size() ; i++)
        {
            BFacePos bp;
            bp.pos = pos;
            bp.owner = boundary[z].faceCells()[i];
            pos++;
            vbp.push_back(bp);
        }
    }

    std::stable_sort(vbp.begin(), vbp.end(),
        [](const BFacePos& a, const BFacePos& b) { return a.owner < b.owner; });

    std::vector<int> rt (vbp.size());

    for (unsigned i = 0; i < rt.size(); i++)
    {
            rt[i] = vbp[i].pos ;
    }

    return rt;
}

BoundaryFaceType bfaceType(const std::string& type_name)
{
    if( "zeroGradient" == type_name ||
        "kqRWallFunction" == type_name)
    {
        return ZERO_GRADIENT;
    }
    else if("fixedValue" == type_name
            || "timeVaryingUniformFixedValue" == type_name
            || "surfaceNormalFixedValue" == type_name)
    {
        return FIXED_VALUE;
    }
    else if("omegaWallFunction" == type_name)
    {
        return OMEGA_WALL_FUNCTION;
    }

    else if("nutkWallFunction" == type_name)
    {
        return NUT_WALL_FUNCTION;
    }

    else if ("slip" == type_name)
    {
        return SLIP_BOUNDARY;
    }

    else if ("inletOutlet" == type_name)
    {
        return INLET_OUTLET;
    }

    else if ("calculated" == type_name)
    {
        return CALCULATED;
    }

    else if("empty" == type_name)
    {
        return UNKNOWN ; // TODO: Maybe should be special type "EMPTY" ?
    }
    else
    {
        throw std::string("Unsupported boundary face type " + type_name
                          + ". The only supported boundary face types are\n"
                          + "zeroGradient \n" +
                          + "fixedValue \n" +
                          + "inletOutlet \n" +
                          + "slip \n" +
                          + "omegaWallFunction\n" +
                          + "kqRWallFunction\n" +
                          + "nutkWallFunction\n" +
                          + "calculated.");
    }
}

template<template<class> class Field, class Type>
int numberOfElements(FieldField<Field, Type> const & field)
{
        int n = 0;

        for (int z=0 ; z<field.size() ; z++)
        {
                n += field[z].size();
        }

        return n;
}

template<template<class> class Field>
void initVec(std::vector< std::vector<double>> & vec,
             const FieldField<Field, scalar>& field)
{
    int nbfaces = numberOfElements(field);

    vec.resize(1);
    vec[0].resize(nbfaces);
}

template<template<class> class Field>
void initVec(std::vector< std::vector<double>>& vec,
             const FieldField<Field, vector>& field)
{
    int nbfaces = numberOfElements(field);

    vec.resize(3);
    vec[0].resize(nbfaces);
    vec[1].resize(nbfaces);
    vec[2].resize(nbfaces);
}

template<template<class> class Field, class Type>
std::vector< std::vector<double>> fieldValues(
                                      const FieldField<Field, Type>& field)
{
    std::vector<std::vector<double>> vres;
    initVec(vres, field);

    int idx = 0;
    for(int z = 0; z < field.size(); z++)
    {
        for (int i = 0; i < field[z].size(); i++)
        {
            setVecVal(vres, field[z][i], idx);
            idx++;
        }
    }

    return vres;
}

template<template<class> class Field, class Type>
std::vector<std::vector<double>> fieldValuesReordered(
                                     const FieldField<Field, Type>& field,
                                     std::vector<int>& reorder_table)
{
    return vector_reorder(fieldValues(field), reorder_table);
}

template<class TField>
BFaceValues boundaryValues(const TField& field, bool values_only)
{
    if(typeid(TField) != typeid(volScalarField)
            && typeid(TField) != typeid(volVectorField)
            && typeid(TField) != typeid( surfaceScalarField)
            && typeid(TField) != typeid( surfaceVectorField))
    {
        throw std::string("Unsupported boundary field data type");
    }

    std::vector<int> vbp = bfacesReorderTable(field.mesh());

    int nbfaces = vbp.size();

    BFaceValues result;

    result.vals = fieldValuesReordered(field.boundaryField(), vbp);

    if (false == values_only)
    {
        std::vector<int> types(nbfaces);

        int idx = 0;
        for (int z = 0; z < field.boundaryField().size(); z++)
        {
            BoundaryFaceType ft = bfaceType(
                        std::string(field.boundaryField().types()[z]));

            for (int i = 0; i < field.boundaryField()[z].size(); i++)
            {
                types[idx] = ft;
                idx++;
            }
        }

        result.types.resize(1);
        result.types[0] = vectorReorder(types, vbp);
    }

    return result;
}

BFaceValues boundaryValues(const volScalarField & field, bool values_only)
{
    return boundaryValues<volScalarField>(field, values_only);
}

BFaceValues boundaryValues(const volVectorField & field, bool values_only)
{
    return boundaryValues<volVectorField>(field, values_only);
}

BFaceValues boundaryValues(const surfaceScalarField& field, bool values_only)
{
    return boundaryValues<surfaceScalarField>(field, values_only);
}

BFaceValues boundaryValues(const surfaceVectorField& field, bool values_only)
{
    return boundaryValues<surfaceVectorField>(field, values_only);
}

} // namespace internal
} // namespace io
} // namespace speeditcl
