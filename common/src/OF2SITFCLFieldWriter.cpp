#include "OF2SITFCLFieldWriter.h"

#include <string>

#include "OF2SITFCLBoundary.h"
#include "OF2SITFCLFieldInternals.h"
#include "OF2SITFCLUtils.h"
#include "OF2SITFCLErrorHandling.h"

namespace speeditcl
{
namespace io
{
namespace internal
{

} // namespace internal


template<class TField>
void OF2SpeeditclFieldWriter::saveField(TField& field, const Foam::fvMesh& mesh,
                                    const std::string& dir_path,
                                    Foam::Time* runTime, bool save_data)
{
    std::string field_name(field.name());
    std::string filename = dir_path + "/" + field_name;

    try
    {
        std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

        if (!file.is_open())
            throw "Can not open file " + filename;

        file << FILE_HEADER  << "\n";
        file << FIELD_HEADER << "\n" << field_name << "\n";

        if(typeid(TField) == typeid(volScalarField)
                || typeid(TField) == typeid(surfaceScalarField )
                || typeid(TField) == typeid(volVectorField)
                || typeid(TField) == typeid(wallDist) )
        {
        }
        else
        {
            throw std::string("Unsupported field type");
        }

        if(save_data)
        {
            vectorWriter.saveVec(file, internal::internalValues(field),
                                 CELL_VALS_HEADER);

            internal::BFaceValues bfv = internal::boundaryValues(field);

            vectorWriter.saveVec(file, bfv.vals , BFACE_VALS_HEADER );
            vectorWriter.saveVec(file, bfv.types, BFACE_TYPES_HEADER);
        }

        {
            OStringStream hss;

            field.writeHeader(hss);
            hss.writeKeyword("dimensions")
                    << field.dimensionedInternalField().dimensions()
                    << token::END_STATEMENT << nl << nl;

            if (hss.bad())
                throw std::string("Can not generate OpenFOAM file header");

            file << "\n" << FIELD_OF_HEADER << "\n";
            file << hss.str().size() << "\n";
            file << hss.str();
        }

        if(save_data)
        {
            OStringStream bss;

            field.boundaryField().writeEntry("boundaryField", bss);
            field.writeEndDivider(bss);

            if (bss.bad())
            throw std::string("Can not generate OpenFOAM file footer");

            file << "\n" << FIELD_OF_FOOTER << "\n";
            file << bss.str().size() << "\n";
            file << bss.str();
        }
        else // Used only for "phi" field
        {
            file << "\n" << FIELD_OF_BOUNDARY_PATCHES << "\n";
            file << field.boundaryField().size() << "\n";

            auto& bfield = field.boundaryField();

            for(int i = 0; i < bfield.size(); i++)
            {
                if(("empty" != bfield.types()[i])
                        && ("calculated" != bfield.types()[i]))
                    throw string("Unsupported boundary patch type \"")
                        + bfield.types()[i] + "\", the only supported types "
                        + "are \"empty\" \"calculated\"" ;

                file << bfield[i].patch().name() << " ";
                file << bfield.types()[i] << " ";
                file << bfield[i].patch().size();
                file << "\n";
            }
        }

        internal::testIO(file);

        if (runTime != nullptr)
        {
            saveTimeVaryingBC(field, mesh, dir_path, *runTime);
            saveGroovyBC(field, mesh, dir_path, *runTime);
            saveTotalPressureBC(field, mesh, dir_path, *runTime);
        }

    }
    catch (std::string & err_str)
    {
        handleError("save_field ERROR when saving " + field.name()
                     + " : " + err_str, dir_path);
    }
}

void OF2SpeeditclFieldWriter::saveField(volScalarField& field,
                                    const Foam::fvMesh& mesh,
                                    const std::string& dir_path,
                                    Foam::Time* runTime)
{
    saveField<volScalarField>(field, mesh, dir_path, runTime);
}

void OF2SpeeditclFieldWriter::saveField(volVectorField& field,
                                    const Foam::fvMesh& mesh,
                                    const std::string& dir_path,
                                    Foam::Time* runTime)
{
    saveField<volVectorField>(field, mesh, dir_path, runTime);
}

void OF2SpeeditclFieldWriter::saveField(surfaceScalarField & field,
                                    const Foam::fvMesh& mesh,
                                    const std::string& dir_path,
                                    Foam::Time* runTime)
{
    /** used only for field "phi", where data is not needed */
    saveField<surfaceScalarField>(field, mesh, dir_path, runTime, false);
}

void OF2SpeeditclFieldWriter::saveField(wallDist& field,
                                    const Foam::fvMesh& mesh,
                                    const std::string& dir_path,
                                    Foam::Time* runTime)
{
    saveField<wallDist>(field, mesh, dir_path, runTime);
}

bool OF2SpeeditclFieldWriter::isTimeVarying(const std::string& type)
{
    if("timeVaryingUniformFixedValue" == type)
    {
        return true;
    }
    return false;
}

bool OF2SpeeditclFieldWriter::isGroovy(const std::string& type)
{
    if("groovyBC" == type)
    {
        return true;
    }
    return false;
}

bool OF2SpeeditclFieldWriter::isTotalPressure(const std::string& type)
{
    return "totalPressure" == type;
}

template<class TField>
void OF2SpeeditclFieldWriter::buildTimeVaryingBoundaryPatch(
    TimeVaryingBoundaryPatch& patch, TField& field, const Foam::fvMesh& mesh,
    Foam::Time& runTime, int bfpatch_idx)
{
    patch.initialized = false;

    if(!isTimeVarying(field.boundaryField().types()[bfpatch_idx]))
        throw std::string("Can not build time varying boundary"
                          " patch from OpenFOAM patch of type \"")
            + field.boundaryField().types()[bfpatch_idx]
            + "\", possibly wrong patch index passed";

    if((0 > bfpatch_idx) || (field.boundaryField().size() <= bfpatch_idx))
        throw std::string("Invalid patch index passed") ;

    buildTimeVaryingPatchFaceIndices(patch, field, mesh, bfpatch_idx);
    buildTimeVaryingPatchFaceValues(patch, field, mesh, runTime, bfpatch_idx);

    patch.initialized = true ;
}

template<class TField>
void OF2SpeeditclFieldWriter::buildTimeVaryingPatchFaceIndices(
    TimeVaryingBoundaryPatch& patch, const TField& field,
    const Foam::fvMesh& mesh, int bfpatch_idx)
{
    // Obtain indices into the aggregated, unsorted
    // array containing consecutive OpenFOAM patches.
    std::vector<int> unsorted_face_positions;
    int idx = 0;
    // Find start position (skip all previous patches)
    for (int i = 0; i < bfpatch_idx; i++)
    {
        idx += field.boundaryField()[i].size();
    }

    // Add all faces from current patch
    for (int i = 0; i < field.boundaryField()[bfpatch_idx].size(); i++)
    {
        unsorted_face_positions.push_back(idx + i);
    }

    std::vector<int> bf_reorder_table = internal::bfacesReorderTable(mesh) ;

    std::vector<int> bf_reorder_table_2;
    bf_reorder_table_2.assign( bf_reorder_table.size(), -1);

    for(unsigned i = 0; i<bf_reorder_table.size(); i++)
    {
        bf_reorder_table_2[bf_reorder_table[i]] = i;
    }

    for (unsigned i = 0; i < unsorted_face_positions.size(); i++)
    {
        patch.face_positions.push_back(
                    bf_reorder_table_2[unsorted_face_positions[i]]);
    }
}


template<class TField>
void OF2SpeeditclFieldWriter::buildTimeVaryingPatchFaceValues(
    TimeVaryingBoundaryPatch& patch, TField& field, const Foam::fvMesh& mesh,
    Foam::Time& runTime, int bfpatch_idx)
{
    int n_time_steps = internal::computeNTimeSteps(runTime);
//runTime.deltaT()
    patch.face_values.resize(n_time_steps);
    patch.timeline.resize(n_time_steps);

    dimensionedScalar start_time = runTime.startTime();
    label start_time_index = runTime.startTimeIndex();

    field.correctBoundaryConditions();

    unsigned time_idx = 0;

    double timePointValue = runTime.startTime().value();

    while (runTime.loop())
    {
        patch.timeline[time_idx] = timePointValue;
        timePointValue += runTime.deltaT0().value();

        std::cout << "<<<<<<<<< Saving time step " << time_idx << std::endl;

        field.correctBoundaryConditions() ;

        patch.face_values[time_idx].resize(n_components(field));

        int nfaces_in_patch = field.boundaryField()[bfpatch_idx].size() ;

        Info << gMax(field.boundaryField()[bfpatch_idx]) << endl;

        // Each component is treated separately
        for (int c = 0; c < n_components(field); c++)
        {
            // Test, if all values are the same.
            // Allows for optimisation of memory usage.
            bool identical = true;
            for(int fidx = 0; fidx < nfaces_in_patch; fidx++)
            {
                if (getComponent(field.boundaryField()[bfpatch_idx][fidx], c)
                        != getComponent(field.boundaryField()[bfpatch_idx][0],
                                        c))
                {
                    identical = false;
                    break;
                }
            }

            // Extract values from OpenFOAM boundary faces
            if (identical)
            {
                patch.face_values[time_idx][c].resize(1);
                patch.face_values[time_idx][c][0] =
                        getComponent(field.boundaryField()[bfpatch_idx][0], c);

            }
            else
            {
                patch.face_values[time_idx][c].resize(nfaces_in_patch);
                for(int fidx = 0; fidx < nfaces_in_patch; fidx++)
                {
                    patch.face_values[time_idx][c][fidx] =
                        getComponent(field.boundaryField()[bfpatch_idx][fidx],
                                      c);
                }
            }
        }

        time_idx++;
    }

    runTime.setTime(start_time, start_time_index);
    field.correctBoundaryConditions();
}


void OF2SpeeditclFieldWriter::saveTimeVaryingBoundaryPatch(
    const TimeVaryingBoundaryPatch& patch, std::ostream& file)
{
    if(!patch.initialized)
        throw std::string("Can not save uninitialized time varying "
                          "boundary patch");

    file << TIME_VARYING_BOUNDARY_PATCH_HEADER << "\n" ;
    vectorWriter.saveVec(file, patch.face_positions,
                         TIME_VARYING_PATCH_FACE_POSITIONS_HEADER);

    file << "\n" << TIME_VARYING_PATCH_VALUES_HEADER;

    for(unsigned t = 0; t < patch.face_values.size(); t++)
    {
        file << "\n\nTime " << t << "\n" ;
        for(unsigned c = 0; c < patch.face_values[t].size(); c++)
        {
            std::stringstream cname;
            cname << "component_" << c;
            vectorWriter.saveVec(file, patch.face_values[t][c],
                                 cname.str().c_str());
        }
    }

    file << "\n";
    file << TIME_VARYING_TIMELINES_TIME_HEADER << "\n";
    file << "time" << "\n";
    file << VECTOR_HEADER << "\n";
    file << internal::dataTypeHeader<double>() << "\n";
    file << patch.timeline.size() << " 1" << "\n";

    for (int i = 0; i < patch.timeline.size(); i++)
        file << patch.timeline[i] << "\n";
}


template<class TField>
std::vector<int> OF2SpeeditclFieldWriter::timeVaryingPatchIndices(
    const TField& field)
{
    std::vector<int> result;

    for(int i = 0; i < field.boundaryField().size(); i++)
    {
        if(isTimeVarying(field.boundaryField().types()[i]))
        {
            result.push_back(i);
        }
    }

    return result;
}

template<class TField>
void OF2SpeeditclFieldWriter::buildGroovyBoundaryPatch(GroovyBoundaryPatch& patch,
                                                   TField& field,
                                                   const Foam::fvMesh& mesh,
                                                   Foam::Time& runTime,
                                                   int  bfpatch_idx)
{
    patch.initialized = false;

    if(!isGroovy(field.boundaryField().types()[bfpatch_idx]))
        throw std::string("Can not build groovy boundary patch from OpenFOAM "
                          "patch of type \"")
            + field.boundaryField().types()[bfpatch_idx]
            + "\", possibly wrong patch index passed";

    if((0 > bfpatch_idx) || (field.boundaryField().size() <= bfpatch_idx))
        throw std::string("Invalid patch index passed");

    buildGroovyPatchFaceIndices(patch, field, mesh, bfpatch_idx);
    buildGroovyPatchFaceValues(patch, field, mesh, runTime, bfpatch_idx);

    patch.initialized = true ;
}

template<class TField>
void OF2SpeeditclFieldWriter::buildGroovyPatchFaceIndices(
    GroovyBoundaryPatch& patch, const TField& field, const Foam::fvMesh& mesh,
    int bfpatch_idx)
{
    // Obtain indices into the aggregated, unsorted array containing
    // consecutive OpenFOAM patches.
    std::vector<int> unsorted_face_positions;
    int idx = 0;

    // Find start position (skip all previous patches)
    for(int i = 0 ; i < bfpatch_idx; i++)
    {
        idx += field.boundaryField()[i].size();
    }

    // Add all faces from current patch
    for(int i = 0; i < field.boundaryField()[bfpatch_idx].size() ; i++)
    {
        unsorted_face_positions.push_back(idx + i);
    }

    std::vector<int> bf_reorder_table = internal::bfacesReorderTable(mesh);

    std::vector<int> bf_reorder_table_2;
    bf_reorder_table_2.assign( bf_reorder_table.size(), -1);

    for(unsigned i = 0; i < bf_reorder_table.size(); i++)
    {
        bf_reorder_table_2[bf_reorder_table[i]] = i;
    }

    for(unsigned i = 0; i < unsorted_face_positions.size(); i++)
    {
        patch.face_positions.push_back(
                    bf_reorder_table_2[unsorted_face_positions[i]]);
    }
}

template<class TField>
void OF2SpeeditclFieldWriter::buildGroovyPatchFaceValues(
    GroovyBoundaryPatch & patch, TField& field, const Foam::fvMesh& mesh,
    Foam::Time& runTime, int bfpatch_idx)
{
    dimensionedScalar start_time = runTime.startTime();
    label start_time_index = runTime.startTimeIndex();

    field.correctBoundaryConditions();

    unsigned int time_idx = 0;

    cout << "<<<<<<<<< Saving time step " << time_idx << "\n";

    field.correctBoundaryConditions();

    patch.face_values.resize(n_components(field));

    int nfaces_in_patch = field.boundaryField()[ bfpatch_idx ].size() ;

    // Each component is treated separately
    for(int c = 0; c < n_components(field); c++)
    {
        // Test, if all values are the same.
        // Allows for optimisation of memory usage.
        bool identical = true;
        for(int fidx = 0; fidx < nfaces_in_patch; fidx++)
        {
            if(getComponent(field.boundaryField()[bfpatch_idx][fidx], c)
                    != getComponent(field.boundaryField()[bfpatch_idx][0], c))
            {
                identical = false;
                break;
            }
        }

        // Extract values from OpenFOAM boundary faces
        if(identical)
        {
            patch.face_values[c].resize(1);
            patch.face_values[c][0]
                    = getComponent(field.boundaryField()[bfpatch_idx][0], c);
        }
        else
        {
            patch.face_values[c].resize(nfaces_in_patch);

            for(int fidx = 0; fidx < nfaces_in_patch; fidx++)
            {
                patch.face_values[ c ][ fidx ]
                        = getComponent(
                            field.boundaryField()[bfpatch_idx][fidx], c);
            }
        }
    }

    runTime.setTime( start_time, start_time_index ) ;
    field.correctBoundaryConditions();
}

void OF2SpeeditclFieldWriter::saveGroovyBoundaryPatch(
    const GroovyBoundaryPatch& patch, std::ostream& file)
{
    if(!patch.initialized)
        throw std::string("Can not save uninitialized "
                          "time varying boundary patch");

    file << GROOVY_BOUNDARY_PATCH_HEADER << "\n";

    vectorWriter.saveVec(file, patch.face_positions,
                         GROOVY_PATCH_FACE_POSITIONS_HEADER);

    file << "\n" << GROOVY_PATCH_VALUES_HEADER;

    for(unsigned c = 0; c < patch.face_values.size(); c++)
    {
        std::stringstream cname;
        cname << "component_" << c;
        vectorWriter.saveVec(file, patch.face_values[c], cname.str().c_str());
    }
}

template<class TField>
std::vector<int> OF2SpeeditclFieldWriter::groovyPatchIndices(const TField& field)
{
    std::vector<int> result;

    for(int i = 0; i < field.boundaryField().size(); i++)
    {
        if(isGroovy(field.boundaryField().types()[i]))
        {
            result.push_back(i);
        }
    }

    return result;
}


template<class TField>
void OF2SpeeditclFieldWriter::buildTotalPressureBoundaryPatch(
    totalPressureBoundaryPatch& patch, TField& field, const Foam::fvMesh& mesh,
    Foam::Time& runTime, int bfpatch_idx)
{
    patch.initialized = false ;

    if (!isTotalPressure( field.boundaryField().types()[bfpatch_idx]))
        throw std::string("Can not build totalPressure boundary patch from "
                          "OpenFOAM patch of type \"")
            + field.boundaryField().types()[bfpatch_idx]
            + "\", possibly wrong patch index passed";

    if((0 > bfpatch_idx) || (field.boundaryField().size() <= bfpatch_idx))
        throw std::string("Invalid patch index passed");

    buildTotalPressurePatchFaceIndices(patch, field, mesh, bfpatch_idx);
    buildTotalPressurePatchFaceValues(patch, field, mesh, runTime, bfpatch_idx);

    patch.initialized = true;
}

template<class TField>
void OF2SpeeditclFieldWriter::buildTotalPressurePatchFaceIndices(
    totalPressureBoundaryPatch& patch, const TField& field,
    const Foam::fvMesh& mesh, int bfpatch_idx)
{
    // Obtain indices into the aggregated, unsorted array containing consecutive
    // OpenFOAM patches.
    std::vector<int> unsorted_face_positions;
    int idx = 0;

    // Find start position (skip all previous patches)
    for(int i = 0; i < bfpatch_idx; i++)
    {
        idx += field.boundaryField()[i].size();
    }

    // Add all faces from current patch
    for(int i = 0; i < field.boundaryField()[bfpatch_idx].size(); i++)
    {
        unsorted_face_positions.push_back(idx + i);
    }

    std::vector<int> bf_reorder_table = internal::bfacesReorderTable(mesh);

    std::vector<int> bf_reorder_table_2;
    bf_reorder_table_2.assign(bf_reorder_table.size(), -1);

    for(unsigned i = 0; i < bf_reorder_table.size(); i++)
    {
        bf_reorder_table_2[bf_reorder_table[i]] = i;
    }

    for(unsigned i = 0; i < unsorted_face_positions.size(); i++)
    {
        patch.face_positions.push_back(
                    bf_reorder_table_2[unsorted_face_positions[i]]);
    }
}

template<class TField>
void OF2SpeeditclFieldWriter::buildTotalPressurePatchFaceValues(
    totalPressureBoundaryPatch& patch, TField& field, const Foam::fvMesh& mesh,
    Foam::Time& runTime, int bfpatch_idx)
{
    dimensionedScalar start_time = runTime.startTime();
    label start_time_index = runTime.startTimeIndex();

    field.correctBoundaryConditions();

    unsigned time_idx = 0;

    cout << "<<<<<<<<< Saving time step " << time_idx << "\n";

    field.correctBoundaryConditions();

    patch.face_values.resize(n_components(field));

    int nfaces_in_patch = field.boundaryField()[bfpatch_idx].size();

    // Each component is treated separately
    for(int c = 0; c < n_components(field); c++)
    {
        // Test, if all values are the same.
        // Allows for optimisation of memory usage.
        bool identical = true;
        for(int fidx = 0; fidx < nfaces_in_patch; fidx++)
        {
            if(getComponent(field.boundaryField()[bfpatch_idx][fidx], c)
                    != getComponent(field.boundaryField()[bfpatch_idx][0], c))
            {
                identical = false;
                break;
            }
        }

        // Extract values from OpenFOAM boundary faces
        if(identical)
        {
            patch.face_values[c].resize(1);
            patch.face_values[c][0]
                    = getComponent(field.boundaryField()[bfpatch_idx][0], c);
        }
        else
        {
            patch.face_values[c].resize(nfaces_in_patch);
            for(int fidx = 0; fidx < nfaces_in_patch; fidx++)
            {
                patch.face_values[c][fidx] = getComponent(
                            field.boundaryField()[bfpatch_idx][fidx], c);
            }
        }
    }

    runTime.setTime(start_time, start_time_index);
    field.correctBoundaryConditions();
}

void OF2SpeeditclFieldWriter::saveTotalPressureBoundaryPatch(
    const totalPressureBoundaryPatch& patch, std::ostream& file)
{
    if(!patch.initialized)
        throw std::string("Can not save uninitialized time varying "
                          "boundary patch");

    file << TOTALPRESSURE_BOUNDARY_PATCH_HEADER << "\n";
    vectorWriter.saveVec(file, patch.face_positions,
                         TOTALPRESSURE_PATCH_FACE_POSITIONS_HEADER);

    file << "\n" << TOTALPRESSURE_PATCH_VALUES_HEADER;

    for(unsigned c = 0; c < patch.face_values.size(); c++)
    {
        std::stringstream cname;
        cname << "component_" << c;
        vectorWriter.saveVec(file, patch.face_values[c], cname.str().c_str());
    }
}

template<class TField>
void OF2SpeeditclFieldWriter::saveTimeVaryingBC(TField& field,
                                            const Foam::fvMesh & mesh,
                                            const std::string& dir_path,
                                            Foam::Time& runTime)
{
    std::vector<int> tv_OF_patch_indices = timeVaryingPatchIndices(field);

    if(tv_OF_patch_indices.size() == 0)
        return;

    std::string field_name (field.name());
    std::string filename = std::string(dir_path) + "/" + field_name
            + "_time_varying_boundary";

    Info << "  Saving time varying boundary conditions\n";

    std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

    if (!file.is_open())
        throw "Can not open file " + filename ;

    file << FILE_HEADER                  << "\n" ;
    file << TIME_VARYING_BOUNDARY_HEADER << "\n" ;
    file << field.name()                 << "\n" ;
    file << n_components( field )        << "\n" ; // number of components for single face (scalar/vector)
    file << internal::computeNTimeSteps(runTime) << "\n" ; // number of time steps
    file << tv_OF_patch_indices.size()   << "\n" ; // number of boundary patches

    for (unsigned int i=0 ; i<tv_OF_patch_indices.size() ; i++)
    {
        TimeVaryingBoundaryPatch tvb_patch ;

        buildTimeVaryingBoundaryPatch(tvb_patch, field, mesh, runTime,
                                      tv_OF_patch_indices[i]);
        saveTimeVaryingBoundaryPatch(tvb_patch, file);
    }
}

template<class TField>
void OF2SpeeditclFieldWriter::saveGroovyBC(TField& field, const Foam::fvMesh & mesh,
                                       const std::string& dir_path,
                                       Foam::Time& runTime)
{
    std::vector<int> g_OF_patch_indices = groovyPatchIndices(field);

    if (g_OF_patch_indices.size() == 0)
        return;

    std::string field_name(field.name()) ;
    std::string filename = std::string(dir_path) + "/" + field_name
            + "_groovy_boundary";
    std::string timelinesfilename = std::string(dir_path) + "/" + field_name
            + "_groovy_timelines";

    Info << "  Saving groovy boundary conditions\n";

    std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary );

    if(!file.is_open())
        throw "Can not open file " + filename;

    std::ofstream timelinesfile(
                timelinesfilename.c_str(), std::ios::out | std::ios::binary);

    if(!timelinesfile.is_open())
        throw "Can not open file " + timelinesfilename ;

    // Finding the name of file with timelines
    std::string tempname;
    std::string timefilename;
    std::string case_dir_name = getenv("FOAM_CASE");

    Info << "FOAM_CASE" << case_dir_name << endl;

    std::string lookfortimefilename = case_dir_name + "/0/" + field_name;

    std::ifstream lookfortimefile(lookfortimefilename.c_str());
    if (!lookfortimefile.is_open())
        throw "Can not open file " + lookfortimefilename ;



    file << FILE_HEADER                  << "\n" ;
    file << GROOVY_BOUNDARY_HEADER       << "\n" ;
    file << field.name()                 << "\n" ;
    file << n_components( field )        << "\n" ; // number of components for single face (scalar/vector)
    //file << N_TIME_STEPS               << "\n" ; // number of time steps
    file << g_OF_patch_indices.size()    << "\n" ; // number of boundary patches

    timelinesfile << FILE_HEADER                  << "\n" ;
    timelinesfile << GROOVY_TIMELINES_HEADER      << "\n" ;
    timelinesfile << field.name()                 << "\n" ;
    timelinesfile << "2"                          << "\n" ; // number of components for single face (scalar/vector)
    //timelinesfile << N_TIME_STEPS               << "\n" ; // number of time steps
    timelinesfile << g_OF_patch_indices.size()    << "\n" ; // number of boundary patches

    for (unsigned i = 0; i < g_OF_patch_indices.size(); i++)
    {
        GroovyBoundaryPatch gb_patch;

        buildGroovyBoundaryPatch(gb_patch, field, mesh, runTime,
                                 g_OF_patch_indices[i]);
        saveGroovyBoundaryPatch(gb_patch, file);

        while (lookfortimefile >> tempname)
        {
            if (tempname == "fileName")
            {
                lookfortimefile >> timefilename;
                break;
            }
        }

        if (!timefilename.empty())
        {
            timefilename.erase(std::remove(timefilename.begin(),
                                           timefilename.end(), '\"'),
                               timefilename.end());
            timefilename.erase(std::remove(timefilename.begin(),
                                           timefilename.end(), ';'),
                               timefilename.end());

            timefilename.erase(0,11);

            timefilename.insert(0, "/");
            timefilename.insert(0, case_dir_name);

            Info << "File with timelines : " << timefilename << endl;
        }

        //std::ifstream timefile (timefilename.c_str() );
        //if ( !timefile.is_open() ) throw "Can not open file " + timefilename ;

        IFstream timefile(timefilename);

        //save_groovy_timelines (timefile, timelinesfile);

        vector2DField values(timefile);

        //Info << "values : " << values[0].x() << " " << values[0].y()  << "size = " << values.size() << endl;

        //timelinesfile << "" << "\n";
        timelinesfile << GROOVY_TIMELINES_PATCH_HEADER << "\n";
        timelinesfile << "" << "\n";
        timelinesfile << GROOVY_TIMELINES_TIME_HEADER << "\n";
        timelinesfile << "time" << "\n";
        timelinesfile << VECTOR_HEADER << "\n";
        timelinesfile << internal::dataTypeHeader<double>() << "\n";
        timelinesfile << values.size() << " 1" << "\n";

        for (int i = 0; i < values.size(); i++)
            timelinesfile << values[i].x() << "\n";

        timelinesfile << "" << "\n";
        timelinesfile << GROOVY_TIMELINES_VALUES_HEADER << "\n";
        timelinesfile << "values" << "\n";
        timelinesfile << VECTOR_HEADER << "\n";
        timelinesfile << internal::dataTypeHeader<double>() << "\n";
        timelinesfile << values.size() << " 1" << "\n";

        for (int i = 0; i < values.size(); i++)
            timelinesfile << values[i].y() << "\n";
    }
}

template<class TField>
std::vector<int> OF2SpeeditclFieldWriter::totalPressurePatchIndices(
    const TField& field)
{
    std::vector<int> result;

    for (int i = 0; i<field.boundaryField().size(); i++)
    {
        if(isTotalPressure(field.boundaryField().types()[i]))
        {
            result.push_back(i);
        }
    }

    return result;
}

template<class TField>
void OF2SpeeditclFieldWriter::saveTotalPressureBC(TField& field,
                                              const Foam::fvMesh& mesh,
                                              const std::string  & dir_path,
                                              Foam::Time& runTime)
{
    std::vector<int> tP_OF_patch_indices = totalPressurePatchIndices(field);

    if(tP_OF_patch_indices.size() == 0)
        return;

    std::string field_name(field.name());
    std::string filename = std::string(dir_path) + "/" + field_name
            + "_totalPressure_boundary";

    Info << "  Saving totalPressure boundary conditions\n";

    std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

    if(!file.is_open())
        throw "Can not open file " + filename;

    file << FILE_HEADER                   << "\n" ;
    file << TOTALPRESSURE_BOUNDARY_HEADER << "\n" ;
    file << field.name()                  << "\n" ;
    file << n_components( field )         << "\n" ; // number of components for single face (scalar/vector)
    //file << N_TIME_STEPS                << "\n" ; // number of time steps
    file << tP_OF_patch_indices.size()    << "\n" ; // number of boundary patches

    for (unsigned i = 0; i < tP_OF_patch_indices.size(); i++)
    {
        totalPressureBoundaryPatch tPb_patch;

        buildTotalPressureBoundaryPatch(tPb_patch, field, mesh, runTime,
                                        tP_OF_patch_indices[i]);

        saveTotalPressureBoundaryPatch(tPb_patch, file);
    }
}

} // namespace io
} // namespace speeditcl
