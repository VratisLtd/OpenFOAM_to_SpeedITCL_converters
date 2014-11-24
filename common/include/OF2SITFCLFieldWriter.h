#ifndef OF_2_SITFCL_FIELD_WRITER_H
#define OF_2_SITFCL_FIELD_WRITER_H

/* OpenFoam headers */
#include "fvCFD.H"
#include "vector2DField.H"
#include "IFstream.H"
#include "wallDist.H"

#include "OF2SITFCLVectorWriter.h"

namespace speeditcl
{
namespace io
{

class OF2SpeeditclFieldWriter
{
public:

    explicit OF2SpeeditclFieldWriter(OF2SpeeditclVectorWriter& vectorWriter)
        : vectorWriter(vectorWriter) {}

    void saveField(volScalarField& field, const Foam::fvMesh& mesh,
                    const std::string& dir_path, Foam::Time* runTime = nullptr);

    void saveField(volVectorField& field, const Foam::fvMesh& mesh,
                   const std::string& dir_path, Foam::Time* runTime = nullptr);

    void saveField(surfaceScalarField& field, const Foam::fvMesh& mesh,
                   const std::string& dir_path, Foam::Time* runTime = nullptr);

    void saveField(wallDist& field, const Foam::fvMesh& mesh,
                   const std::string& dir_path, Foam::Time* runTime = nullptr);

private:

    struct TimeVaryingBoundaryPatch
    {
        bool initialized = false;

        // Indices of boundary faces in GPU memory, which should be changed
        std::vector<int> face_positions;

        std::vector<        // consecutive time steps
        std::vector<        // boundary face values for single time step - each component (x, y, z)
        std::vector<double> // is a vector on current level
        > > face_values;


        std::vector<double> timeline;

        // Example: for 10 time steps and 12 faces values, which are 3D vectors
        // face_values is an array of dimensions [10][3][12]
        // This way allows to save memory and disc space, when all face values
        // for given time step are the same.
    };

    struct GroovyBoundaryPatch
    {
        bool initialized = false;


        std::vector<int> face_positions;

        std::vector<            // boundary face values for single time step - each component (x, y, z)
        std::vector<double>   // is a vector on current level
        > face_values;

        // This way allows to save memory and disc space, when all face values
        // for given time step are the same.
    };

    struct totalPressureBoundaryPatch
    {
        bool initialized = false;
        std::vector<int> face_positions ; // Indices of boundary faces in GPU memory, which should be changed

         std::vector< // boundary face values for single time step - each component (x, y, z)
         std::vector<double>   // is a vector on current level
         > face_values;

         // This way allows to save memory and disc space, when all face values
         //for given time step are the same.
    };

    inline Foam::scalar getComponent(const Foam::vector& v, int cidx)
    {
        return v.component(cidx);
    }

    inline Foam::scalar getComponent(double d, int cidx)
    {
        if (0 != cidx)
            throw std::string("Can not extract component greater "
                              "than 0 from scalar");
        return d;
    }

    inline int n_components(const volScalarField&)
    {
        return 1;
    }

    inline int n_components(const volVectorField&)
    {
        return 3;
    }

    inline int n_components(const surfaceScalarField&)
    {
        return 1;
    }

    template<class TField>
    void saveField(TField& field, const Foam::fvMesh& mesh,
                   const std::string& dir_path,
                   Foam::Time* runTime, bool save_data = true);


    bool isTimeVarying(const std::string& type);
    bool isGroovy(const std::string& type);
    bool isTotalPressure(const std::string& type);

    template<class TField>
    void buildTimeVaryingBoundaryPatch(TimeVaryingBoundaryPatch& patch,
                                       TField& field,
                                       const Foam::fvMesh& mesh,
                                       Foam::Time& runTime,
                                       int bfpatch_idx);

    template<class TField>
    void buildTimeVaryingPatchFaceIndices(TimeVaryingBoundaryPatch& patch,
                                          const TField& field,
                                          const Foam::fvMesh& mesh,
                                          int bfpatch_idx);

    template<class TField>
    void buildTimeVaryingPatchFaceValues(TimeVaryingBoundaryPatch& patch,
                                         TField& field,
                                         const Foam::fvMesh& mesh,
                                         Foam::Time& runTime,
                                         int bfpatch_idx);

    void saveTimeVaryingBoundaryPatch(const TimeVaryingBoundaryPatch & patch,
                                      std::ostream& file);

    template<class TField>
    std::vector<int> timeVaryingPatchIndices(const TField & field);

    template<class TField>
    void buildGroovyBoundaryPatch(GroovyBoundaryPatch & patch, TField& field,
                                  const Foam::fvMesh& mesh, Foam::Time& runTime,
                                  int  bfpatch_idx);

    template<class TField>
    void buildGroovyPatchFaceIndices(GroovyBoundaryPatch & patch,
                                     const TField& field,
                                     const Foam::fvMesh& mesh, int bfpatch_idx);

    template<class TField>
    void buildGroovyPatchFaceValues(GroovyBoundaryPatch & patch,
                                    TField& field, const Foam::fvMesh& mesh,
                                    Foam::Time& runTime, int bfpatch_idx);


    void saveGroovyBoundaryPatch(const GroovyBoundaryPatch& patch,
                                 std::ostream& file);

    template<class TField>
    std::vector<int> groovyPatchIndices(const TField& field);

    template<class TField>
    std::vector<int> totalPressurePatchIndices(const TField& field);

    template<class TField>
    void buildTotalPressureBoundaryPatch(totalPressureBoundaryPatch& patch,
                                         TField& field,
                                         const Foam::fvMesh& mesh,
                                         Foam::Time& runTime,
                                         int bfpatch_idx);

    template<class TField>
    void buildTotalPressurePatchFaceIndices(totalPressureBoundaryPatch& patch,
                                            const TField& field,
                                            const Foam::fvMesh& mesh,
                                            int bfpatch_idx);

    template<class TField>
    void buildTotalPressurePatchFaceValues(totalPressureBoundaryPatch& patch,
                                           TField& field,
                                           const Foam::fvMesh& mesh,
                                           Foam::Time& runTime,
                                           int bfpatch_idx);

    void saveTotalPressureBoundaryPatch(const totalPressureBoundaryPatch& patch,
                                        std::ostream& file);

    template<class TField>
    void saveTimeVaryingBC(TField& field, const Foam::fvMesh & mesh,
                           const std::string  & dir_path, Foam::Time& runTime);

    template<class TField>
    void saveGroovyBC(TField& field, const Foam::fvMesh & mesh,
                      const std::string& dir_path, Foam::Time& runTime);

    template<class TField>
    void saveTotalPressureBC(TField& field, const Foam::fvMesh& mesh,
                             const std::string  & dir_path,
                             Foam::Time& runTime);

    // Empty function used for phi field, where time varying boundary conditions
    // should not occur.
    inline void saveTimeVaryingBC(surfaceScalarField &, const Foam::fvMesh&,
                                  const std::string&, Foam::Time&){}

    inline void saveGroovyBC(surfaceScalarField&, const Foam::fvMesh&,
                             const std::string&, Foam::Time&){}

    inline void saveTotalPressureBC(surfaceScalarField&, const Foam::fvMesh&,
                                    const std::string&, Foam::Time&){}

    OF2SpeeditclVectorWriter& vectorWriter;
};

} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_FIELD_WRITER_H
