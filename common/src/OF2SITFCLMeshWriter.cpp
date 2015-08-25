#include "OF2SITFCLMeshWriter.h"
#include "OF2SITFCLFieldInternals.h"
#include "OF2SITFCLBoundary.h"
#include "OF2SITFCLErrorHandling.h"

namespace speeditcl
{
namespace io
{

using namespace internal;

void OF2SpeeditclMeshWriter::saveMesh(fvMesh& mesh, const std::string& dir_path)
{
    std::string filename = dir_path + "/" + FILE_NAME_MESH ;

    try
    {
        std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);

        if (!file.is_open()) throw "Can not open file " + filename;

        file << FILE_HEADER << "\n";
        file << MESH_HEADER << "\n";

        saveMeshDescription(mesh, file) ;

        vectorWriter.saveVec(file, internalValues(mesh.Sf()),
                         IFACE_AREA_VECTORS_HEADER);

        vectorWriter.saveVec(file, mesh.weights().cdata(),
                             mesh.weights().size(), IFACE_WEIGHTS_HEADER);

        vectorWriter.saveVec(file, mesh.deltaCoeffs().cdata(),
                         mesh.deltaCoeffs().size(), IFACE_DELTA_COEFFS_HEADER);

#if (2 == OF_VERSION && 1 > OF_VERSION_MINOR ) || 1 == OF_VERSION
        if(!mesh.orthogonal())
        {
            vectorWriter.saveVec(file,
                                 internalValues(mesh.correctionVectors()),
                                 IFACE_CORRECTION_VEC_HEADER ) ;
        }
#endif
#if (2 == OF_VERSION && 1 <= OF_VERSION_MINOR )
        if(!mesh.checkFaceOrthogonality())
        {
            vectorWriter.saveVec(file,
                                 internalValues(mesh.nonOrthCorrectionVectors()),
                                 IFACE_CORRECTION_VEC_HEADER ) ;
        }
#endif
        vectorWriter.saveVec(file, mesh.owner().cdata(), mesh.owner().size(),
                             IFACE_OWNERS_HEADER);

        vectorWriter.saveVec(file, mesh.neighbour().cdata(),
                             mesh.neighbour().size(),
                             IFACE_NEIGHBOURS_HEADER);

        vectorWriter.saveVec(file, boundaryValues(mesh.Sf(), true).vals,
                             BFACE_AREA_VECTORS_HEADER);

        vectorWriter.saveVec(file, boundaryValues(mesh.weights(), true).vals,
                             BFACE_WEIGHTS_HEADER);

        vectorWriter.saveVec(file,
                             boundaryValues(mesh.deltaCoeffs(), true).vals,
                             BFACE_DELTA_COEFFS_HEADER);

#if (2 == OF_VERSION && 1 > OF_VERSION_MINOR ) || 1 == OF_VERSION
        if(!mesh.orthogonal())
        {

            auto& boundaryField = mesh.correctionVectors().boundaryField();

            for (int i = 0; i < boundaryField.size(); i++)
                for (int j=0 ; j < boundaryField[i].size() ; j++)
                {
                    if (0 != boundaryField[i][j].x()
                            || 0 != boundaryField[i][j].y()
                            || 0 != boundaryField[i][j].z())
                        throw std::string("boundary correction vectors "
                                          "different than 0 are not supported"
                                          " - do you try to use \"coupled\" "
                                          "boundary faces ?") ;
                }

            vectorWriter.saveVec(file, boundaryValues(mesh.correctionVectors(),
                                                   true).vals,
                             BFACE_CORRECTION_VEC_HEADER);
        }
#endif
#if (2 == OF_VERSION && 1 <= OF_VERSION_MINOR )
        if(!mesh.checkFaceOrthogonality())
        {

            auto& boundaryField = mesh.nonOrthCorrectionVectors().boundaryField();

            for (int i = 0; i < boundaryField.size(); i++)
                for (int j=0 ; j < boundaryField[i].size() ; j++)
                {
                    if (0 != boundaryField[i][j].x()
                            || 0 != boundaryField[i][j].y()
                            || 0 != boundaryField[i][j].z())
                        throw std::string("boundary correction vectors "
                                          "different than 0 are not supported"
                                          " - do you try to use \"coupled\" "
                                          "boundary faces ?") ;
                }

            vectorWriter.saveVec(file, boundaryValues(mesh.nonOrthCorrectionVectors(),
                                                   true).vals,
                             BFACE_CORRECTION_VEC_HEADER);
        }
#endif

        vectorWriter.saveVec(file, boundaryOwners(mesh), BFACE_OWNERS_HEADER);

        vectorWriter.saveVec(file, bfacesReorderTable(mesh),
                         BFACE_REORDER_TABLE_HEADER);

        vectorWriter.saveVec(file, mesh.V().cdata(), mesh.V().size(),
                         CELL_VOLUMES_HEADER);
    }
    catch (std::string & err_str)
    {
        handleError("save_mesh ERROR : " + err_str, dir_path);
    }
}


void OF2SpeeditclMeshWriter::saveMeshDescription(Foam::fvMesh& mesh,
                                             std::ostream& file)
{
    file << mesh.nCells() << " ";  // number of cells
    file << mesh.nInternalFaces() << " "; // number of internal faces
    std::vector<int> bowners = boundaryOwners(mesh);
    file << bowners.size();

    file << "\n";

#if (2 == OF_VERSION && 1 > OF_VERSION_MINOR ) || 1 == OF_VERSION
        if ( mesh.orthogonal() )
        {
            file << "orthogonal\n" ;
        } else {
            file << "nonorthogonal\n" ;
        } ;
#endif
#if (2 == OF_VERSION && 1 <= OF_VERSION_MINOR )
        if ( mesh.checkFaceOrthogonality() )
        {
            file << "orthogonal\n" ;
        } else {
            file << "nonorthogonal\n" ;
        } ;
#endif

    // Some kind of magic...
    Vector<scalar>::labelType validComponents(
                pow(mesh.solutionD(),
                    pTraits<powProduct<Vector<label>,
                    Vector<scalar>::rank>::type>::zero));

    file << "X  " << (1 == validComponents[0]) << "\n";
    file << "Y  " << (1 == validComponents[1]) << "\n";
    file << "Z  " << (1 == validComponents[2]);

    testIO(file);
}

std::vector<int> OF2SpeeditclMeshWriter::boundaryOwners(Foam::fvMesh const & mesh)
{
    std::vector<int> v;
    auto& boundary = mesh.boundary();

    for(int z = 0; z < boundary.size(); z++)
    {
        for (int i = 0; i < boundary[z].faceCells().size(); i++)
        {
            v.push_back(boundary[z].faceCells()[i]);
        }
    }

    std::vector<int> vbp = bfacesReorderTable(mesh);
    return vectorReorder(v, vbp);
}

} // namespace io
} // namespace speeditcl
