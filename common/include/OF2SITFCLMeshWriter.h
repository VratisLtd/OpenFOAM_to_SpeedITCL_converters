#ifndef OF_2_SITFCL_MESH_WRITER_H
#define OF_2_SITFCL_MESH_WRITER_H

#include <fstream>
#include <vector>

#include "OF2SITFCLVectorWriter.h"

namespace speeditcl
{
namespace io
{

class OF2SpeeditclMeshWriter
{
public:

    explicit OF2SpeeditclMeshWriter(OF2SpeeditclVectorWriter& vectorWriter)
        : vectorWriter(vectorWriter) {}

    void saveMesh(fvMesh& mesh, const std::string& dir_path);

private:

    void saveMeshDescription(Foam::fvMesh& mesh, std::ostream& file);
    std::vector<int> boundaryOwners(Foam::fvMesh const& mesh);

    OF2SpeeditclVectorWriter& vectorWriter;
};

} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_MESH_WRITER_H
