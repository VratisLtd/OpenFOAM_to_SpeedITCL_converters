#ifndef OF_2_SITFCL_SOLVER_LOADER_H
#define OF_2_SITFCL_SOLVER_LOADER_H

#include "fvCFD.H"
#include "OF2SITFCLConstants.h"
#include "OF2SITFCLVectorWriter.h"

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;


void loadOF_Stationary(Foam::Time& runTime,
                   Foam::fvMesh& mesh,
                   volVectorField& U,
                   volScalarField& p,
                   surfaceScalarField& phi,
                   const std::string& outputPath,
                   speeditcl::io::OF2SpeeditclVectorWriter& vectorWriter,
                   pt::ptree& settings);

void loadOF_Transient(Foam::Time& runTime,
                Foam::fvMesh& mesh,
                volVectorField& U,
                volScalarField& p,
                surfaceScalarField& phi,
                const std::string& outputPath,
                speeditcl::io::OF2SpeeditclVectorWriter& vectorWriter,
                pt::ptree& settings);

void loadOF_Turbulence(Foam::Time& runTime,
                       Foam::fvMesh& mesh,
                       const std::string& outputPath,
                       speeditcl::io::OF2SpeeditclVectorWriter& vectorWriter,
                       pt::ptree& settings);


#endif //OF_2_SITFCL_SOLVER_LOADER_H
