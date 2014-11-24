/*---------------------------------------------------------------------------*\
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\*---------------------------------------------------------------------------*/

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#ifndef OF_VERSION
	#error "Undefined OpenFOAM version"
#endif

#if OF_VERSION != 1 && OF_VERSION != 2
	#error "Unsupported OpenFOAM version"
#endif

#include <string>
#include <iostream>
#include <array>

#include "OF2SITFCLSettings.h"
#include "OF2SITFCLUtils.h"
#include "OF2SITFCLMeshWriter.h"
#include "OF2SITFCLFieldWriter.h"
#include "OF2SITFCLErrorHandling.h"
#include "OF2SITFCLSolverLoader.h"


//#if OF_VERSION == 2
//	#include "simpleControl.H"
//#endif


namespace po = boost::program_options;
namespace pt = boost::property_tree;


void getCellReferenceValues(const word app_name,
                            const Foam::fvMesh& mesh,
                            const volScalarField& p,
                            label& pRefCell,
                            scalar& pRefValue);



//returns true if RASproperties says to use turbulence model
//currently only kOmegaSST
bool Turbulence(Foam::Time& runTime, Foam::fvMesh& mesh);

int main(int argc, char *argv[])
{
    std::cout << "OPEN FOAM TO SPEEDIT FLOWCL CONVERTER" << std::endl;

    std::string casePath = "./";
    std::string outputPath = speeditcl::io::DATA_DIR_NAME;
    bool binaryMode = false;
    int precision = speeditcl::io::OF2SpeeditclVectorWriter::defaultPrecision;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help",
         "produce help message")
        ("case", po::value(&casePath)->default_value(casePath),
         "path to case")
        ("out", po::value(&outputPath)->default_value(outputPath),
         "output path")
        ("binary", po::value(&binaryMode)->default_value(binaryMode),
         "use binary mode")
        ("precision", po::value(&precision)->default_value(precision),
         "number of digits in non-binary mode")
    ;

    /** parsing arguments */
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    /** printing help message */
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }

    /** creating output directory if doesn't exist */
    speeditcl::io::createDirectory(outputPath.c_str());

    /** creating writers */
    speeditcl::io::OF2SpeeditclVectorWriter vectorWriter(binaryMode, precision);

    /** preparing parameters for OpenFoam argument processor */
    const char* arguments[3] = {argv[0], "-case", casePath.c_str()};
    int argumentsSize = sizeof(arguments) / sizeof(char*);

    try
    {
        /** disabling stadnard OpenFoam's banner */
        Foam::argList::noBanner();

        /** passing prepraed arguments */
        char** m_args = const_cast<char**>(arguments);
        Foam::argList args(argumentsSize, m_args);

        /** Parallel run is not supported */
        if(args.optionFound("parallel"))
        {
            throw std::string("Can not convert parallel case, please reorganize"
                              " your case to sigle-domain");
        }

        /** checks if -case directory can be opened */
        if(!args.checkRootCase())
        {
            throw std::string("Wrong parameters passed,"
                              "OpenFOAM check failed");
        }

        Foam::Info << "Create time" << Foam::endl;
        Foam::Time runTime(Foam::Time::controlDictName, args);

        ///** checks consistency with 'application' in 'controlDict' file */
        //check_application_name (runTime, "simpleFoam") ;

        word app_name = runTime.controlDict().lookup("application");
        Foam::Info << "Application name: " << app_name << endl;

        std::string solver_name = speeditcl::appToSolver(app_name);

        Foam::Info << "Solver name: " << solver_name << endl;

        /** from createMesh.H */
        Foam::Info
            << "Create mesh for time = "
            << runTime.timeName() << Foam::endl;

        Foam::fvMesh mesh(Foam::IOobject(Foam::fvMesh::defaultRegion,
                                         runTime.timeName(), runTime,
                                         Foam::IOobject::MUST_READ));
        /** createMesh.H end */


        Info << "Creating settings file \n" ;
        std::string filename = outputPath + "/" + speeditcl::io::FILE_NAME_SETTINGS;
        std::ofstream settingsFile(filename.c_str());

        if (!settingsFile.is_open())
            throw std::string("Can not open file " + filename);


        Info<< "Reading transportProperties" << endl;
        IOdictionary transportProperties(
            IOobject("transportProperties", runTime.constant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE));

        dimensionedScalar nu(transportProperties.lookup("nu"));


        volScalarField p(IOobject(
                             "p",
                             runTime.timeName(),
                             mesh,
                             IOobject::MUST_READ,
                             IOobject::AUTO_WRITE
                             ),
                         mesh);

        volVectorField U(IOobject(
                             "U",
                             runTime.timeName(),
                             mesh,
                             IOobject::MUST_READ,
                             IOobject::AUTO_WRITE
                             ),
                         mesh);

        surfaceScalarField phi(IOobject(
                               "phi",
                               runTime.timeName(),
                               mesh,
                               IOobject::READ_IF_PRESENT,
                               IOobject::AUTO_WRITE
                               ),
                linearInterpolate(U) & mesh.Sf());

        //READ RAS Properties if they are present, there we can deterimine model used in solver
        // if not present the laminar model is used;

        IOdictionary RASProperties(
            IOobject("RASProperties", runTime.constant(), mesh,
            IOobject::READ_IF_PRESENT, IOobject::NO_WRITE));

        word RASModel(RASProperties.lookupOrDefault<word>("RASModel", speeditcl::io::TURBULENCE_LAMINAR));
        bool RASturbulence(RASProperties.lookupOrDefault<bool>("turbulence", false));

        Info << "RASModel: " << RASModel << " turbulence: " << RASturbulence << nl;

        // these two conditions needs to be fulfilled to use turbulence model in solver
        bool turbulence = (RASModel == speeditcl::io::TURBULENCE_K_OMEGA_SST) &&
                (RASturbulence == true);

        label pRefCell = 0;
        scalar pRefValue = 0.0;

        getCellReferenceValues(app_name, mesh, p, pRefCell, pRefValue);


        pt::ptree root;
        speeditcl::io::saveSettingsCommon(runTime, speeditcl::appToDictName(app_name),
                                      root, pRefCell, pRefValue, nu, turbulence);

        //each solver saves his required fields;
        if(solver_name == speeditcl::io::SOLVER_SIMPLE)
        {

            loadOF_Stationary(runTime, mesh, U, p, phi, outputPath, vectorWriter,
                          root);

            speeditcl::io::saveRelaxationParameter<volVectorField>(U, root);
            speeditcl::io::saveRelaxationParameter<volScalarField>(p, root);

            speeditcl::io::saveIterSolversSettings(root, mesh, "U");
            speeditcl::io::saveIterSolversSettings(root, mesh, "p");

            if(turbulence)
            {
                speeditcl::io::saveIterSolversSettings(root, mesh, "k");
                speeditcl::io::saveIterSolversSettings(root, mesh, "omega");
                speeditcl::io::saveTurbulenceParameters(RASProperties, root);

                loadOF_Turbulence(runTime, mesh, outputPath, vectorWriter, root);
            }

        }
        else if(solver_name == speeditcl::io::SOLVER_PISO)
        {
            loadOF_Transient(runTime, mesh, U, p, phi, outputPath, vectorWriter,
                       root);

            speeditcl::io::saveRelaxationParameter<volVectorField>(U, root);
            speeditcl::io::saveRelaxationParameter<volScalarField>(p, root);

            speeditcl::io::saveIterSolversSettings(root, mesh, "U");
            speeditcl::io::saveIterSolversSettings(root, mesh, "p");

            if(turbulence)
            {
                speeditcl::io::saveIterSolversSettings(root, mesh, "k");
                speeditcl::io::saveIterSolversSettings(root, mesh, "omega");
                speeditcl::io::saveTurbulenceParameters(RASProperties, root);

                loadOF_Turbulence(runTime, mesh, outputPath, vectorWriter, root);
            }
        }
        else if(solver_name == speeditcl::io::SOLVER_ICO)
        {
            loadOF_Transient(runTime, mesh, U, p, phi, outputPath, vectorWriter,
                       root);

            speeditcl::io::saveRelaxationParameter<volVectorField>(U, root);
            speeditcl::io::saveRelaxationParameter<volScalarField>(p, root);

            speeditcl::io::saveIterSolversSettings(root, mesh, "U");
            speeditcl::io::saveIterSolversSettings(root, mesh, "p");

        }

    pt::write_ini(settingsFile, root);

    }
    catch(std::string& err_str)
    {
        speeditcl::io::handleError(err_str, outputPath);
    }


    Foam::Info << endl << "DONE!" << endl;

    return 0;
}

void getCellReferenceValues(const word app_name,
                            const Foam::fvMesh& mesh,
                            const volScalarField& p,
                            label& pRefCell,
                            scalar& pRefValue)
{
    /** GeometricField.C
        Search all boundary conditions, if any are
        fixed-value or mixed (Robin) do not set reference level for solution.
    */
    if(p.needReference())
    {
        const dictionary& Dict =
                mesh.solutionDict().subDict(speeditcl::appToDictName(app_name));

        if (!Dict.readIfPresent(p.name() +  "RefCell", pRefCell))
        {
            throw std::string(p.name() +  "RefCell is required for currend boundary conditions configuration");
        }
        if (!Dict.readIfPresent(p.name() +  "RefValue", pRefValue))
        {
            throw std::string(p.name() +  "RefValue is required for currend boundary conditions configuration");
        }
    }
}


