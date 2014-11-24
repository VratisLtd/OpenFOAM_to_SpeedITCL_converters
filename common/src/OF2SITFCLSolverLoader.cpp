#include "OF2SITFCLSolverLoader.h"
#include "OF2SITFCLMeshWriter.h"
#include "OF2SITFCLFieldWriter.h"
#include "OF2SITFCLSettings.h"
#include "wallDist.H"
#include "nearWallDist.H"

#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"


#if OF_VERSION == 2
        #include "simpleControl.H"
#endif

void loadOF_Stationary(Foam::Time& runTime,
                   Foam::fvMesh& mesh,
                   volVectorField& U,
                   volScalarField& p,
                   surfaceScalarField& phi,
                   const std::string& outputPath,
                   speeditcl::io::OF2SpeeditclVectorWriter& vectorWriter,
                   pt::ptree& settings)
{
#if 2 == OF_VERSION
    simpleControl simple(mesh);
    #endif
#if 1 == OF_VERSION
    #include "readSIMPLEControls.H"
#endif

    pt::ptree simpleSettings;

#if OF_VERSION == 1
    simpleSettings.put(speeditcl::io::SETTING_N_NONORTH_CORR_ITER, nNonOrthCorr);
#elif OF_VERSION == 2
    simpleSettings.put(speeditcl::io::SETTING_N_NONORTH_CORR_ITER, simple.nNonOrthCorr());
#endif

    //read resudualControl;
    dictionary residualControls  = mesh.solutionDict().
            subDict("SIMPLE").subOrEmptyDict("residualControl");

    scalar pRequiredResidual = 0.0;
    scalar URequiredResidual = 0.0;

    residualControls.readIfPresent<scalar>("p", pRequiredResidual);
    residualControls.readIfPresent<scalar>("U", URequiredResidual);

    simpleSettings.put(speeditcl::io::SETTING_P_REQUIRED_RESIDUAL, pRequiredResidual);

    simpleSettings.put(speeditcl::io::SETTING_U_REQUIRED_RESIDUAL, URequiredResidual);

    settings.push_front(pt::ptree::value_type("SIMPLE", simpleSettings));

    //TODO:: add residuals for k, omega in case of turbulence


    speeditcl::io::OF2SpeeditclMeshWriter meshWriter(vectorWriter);
    speeditcl::io::OF2SpeeditclFieldWriter fieldWriter(vectorWriter);

    Info << "Saving mesh in SpeedITCL format" << endl;
    meshWriter.saveMesh(mesh, outputPath);
    Info << "Saving p in SpeedITCL format" << endl;
    fieldWriter.saveField(p, mesh, outputPath);
    Info << "Saving U in SpeedITCL format" << endl;
    fieldWriter.saveField(U, mesh, outputPath);
    Info << "Saving phi in SpeedITCL format" << endl;
    fieldWriter.saveField(phi, mesh, outputPath);
}

void loadOF_Transient(Foam::Time& runTime,
                Foam::fvMesh& mesh,
                volVectorField& U,
                volScalarField& p,
                surfaceScalarField& phi,
                const std::string& outputPath,
                speeditcl::io::OF2SpeeditclVectorWriter& vectorWriter,
                pt::ptree& settings)
{
    const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");

    const int nCorrectors =
        pisoDict.lookupOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorrectors =
        pisoDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    pt::ptree icoSettings;

    icoSettings.put(speeditcl::io::SETTING_N_PISO_CORR_ITER, nCorrectors);

    icoSettings.put(speeditcl::io::SETTING_N_NONORTH_CORR_ITER, nNonOrthCorrectors);

    settings.push_front(pt::ptree::value_type("PISO", icoSettings));

    speeditcl::io::OF2SpeeditclMeshWriter meshWriter(vectorWriter);
    speeditcl::io::OF2SpeeditclFieldWriter fieldWriter(vectorWriter);

    Info << "Saving mesh in SpeedITCL format" << endl;
    meshWriter.saveMesh(mesh, outputPath);
    Info << "Saving p in SpeedITCL format" << endl;
    fieldWriter.saveField(p, mesh, outputPath, &runTime);
    Info << "Saving U in SpeedITCL format" << endl;
    fieldWriter.saveField(U, mesh, outputPath, &runTime);
    Info << "Saving phi in SpeedITCL format" << endl;
    fieldWriter.saveField(phi, mesh, outputPath, &runTime);
}


void loadOF_Turbulence(Foam::Time& runTime,
                       Foam::fvMesh& mesh,
                       const std::string& outputPath,
                       speeditcl::io::OF2SpeeditclVectorWriter& vectorWriter,
                       pt::ptree& settings)
{

    Info<< "Reading field k" << endl;

    volScalarField k(IOobject("k", runTime.timeName(), mesh,
                              IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);

    Info<< "Reading field nut" << endl;
    volScalarField nut(IOobject("nut", runTime.timeName(), mesh,
                                IOobject::MUST_READ, IOobject::AUTO_WRITE),
                       mesh);

    Info<< "Reading field omega\n" << endl;
    volScalarField omega(IOobject("omega", runTime.timeName(), mesh,
                                  IOobject::MUST_READ, IOobject::AUTO_WRITE),
                         mesh);

    // obliczanie odleglosci komorek od scianek zew. ?
    Info<< "Calculating y (wallDistance)  field\n" << endl;
    wallDist wallDistance(mesh, true);

    nearWallDist nearWallDistance(mesh); //for boundary

    //to distinct field name from wallDistance we name it yNW
    volScalarField nearWallDistanceField ("yNW", wallDistance);
    nearWallDistanceField.boundaryField() = nearWallDistance.y();

    speeditcl::io::OF2SpeeditclFieldWriter fieldWriter(vectorWriter);

    Info << "Saving k in SpeedITCL format\n";
    fieldWriter.saveField(k, mesh, outputPath, &runTime);
    Info << "Saving omega in SpeedITCL format\n";
    fieldWriter.saveField(omega, mesh, outputPath, &runTime);
    Info << "Saving wallDistance(y) in SpeedITCL format\n";
    fieldWriter.saveField(wallDistance, mesh, outputPath, &runTime);
    Info << "Saving nearWallDistance(yNW) in SpeedITCL format\n" ;
    fieldWriter.saveField(nearWallDistanceField, mesh, outputPath, &runTime );
    Info << "Saving nut in SpeedITCL format\n";
    fieldWriter.saveField(nut, mesh, outputPath, &runTime);

    //save turbulence settings
    speeditcl::io::saveRelaxationParameter<volScalarField>(k, settings);
    speeditcl::io::saveRelaxationParameter<volScalarField>(omega, settings);


}
