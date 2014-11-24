#include "OF2SITFCLSettings.h"

#include "OF2SITFCLUtils.h"


namespace speeditcl
{
namespace io
{

namespace internal
{

bool testIterativeSolver(const Foam::dictionary& solver)
{
        word s(solver.lookup("solver"));

        if("PCG" != s && "PBiCG" != s )
            throw std::string("Unsupported solver ") + s
                + ", the only supported solvers are: PCG PBiCG";

        word prec( solver.lookup("preconditioner"));

        if ("diagonal" != prec && "none" != prec && "AMG" != prec)
            throw std::string("Unsupported preconditioner ") + prec
                + " for solver " + s +  ", the only supported preconditioners"
                " are: diagonal AMG";

        return true;
}

struct IterSolverSettings
{
    IterSolverSettings(const dictionary & dict)
    {
        testIterativeSolver(dict);

        max_iterations = dict.lookupOrDefault<label> ("maxIter", 1000 ) ;
        required_tolerance = dict.lookupOrDefault<scalar>("tolerance", 1e-6);
        relative_tolerance = dict.lookupOrDefault<scalar>("relTol", 0.0);


        word prec(dict.lookup("preconditioner"));
        preconditioner = std::string(prec);
    }

    int max_iterations;
    double required_tolerance;
    double relative_tolerance;
    std::string preconditioner;
    unsigned amg_rebuild_couter = 550;
};

} // namespace internal



void saveIterSolversSettings(pt::ptree& settings, Foam::fvMesh& mesh, const std::string& eqn)
{
    const dictionary& fvSolution = mesh.solutionDict();

    dictionary solvers(fvSolution.subDict("solvers"));

    {
        pt::ptree SolverSettings;

        internal::IterSolverSettings iterSolverSettings(solvers.subDict(eqn)) ;

        SolverSettings.put(SETTING_EQN_MAX_ITER,
                            iterSolverSettings.max_iterations);
        SolverSettings.put(SETTING_EQN_REQ_TOLERANCE,
                            iterSolverSettings.required_tolerance);
        SolverSettings.put(SETTING_EQN_REL_TOLERANCE,
                            iterSolverSettings.relative_tolerance);
        SolverSettings.put(SETTING_EQN_PRECOND,
                            iterSolverSettings.preconditioner);

        if(eqn == "p")
            SolverSettings.put(SETTING_AMG_PRECOND_REBUILD_COUNTER,
                                iterSolverSettings.amg_rebuild_couter);

        settings.push_back(pt::ptree::value_type(eqn, SolverSettings));
    }


}

void saveSettingsCommon(Foam::Time &runTime, const std::string &solver_name,
                        pt::ptree& settings, label pRefCell, scalar pRefValue,
                        dimensionedScalar& nu, const bool turbulence)
{
    pt::ptree commons;

    commons.put(SETTING_SOLVER_NAME, solver_name);
//    commons.put(SETTING_N_TIME_STEPS,
//                 static_cast<int>(internal::computeNTimeSteps(runTime)));
    commons.put(SETTING_P_FIXED_CELL_IDX, static_cast<int>(pRefCell));
    commons.put(SETTING_P_FIXED_CELL_VALUE, pRefValue);
    commons.put(SETTING_NU, nu.value());
    commons.put(SETTING_CASE_DIR, runTime.path());
    commons.put(SETTING_TURBULENCE_MODE, turbulence);

    settings.push_back(pt::ptree::value_type("commons", commons));

    //TODO: here put it to external function
    // read time table in case of timedependent sim
    pt::ptree runControl;
    const dictionary& controlDict = runTime.controlDict();

    scalar startTime = runTime.startTime().value();
    scalar endTime = 0.0;

    if(startTime != scalar(0))
    {
        Foam::Warning << "Starting from time different than 0 is not supported.\n"
                       <<"Starting from 0" << Foam::nl;
        startTime = scalar(0);
    }
    endTime = runTime.endTime().value();

    scalar deltaT;
    controlDict.lookup("deltaT") >> deltaT;

    runControl.put(speeditcl::io::SETTING_START_TIME, startTime);
    runControl.put(speeditcl::io::SETTING_END_TIME, endTime);
    runControl.put(speeditcl::io::SETTING_DELTA_T, deltaT);

    //correct are timeStep and adjustableRunTime;
    word writeControl(controlDict.lookup("writeControl"));

    if(writeControl != "timeStep" && writeControl != "adjustableRunTime")
    {
        std::stringstream err;
        err << "writeControl value other than"
            << "[timeStep] and [adjustableRunTime]"
            << " are not supported.";
        throw err.str();

    }

    runControl.put(speeditcl::io::SETTING_WRITE_CONTROL, writeControl);

    scalar writeInterval;

    controlDict.lookup("writeInterval") >> writeInterval;

    runControl.put(speeditcl::io::SETTING_WRITEINTERVAL, writeInterval);

    bool adjustableRunTime(controlDict.lookupOrDefault("adjustableRunTime", false));

    if(adjustableRunTime)
    {
        runControl.put(speeditcl::io::SETTING_ADJUSTABLE_RUN_TIME, adjustableRunTime);

        scalar maxCo = controlDict.lookupOrDefault("maxCo", 1.0);
        runControl.put(speeditcl::io::SETTING_MAX_CO, maxCo);
    }

    settings.push_back(pt::ptree::value_type("runControl", runControl));
}

void saveTurbulenceParameters(const IOdictionary& RASProperties,
                              pt::ptree& settings)
{
    std::string dictName = "kOmegaSSTCoeffs";
    dictionary kOmegaSSTCoeffsDict
            (RASProperties.subOrEmptyDict(dictName));

    scalar alphaK1 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("alphaK1", 0.85034);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_ALPHA_K1, alphaK1);

    scalar alphaK2 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("alphaK2", 1.0);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_ALPHA_K2, alphaK2);

    scalar alphaOmega1 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("alphaOmega1", 0.5);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_ALPHA_OMEGA1, alphaOmega1);

    scalar alphaOmega2 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("alphaOmega2", 0.85616);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_ALPHA_OMEGA2, alphaOmega2);

    scalar gamma1 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("gamma1", 0.5532);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_GAMMA1, gamma1);

    scalar gamma2 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("gamma2", 0.4403);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_GAMMA2, gamma2);

    scalar beta1 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("beta1", 0.0750);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_BETA1, beta1);

    scalar beta2 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("beta2", 0.0828);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_BETA2, beta2);

    scalar betaStar = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("betaStar", 0.09);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_BETASTAR, betaStar);

    scalar a1 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("a1", 0.31);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_A1, a1);

    scalar c1 = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("c1", 10);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_C1, c1);

    scalar Cmu = kOmegaSSTCoeffsDict.lookupOrDefault<scalar>("Cmu", 0.09);
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_CMU, Cmu);

    scalar kappa = 0.41;
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_KAPPA, kappa);

    scalar E = 9.8;
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_E, E);

    scalar ypl = 11.0;
    for (int i = 0; i < 10; i++)
    {
        ypl = ::log(E*ypl)/kappa;
    }
    scalar yPlusLam = ypl;
    settings.put(dictName + "." + speeditcl::io::TURBULENCE_KOMEGASST_YPLUS_LAM, yPlusLam);

}

} // namespace io
} // namespace speeditcl
