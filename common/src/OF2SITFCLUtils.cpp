#include "OF2SITFCLUtils.h"

#include <sys/stat.h>
#include <sys/types.h>

namespace speeditcl
{
namespace io
{

void createDirectory(const char* path)
{
    int status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (0 > status)
    {
        if(errno == EEXIST)
        {
            std::cout << "\nDirectory " << path ; //Foam::FatalError << "Directory " << path ;
            std::cout << " already exists, files will be overwritten. \n" ; //Foam::FatalError << " already exists, please remove it." ;
            //Foam::FatalError << Foam::exit(Foam::FatalError, -1) ;
        }
        if (errno != EEXIST)
        {
            Foam::FatalError << "Can not create directory " << path;
            Foam::FatalError << Foam::exit(Foam::FatalError, -1);
        }
    }
}

void removeDirectory(const char* dir_path)
{
    std::string command = std::string("rm -rf ") + dir_path;
    if(system(command.c_str()))
        return;
}


namespace internal
{

int computeNTimeSteps(Foam::Time& runTime)
{
    dimensionedScalar start_time = runTime.startTime();
    label start_time_index = runTime.startTimeIndex();

    int time_idx = 0;

    while (runTime.loop())
    {
        time_idx++;
    }

    runTime.setTime(start_time, start_time_index);

    return time_idx;
}

void setVecVal(std::vector< std::vector<double>>& vout, scalar val, int idx)
{
    vout[0][idx] = val;
}

void setVecVal(std::vector< std::vector<int>>& vout, int val, int idx)
{
    vout[0][idx] = val;
}

void setVecVal(std::vector< std::vector<double>>& vout, vector val, int idx)
{
    vout[0][idx] = val.x();
    vout[1][idx] = val.y();
    vout[2][idx] = val.z();
}

bool testIO(std::ostream& f)
{
    if (f.good())
    {
    }
    else
    {
        throw std::string("Can not write data to file") ;
    }

    return true;
}

void saveString(std::ostream& f, const std::string& str)
{
    f << str << "\n";
    testIO(f);
}

} // namespace internal
} // namespace io

std::string appToSolver(const std::string& app)
{
    if(app == "simpleFoam")
        return speeditcl::io::SOLVER_SIMPLE;
    else if(app == "icoFoam")
        return speeditcl::io::SOLVER_ICO;
    else if(app == "pisoFoam")
        return speeditcl::io::SOLVER_PISO;

    throw std::string("Unknown app name: ") + app;
}

std::string appToDictName(const std::string& app)
{
    if(app == "simpleFoam")
        return "SIMPLE";
    else if(app == "icoFoam" || app == "pisoFoam")
        return "PISO";

    throw std::string("Unknown app name: ") + app;
}
} // namespace speeditcl
