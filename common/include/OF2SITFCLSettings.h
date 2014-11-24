#ifndef OF_2_SITFCL_SETTINGS_H
#define OF_2_SITFCL_SETTINGS_H

#include "fvCFD.H" /* OpenFoam header */

#include "OF2SITFCLConstants.h"
#include <boost/property_tree/ptree.hpp>

#include <vector>
#include <fstream>
#include <iomanip>


namespace pt = boost::property_tree;


namespace speeditcl
{
namespace io
{

template<class T>
inline std::ostream& settingName(std::ostream& f, T& name)
{
    f << std::left;
    f << std::setw(60) << name;
    return f;
}

void saveIterSolversSettings(pt::ptree& settings, Foam::fvMesh& mesh, const std::string& eqn);

void saveSettingsCommon(Foam::Time &runTime, const std::string &solver_name,
                        pt::ptree& settings, label pRefCell, scalar pRefValue,
                        dimensionedScalar& nu, const bool turbulence);

template <typename TField>
void saveRelaxationParameter(const TField& F,
                             pt::ptree& settings)
{
    std::string filed_name = F.name();
    if(F.mesh().relax(filed_name))
    {
        settings.put(std::string("relaxationFactors.")+filed_name,
                     F.mesh().relaxationFactor(filed_name));
    }
}

void saveTurbulenceParameters(const IOdictionary& RASProperties,
                              pt::ptree& settings);

} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_SETTINGS_H
