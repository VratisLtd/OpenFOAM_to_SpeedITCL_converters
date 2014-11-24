#include "OF2SITFCLErrorHandling.h"

#include "OF2SITFCLUtils.h"

#include <sys/stat.h>

namespace speeditcl
{
namespace io
{

void handleError(const std::string& error_string, const std::string& dir_path)
{
    removeDirectory(dir_path.c_str());

    Foam::FatalError << error_string.c_str() ;
    Foam::FatalError << Foam::exit(Foam::FatalError, -1);
}

} // namespace io
} // namespace speeditcl
