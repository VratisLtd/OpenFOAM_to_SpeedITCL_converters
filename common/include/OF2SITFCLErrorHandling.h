#ifndef OF_2_SITFCL_ERROR_HANDLING_H
#define OF_2_SITFCL_ERROR_HANDLING_H

#include <string>

namespace speeditcl
{
namespace io
{

void handleError(const std::string& error_string, const std::string& dir_name);

} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_ERROR_HANDLING_H
