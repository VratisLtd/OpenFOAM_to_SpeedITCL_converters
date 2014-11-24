#ifndef OF_2_SITFCL_UTILS_H
#define OF_2_SITFCL_UTILS_H

#include "fvCFD.H" /* OpenFoam header */

#include "OF2SITFCLConstants.h"

#include <vector>
#include <fstream>

namespace speeditcl
{
namespace io
{

/** Warning: those 2 functions use system-dependent */
void createDirectory (const char* path);
void removeDirectory(const char* dir_path);

namespace internal
{

int computeNTimeSteps(Foam::Time & runTime);

/**
  * vout - result
  * val - data source
  */
void setVecVal(std::vector<std::vector<double>>& vout, scalar val, int idx);
void setVecVal(std::vector<std::vector<int>>& vout, int val, int idx);
void setVecVal(std::vector<std::vector<double>>& vout, vector val, int idx);

bool testIO(std::ostream& f);
void saveString(std::ostream& f, const std::string& header);

template<class Type>
std::vector<Type> vectorReorder(const std::vector<Type>& in_vec,
                                const std::vector<int>& reorder_table)
{
    std::vector<Type> vres;
    vres.resize(in_vec.size());

    for (unsigned i = 0; i < in_vec.size(); i++)
    {
        vres[i] = in_vec[reorder_table[i]];
    }

    return vres;
}

template<class Type>
std::vector<std::vector<Type>> vector_reorder(
                                   const std::vector<std::vector<Type>>& in_vec,
                                   const std::vector<int>& reorder_table)
{
    std::vector< std::vector<Type> > vres;
    vres.resize(in_vec.size());

    for(unsigned i = 0; i < in_vec.size(); i++)
    {
        vres[i] = vectorReorder(in_vec[i], reorder_table);
    }

    return vres;
}


template<class T>
inline std::string dataTypeHeader()
{
    throw std::string("Unsupported data type in vector : \"")
            + typeid(T).name()
            + "\", the only supported data types are: double int.";
}


template<>
inline std::string dataTypeHeader<double>()
{
    return DATA_DOUBLE;
}

template<>
inline std::string dataTypeHeader<int>()
{
    return DATA_INT;
}

} // namespace internal
} // namespace io

std::string appToSolver(const std::string& app);
std::string appToDictName(const std::string& app);
} // namespace speeditcl

#endif // OF_2_SITFCL_UTILS_H
