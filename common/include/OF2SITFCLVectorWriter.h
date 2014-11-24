#ifndef OF_2_SITFCL_VECTOR_WRITER_H
#define OF_2_SITFCL_VECTOR_WRITER_H

#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>

#include "OF2SITFCLConstants.h"
#include "OF2SITFCLUtils.h"

namespace speeditcl
{
namespace io
{

class OF2SpeeditclVectorWriter
{
public:

    static constexpr int defaultPrecision = 16;

    OF2SpeeditclVectorWriter(bool binary = false, int precision = defaultPrecision)
        : binary(binary), precision(precision)
    {
    }

    template<class T>
    void saveVec(std::ostream& f, const T* data, unsigned int size,
                 const std::string & vec_name)
    {
        saveVecHeader<T>(f, size, 1, vec_name);
        saveVecData(f, data, size);
    }

    template<class T>
    void saveVec(std::ostream & f, const std::vector<T>& vec,
                  const std::string& vec_name)
    {
        saveVec(f, vec.data(), vec.size(), vec_name);
    }

    template<class T>
    void saveVec(std::ostream& f, const std::vector<std::vector<T> >& vec,
                 const std::string& vec_name)
    {
        saveVecHeader<T>(f, vec[0].size(), vec.size(), vec_name);

        for (unsigned int i = 0; i < vec.size(); i++)
            saveVecData(f, &(vec[i][0]), vec[0].size());
    }

private:

    template<class T>
    void saveVecHeader(std::ostream& f, unsigned size, unsigned ncoordinates,
                       const std::string& vec_name)
    {
        internal::saveString(f, "");
        internal::saveString(f, vec_name);
        internal::saveString(f, VECTOR_HEADER);

        f << internal::dataTypeHeader<T>() << "\n";

        internal::saveString(f, binary ? MODE_BINARY : MODE_TEXT);

        f << size << " ";
        f << ncoordinates << "\n";

        internal::testIO(f);
    }

    template<class T>
    void saveVecData (std::ostream & f, const T* data, unsigned size)
    {
        if(binary)
        {
            f.write( reinterpret_cast<const char *>(data), size * sizeof(T));
        }
        else
        {
            f << std::scientific << std::setprecision(precision);
            std::copy(data, data + size, std::ostream_iterator<T>(f, "\n"));
        }

        internal::testIO(f);
    }

    bool binary;
    int precision;
};

} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_VECTOR_WRITER_H
