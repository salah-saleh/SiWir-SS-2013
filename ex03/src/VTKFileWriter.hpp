#ifndef VTKFILEWRITER_HH
#define VTKFILEWRITER_HH

#include "Debug.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>

class VTKFileWriter
{
    public:

        VTKFileWriter(const std::vector< int > &dimension,
                const std::vector< float >& length, const std::string &problem);

        ~VTKFileWriter();

        void WriteVector(const std::vector< std::vector< double > > &vec,
                const std::string &name);

        void WriteScalar(const std::vector<double> &scalar, const std::string &name);

    private:

        FILE *mStream;
        unsigned int mDataCount;
        std::vector< int > mDim;
};
#endif // VTKFILEWRITER_HH

