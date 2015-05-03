#include "VTKFileWriter.hpp"


VTKFileWriter::VTKFileWriter(const std::string& problem)
{
    INFO("*****************************************************");
    INFO("Setting output to file: "<<problem);
    INFO("*****************************************************");

    std::string filename = problem;
    filename += ".vtk";
	mStream = fopen(filename.c_str(), "w");

    if ( mStream == NULL )
    {
       std::cout<<"Error opening outputfile";
    }

    fprintf( mStream, "# vtk DataFile Version 3.0\n");
    fprintf( mStream, "LBM - Lid driven cavity\n");
    fprintf( mStream, "ASCII\n");
    fprintf( mStream, "DATASET UNSTRUCTURED_GRID\n");
	
}



VTKFileWriter::~VTKFileWriter( )
{
    fprintf( mStream, "\n" );
    fclose(mStream);
}


void VTKFileWriter::WriteVectorP(const std::vector<std::vector< double > >& vec,
        const  std::string& name)
{
    INFO("*****************************************************");
    INFO("Register Vector: "<<name<<" with dimension "<<vec.size());
    INFO("*****************************************************");

    fprintf( mStream, "POINTS %s double\n", name.c_str());

    switch( vec.size() )
    {
        case 2:
            for(int i=0, end = vec[0].size(); i<end; i++)
            {
                for(int dim=0, iend=vec.size(); dim<iend; dim++)
                {
                    fprintf( mStream, "%f ", vec[dim][i]);
                }
                fprintf( mStream, " 0\n ");
            }
            fprintf( mStream, "\n");
            break;

        case 3:
            for(int i=0, end = vec[0].size(); i<end; i++)
            {
                for(int dim=0, iend=vec.size(); dim<iend; dim++)
                {
                    fprintf( mStream, "%f ", vec[dim][i]);
                }
                fprintf( mStream, "\n ");
            }
            fprintf( mStream, "\n");
            break;
        default:  
	  std::cout<<"Unsupported Dimension "<<vec.size();
    }
	
	fprintf( mStream, "POINT_DATA %s\n",name.c_str());
}

void VTKFileWriter::WriteVector(const std::vector<std::vector< double > >& vec,
        const  std::string& name)
{
    INFO("*****************************************************");
    INFO("Register Vector: "<<name<<" with dimension "<<vec.size());
    INFO("*****************************************************");

    fprintf( mStream, "VECTORS %s double\n", name.c_str());

    switch( vec.size() )
    {
        case 2:
            for(int i=0, end = vec[0].size(); i<end; i++)
            {
                for(int dim=0, iend=vec.size(); dim<iend; dim++)
                {
                    fprintf( mStream, "%f ", vec[dim][i]);
                }
                fprintf( mStream, " 0\n ");
            }
            fprintf( mStream, "\n");
            break;

        case 3:
            for(int i=0, end = vec[0].size(); i<end; i++)
            {
                for(int dim=0, iend=vec.size(); dim<iend; dim++)
                {
                    fprintf( mStream, "%f ", vec[dim][i]);
                }
                fprintf( mStream, "\n ");
            }
            fprintf( mStream, "\n");
            break;
        default:  
	  std::cout<<"Unsupported Dimension "<<vec.size();
    }

}

void VTKFileWriter::WriteScalar(const std::vector<double>& scalar,
        const  std::string& name)
{
    INFO("*****************************************************");
    INFO("Register Scalar: "<<name);
    INFO("*****************************************************");

    fprintf( mStream, "SCALARS %s double 1\n",name.c_str() );
    fprintf( mStream, "LOOKUP_TABLE default\n");

    for(int i=0, end=scalar.size(); i<end; i++)
    {
        fprintf( mStream, "%f\n", scalar[i] );
    }
    fprintf( mStream, "\n");
}


