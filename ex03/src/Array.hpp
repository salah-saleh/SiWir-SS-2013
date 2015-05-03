#ifndef ARRAY_HH
#define ARRAY_HH

#include "Types.hpp"
#include "FileReader.hpp"
#include <iostream>
#include <fstream>
#include <sstream>


class Array {

    public:

    // Constructor 1D
    Array(int x) 
	{
		Nx=x;
		Ny=0;
		Nz=0;
		myvalues = new double[Nx];
        }


    // Constructor 2D
    Array(int x,int y) 
	{
		Nx=x;
		Ny=y;
		Nz=0;
		myvalues = new double[(Nx)*(Ny)];
	}

    // Constructor 3D
    Array(int x,int y,int z) 
	{
		Nx=x;
		Ny=y;
		Nz=z;
		myvalues = new double[Nx*Ny*Nz];
	}

    //If you need the following depends on your implementation
    // Destructor
    // Default implementation is just fine here
    ~Array()
	{
	  delete[] myvalues;
	}

    //copy Constructor
    // Use compiler generated default implementation, which is just fine here
    //  Array(const Array& s);

    //assignment  Operator 
    // Use compiler generated default implementation, which is just fine here
    // Array& operator= (const Array& s);

    // Operator() 1D
    inline real& operator ()(int i)
    { 
    	//myvalues[i]+=myvalues[i];
  		return myvalues[i];
    }

    // Operator() 2D
    inline real& operator ()(int i,int j)
    {
  		//myvalues[i*Nx+j]+=myvalues[i*Nx+j];
  		return myvalues[j*Nx+i];
    }

    // Operator() 3D
    inline real& operator ()(int i, int j, int k)
    { 
 		 //myvalues[i*Nx+j+k*Nx*Ny]+=myvalues[i*Nx+j+k*Nx*Ny];
  		return myvalues[j*Nx+i+k*Nx*Ny];
    }


    //initialize the whole array with a constant value
    inline void fill(const real in)
    {
		if(Nx==0) 
			std::cout<<"Error: The dimension is '0'(zero)!"<<"\n";
		else if(Ny==0 && Nz==0)
		{
		 	for(int i=0;i<Nx;i++)myvalues[i]=in;
		}
		else if(Nz==0) 
		{	
			for(int i=0;i<Nx*Ny;i++)myvalues[i]=in;	
		}
		else 
		{
			for(int i=0;i<Nx*Ny*Nz;i++) myvalues[i]=in;
		}
    }
    //return the dimension of the array
    inline int getDimension() const
    { 
    	if (Nx==0){std::cout<<"Error: The dimension is '0'(zero)!"<<"\n";return 0;}
		else if (Ny==0) return 1;
			else if(Nz==0) return 2;
				else return 3;
    }
			 	
	//return total size of the array
    inline int getSize() const
    { 
    	if (Nx==0){std::cout<<"Error: The dimension is '0'(zero)!"<<"\n";return 0;}
		else if (Ny==0) return Nx;
			else if(Nz==0) return Nx*Ny;
				else return Nx*Ny*Nz;
    }


    // Print the whole array (for debugging purposes)
    inline void print()
    {	

		if(Nx==0) 
			std::cout<<"Error:The array is empty! The dimension is '0'(zero)!"<<"\n";
		else if(Ny==0 && Nz==0)
		{
		 	for(int i=0;i<Nx;i++)std::cout<<myvalues[i]<<" ";
		}
		else if(Nz==0) 
		{	for(int j=Ny-1;j>=0;j--)
			{
			  for(int i=0;i<Nx;i++)
			    std::cout<<myvalues[j*Nx+i]<<" ";
			  std::cout<<"\n";
			}
		}
		else 
		{
			for(int i=0;i<Nx*Ny*Nz;i++)std::cout<<myvalues[i]<<" ";
		}
    }
        
    int fprint(const char* fname)
    {	
		std::ofstream myfile;
	    myfile.open(fname);
       if (myfile.is_open())
       {
            if(Nx==0) 
		  	  	myfile<<"Error:The array is empty! The dimension is '0'(zero)!"<<"\n";
		  	else if(Ny==0 && Nz==0)
		  	{
		 		for(int i=0;i<Nx;i++)myfile<<myvalues[i]<<" ";
		  	}
		  	else if(Nz==0) 
		  	{	
				for(int j=Ny-1;j>=1;j--)
				{
			  		for(int i=1;i<Nx-1;i++)
			  		{
						myfile<<j<<" "<<i<<" "<<myvalues[j*Nx+i]<<"\n";
		      		}
			   		myfile<<"\n";
				}
		  	}
		  	else 
		  	{
				for(int i=0;i<Nx*Ny*Nz;i++)myfile<<myvalues[i]<<" ";
		  	}
		  	
		  	myfile.close();
		  	return 1;
		}
		else
		  	return 0;   	  
	}
	
        
    private:
    	
    	int Nx,Ny,Nz;
    	double *myvalues;
    	


};
#endif //ARRAY_HH

