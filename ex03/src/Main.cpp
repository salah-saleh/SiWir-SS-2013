#include <iostream>
#include <sstream>
#include <string>
#include "FileReader.hpp"
#include "Solver.hpp"

//////////
#include <sys/time.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>



void timing(double* wcTime, double* cpuTime);
/////////

int main(int argc, char** argv)
{
    std::string infile;
    if(argc == 2){
        infile=argv[1];
    }
    else{
        std::cout<<"WRONG INPUTS!!"<<std::endl;
        exit(-1);
		}
  std::cout<<"\nProgram begined!\n";	
  
  FileReader readFile=FileReader();
  readFile.RegisterDoubleParameter("omega",1);
  readFile.RegisterIntParameter("timesteps",1);
  readFile.RegisterIntParameter("sizex",1);
  readFile.RegisterIntParameter("sizey",1);
  readFile.RegisterStringParameter("vtk_file","-");
  readFile.RegisterIntParameter("vtk_step",1);
  readFile.RegisterStringParameter("geometry","-");
  //std::cout<<"Before reading\n";
  //readFile.PrintParameters();
  //std::cout<<"\n";
  
  std::string rezfile;
  if(readFile.ReadFile(infile)!=true){
    std::cout<<"ERROR:: You introduced: '"<<infile<<"': The file could not be opened!"<<std::endl;
	exit(-1);
  }
  std::cout<<"After reading\n";
  if(readFile.GetStringParameter("geometry")!="-")
     readFile.readPGM(readFile.GetStringParameter("geometry"));
  readFile.PrintParameters();
  std::cout<<"\n";

  std::cout<<"\n Solve!\n";
  
  Solver mysolver = Solver(readFile);
 
  double wct_start,wct_end,cput_start,cput_end;
  double dauer=0.0;

  timing(&wct_start, &cput_start);
  if(readFile.GetStringParameter("geometry")!="-")
    mysolver.Solve_geometry();
  else
    mysolver.Solve();
  
  timing(&wct_end, &cput_end);
  dauer=wct_end-wct_start;

  

  double mlups =(double)readFile.GetIntParameter("timesteps") * readFile.GetIntParameter("sizey") * readFile.GetIntParameter("sizex") / dauer / 1.e6;  
  std::cout<<"MLUPs: "<< mlups<<std::endl;
  //mysolver.WriteResults(); //outputs all data to text files for debugging purpose 
  std::cout<<"\nProgram ended!\n";
  return 0;
	
}

///////////
void timing(double* wcTime, double* cpuTime)
{
   struct timeval tp;
   struct rusage ruse;

   gettimeofday(&tp, NULL);
   *wcTime=(double) (tp.tv_sec + tp.tv_usec/1000000.0); 
  
   getrusage(RUSAGE_SELF, &ruse);
   *cpuTime=(double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec / 1000000.0);
}
//////////
