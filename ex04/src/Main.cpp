#include <iostream>
#include <sstream>
#include <string>
#include "FileReader.hpp"
#include "Solver.hpp"
#include "Particle.hpp"

//////////
#include <sys/time.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>


int main(int argc, char** argv)
{
    std::string param_file,data_file;
    if(argc == 3){
        param_file=argv[1];
	data_file=argv[2];
    }
    else{
        std::cerr<<"WRONG INPUTS!! You should introduce a paratemer and a data file as arguments!"<<std::endl;
        exit(-1);
		}
  std::cout<<"\nProgram begined!\n";	
  
  FileReader readFile=FileReader();
  readFile.RegisterStringParameter("name","-");
  readFile.RegisterIntParameter("vis_space",1);
  readFile.RegisterDoubleParameter("t_start",1);
  readFile.RegisterDoubleParameter("t_end",1);
  readFile.RegisterDoubleParameter("delta_t",1);
  readFile.RegisterDoubleParameter("x_min",1);
  readFile.RegisterDoubleParameter("y_min",1);
  readFile.RegisterDoubleParameter("z_min",1);
  readFile.RegisterDoubleParameter("x_max",1);
  readFile.RegisterDoubleParameter("y_max",1);
  readFile.RegisterDoubleParameter("z_max",1);
  readFile.RegisterDoubleParameter("r_cut",1);
  readFile.RegisterDoubleParameter("epsilon",1);
  readFile.RegisterDoubleParameter("sigma",1);
  
  //std::cout<<"Before reading\n";
  //readFile.PrintParameters();
  //std::cout<<"\n";
  
  std::string rezfile;
  if(readFile.ReadFile(param_file)!=true){
    std::cout<<"ERROR:: For the parameter file you introduced: '"<<param_file<<"': The file could not be opened!"<<std::endl;
	exit(-1);
  }
  std::cout<<"After reading\n";
  if(readFile.read_data(data_file)!=true){
    std::cout<<"ERROR::For the data file you introduced: '"<<data_file<<"': The file could not be opened!"<<std::endl;
	exit(-1);
  }
  std::cout<<"\nProgram begined!\n";
  
  readFile.PrintParameters();
  std::cout<<"\n";

  
  //std::cout<<"\n Solve!\n";
  
  Solver mysolver = Solver(readFile);

  mysolver.Solve();
  

 
  //mysolver.WriteResults(); //outputs all data to text files for debugging purpose 
  std::cout<<"\nProgram ended!\n";
  return 0;
	
}

