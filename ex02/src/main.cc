#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include<vector>
// for posix_memalign:
#include <stdlib.h>
#include "utils.h"
#include "kernels.h"

using namespace std;
// double is defined as double 
int main(int argc, char **argv) {

  double delt,eps;
  int refine;

  if(argc == 4){
    delt= atof(argv[1]);
    eps = atof(argv[2]); 
    refine=atoi(argv[3]);
  }

  else{
    if(argc == 3){
      delt = atof(argv[1]);
      eps = atof(argv[2]);
      refine=0;
    }
    else{
      cout<<"WRONG INPUTS!!"<<endl;
      exit(-1);
    }   
  }

  int nv=0; 
  int nf=0;
  int nn=0;

  double * vertex = NULL;
  int * face = NULL;
  double * kFunc = NULL; 
  double * refvertex = NULL;
  int * refface = NULL;

  // Neighbours of all the vertcies 
  int * neighbList = NULL; 

  // For every vertex the first element has the number of neighbours and the second one holds the position of the first neighb. in neighb. list   
  int * numNeighbList = NULL; 

  //######################## Read file ################################
  if(ReadFile("../unit_circle.txt", vertex, nv, face, nf)==false){

    cout<<"The input file is missing! Check that 'unit_circle.txt' is here : ./ex02_group13  "<<"\n";
    cout<<"Trying to read 'unit_circle_fine.txt' !"<<"\n";
    
    if(ReadFile("../unit_circle_fine.txt", vertex, nv, face, nf)==false){
      cout<<"No input file! !Terminating program! "<<"\n";
      return 0;
    }

  }

  //######################## Refine grid #################################
  for(int i=0;i<refine;i++){

    RefineGrid( refvertex, refface, vertex, face, nv, nf, refine);

    double * tempv=vertex; vertex=refvertex;
    int * tempf=face; face=refface;

    free(tempv);
    free(tempf);
  }

  //###################### Compute k(x,y) ################################
  allocMem(&kFunc, nv);
  solvek(kFunc, vertex, nv, delt);
  
  // Computes variable coefficient that depends on the refractive index of the 
  // material and the wavelength of the propagating beam. 
  kprint(kFunc , vertex, nv);

  //#################### Create neighbours list ##########################
  allocMem(&numNeighbList, 2*nv);
  CreateNeighbourList(numNeighbList, neighbList, face, nv, nf, nn);

  //#################### Build global matrices ########################### 
  vector< map<int, double> > A(nv);
  vector< map<int, double> > M(nv);

  global_matrices( A, M, vertex, face, nf);
  MatrixPrint( A, M, nv);  

  //#################### Inverse power iteration ##########################
  double * u = NULL;
  allocMem(& u, nv);
  cout<<"lambda= "<<power_it( A, M, u, nv, eps)<<endl;
  eigprint(u , vertex, nv);

  free(vertex);
  free(face);
  free(kFunc);
  free(u);
  free(neighbList);
  free(numNeighbList);
  return 0;
}

