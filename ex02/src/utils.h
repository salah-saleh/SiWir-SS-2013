#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include<fstream>
#include <sstream>
#include <map>
#include<vector>
// for posix_memalign:
#include <stdlib.h>

using namespace std;

template <class T>
void allocMem(T ** stencil, const int memSize){ //max size 78658 

    int ok;
    const int align_to = 16;
    ok = posix_memalign((void**)stencil, align_to, (memSize)*sizeof(T));
    if(ok != 0) {std::cout<<"Memory was not allocated!"; exit(EXIT_FAILURE);}

} 

//Initialize array vith value
void initialize(double* vect, const int memSize, const double value); 

//Computes L2 norm 
double L2N(double * const resid, const int NX);

//Global matrices printer  
void MatrixPrint( vector< map<int, double> > &A, vector< map<int, double> > &M, const int NX); 

//Read the file
bool ReadFile(const std::string &name, double* &vertex, int & nv, int* &face, int & nf); 

void faceprint(int * const array, const int n);

void vertexprint(double * const array, const int n);

void solvek(double * k,double * array, const int n,const double delt);

void kprint(double * const k ,double * const array, const int n);

void CreateNeighbourList(int * nneighb, int * &neighb, int * const f, const int nv, const int nf, int & nn);

void getCoord(double * const vertex, const int poz, double& x, double& y);

void getVertexes(int * const face, const int nr, int& v1, int& v2, int& v3);

void getCorners(double * corners,double * vertex,int * face, const int pos);

double comGrad(const double x, const double y);

void copyArray(double * const source, double * const dest, const int n);

void eigprint(double * const eig ,double * const array, const int n);