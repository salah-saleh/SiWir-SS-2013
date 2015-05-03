#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include<vector>
// for posix_memalign:
#include <stdlib.h>
#include "utils.h"

using namespace std;

void global_matrices(vector< map<int, double> > &A, vector< map<int, double> > &M, double * const vertex, int * const face, const int nf);
double power_it(vector< map<int, double> > &A, vector< map<int, double> > &M, double* u, const int nv, const double eps); 
void compRHS(vector< map<int, double> > &M, double* const u, double* rhs, const int nv);
void GS(vector< map<int, double> > &A, double* const u, double* const rhs, const int nv, const double tol);
void compNormU(double* const u, const int nv);
double compLambda(vector< map<int, double> > &A, vector< map<int, double> > &M, double* const u, const int nv);
void RefineGrid( double* &refvertex, int* &refface, double* vertex, int* face, int &nv, int &nf,const int refine);
