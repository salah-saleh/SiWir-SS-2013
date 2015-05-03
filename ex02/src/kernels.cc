#include "./Source/Colsamm.h"
#include "kernels.h"

using namespace ::_COLSAMM_;
using namespace std;

//////////////////////////////////////////////////////////////////
void global_matrices(vector< map<int, double> > &A, vector< map<int, double> > &M, double * const vertex, int * const face, const int nf){   

   for(int t=0;t<nf;t++){
      double * corners = NULL;
      allocMem(& corners, 6);   

      //an example for the 8th face (number 7 in the code)
      getCorners(corners, vertex, face, t); 
      int v0=0;
      int v1=0;
      int v2=0;

      getVertexes(face, t, v0, v1, v2);
      int v[3];
      v[0]=v0;v[1]=v1;v[2]=v2;
      ELEMENTS::Triangle my_triang; 
      
      my_triang(corners);

      free(corners);

      vector< vector <double> > M_l;
      vector< vector <double> > A_l;

      M_l = my_triang.integrate(v_() * w_());
      A_l = my_triang.integrate( grad(v_()) * grad(w_()) - func<double>(comGrad) * v_() * w_() );
      
      pair<map<int,double>::iterator,bool> ret;
      for(int i=0; i<3; ++i)
         for(int j=0; j<3; ++j){		
            ret=A[v[i]].insert(pair<int, double>(v[j], A_l[i][j]));
            if (ret.second==false)
               A[v[i]][v[j]] += A_l[i][j];

         }

         for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j){		
               ret=M[v[i]].insert(pair<int, double>(v[j], M_l[i][j]));
               if (ret.second==false)
                  M[v[i]][v[j]] += M_l[i][j];

            }	
   }
}

///////////////////////////////////////////////////////////////////
double power_it(vector< map<int, double> > &A, vector< map<int, double> > &M, double* u, const int nv, const double tol){

   double l=0.1;
   double l_o=10;

   initialize(u, nv, 1/sqrt(nv));

   double * rhs = NULL;
   allocMem(& rhs, nv);

   while(abs((l-l_o)/l_o) > 10.e-10){

      l_o = l;
      compRHS(M, u, rhs, nv);      
      GS(A, u, rhs, nv, tol);
      compNormU(u, nv);
      l = compLambda(A, M, u, nv);

   }

   free(rhs);
   return l;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void compRHS(vector< map<int, double> > &M, double* const u, double* rhs, const int nv){

   map<int,double>::iterator itgs;

   for(int i=0; i<nv; ++i){
      double temp=0.0;
      for (itgs=M[i].begin(); itgs!=M[i].end(); ++itgs)
         temp+=(itgs->second)*u[itgs->first];
      rhs[i]=temp;
   }

}

//////////////////////////////////////////////////////////////////////////////////////////////
void GS(vector< map<int, double> > &A, double* const u, double* const rhs, const int nv, const double tol){

      //double * u_o = NULL;
      //allocMem(& u_o, nv);

   double * err = NULL;
   allocMem(& err, nv);

   map<int,double>::iterator itgs;     
   double stop=1.0;
   int maxiter=10000;
   int nriter=0;

   // solve Gauss Seidel
   while((stop>tol)&(nriter<maxiter)){                          

      for (int i=0;i<nv;i++){
         double s=0.0;
         double temp=0.0;

         for (itgs=A[i].begin(); itgs!=A[i].end(); ++itgs){

            if(i==itgs->first) temp=itgs->second; 
            else s=s+(itgs->second)*u[itgs->first];
         }
         u[i]=(rhs[i]-s)/temp;  
      }

      for(int i=0; i<nv; ++i){

         double temp=0.0;

         for (itgs=A[i].begin(); itgs!=A[i].end(); ++itgs)
            temp+=(itgs->second)*u[itgs->first];

         err[i]=temp-rhs[i];
      }

      st
      op=L2N(err ,nv);
      nriter++;
   }

      free(err);
}

//////////////////////////////////////////////////////////////////////////////////////////////
void compNormU(double* const u, const int nv){

   double tempnorm;
   tempnorm=L2N(u, nv);
   
   for(int i=0; i<nv; ++i){  
      u[i]=u[i]/tempnorm;
   }

}

//////////////////////////////////////////////////////////////////////////////////////////////
double compLambda(vector< map<int, double> > &A, vector< map<int, double> > &M, double* const u, const int nv){

   double l=0.0;	
   map<int,double>::iterator itgs;
   double denum=0.0;

   for(int i=0; i<nv; ++i){
      double temp=0.0;

      for (itgs=M[i].begin(); itgs!=M[i].end(); ++itgs)
         temp+=(itgs->second)*u[itgs->first];

      denum+=u[i]*temp;
   }

   double num=0.0;

   for(int i=0; i<nv; ++i){

      double temp=0.0;

      for (itgs=A[i].begin(); itgs!=A[i].end(); ++itgs)
         temp+=(itgs->second)*u[itgs->first];
      
      num+=u[i]*temp;
   }

   return l = num/denum;
}


//////////////////////////////////////////////////////////////////////////////////////////////
void RefineGrid( double* &refvertex, int* &refface, double* vertex, int* face, int &nv, int &nf,const int refine){

   vector <double> nvertex;
   vector <int> nface;

   double * fine_cor = NULL;
   int * fine_pos = NULL;
   allocMem(& fine_cor, 6);
   allocMem(& fine_pos, 3);
   
   for (int n=0; n<nv*2; n++)
     nvertex.push_back(vertex[n]);

   for(int j=0;j<nf;j++){

      fine_cor[0]=(vertex[face[j*3]*2]+vertex[face[j*3+1]*2])/2.0;
      fine_cor[1]=(vertex[face[j*3]*2+1]+vertex[face[j*3+1]*2+1])/2.0;
      fine_cor[2]=(vertex[face[j*3]*2]+vertex[face[j*3+2]*2])/2.0;
      fine_cor[3]=(vertex[face[j*3]*2+1]+vertex[face[j*3+2]*2+1])/2.0;
      fine_cor[4]=(vertex[face[j*3+1]*2]+vertex[face[j*3+2]*2])/2.0;
      fine_cor[5]=(vertex[face[j*3+1]*2+1]+vertex[face[j*3+2]*2+1])/2.0;

      for(int k=0;k<3;k++){
         bool found=false;
         double x=fine_cor[k*2];
         double y=fine_cor[k*2+1];

         for(unsigned int t=nv;t<nvertex.size()/2;t++){  
          if(x==nvertex[t*2]) 
            if(y==nvertex[t*2+1]){ 
             found=true;
             fine_pos[k]=t;
             break;
          }
            }//t 
            
            if(found==false) {
             fine_pos[k]=nvertex.size()/2;
             nvertex.push_back(fine_cor[k*2]);
             nvertex.push_back(fine_cor[k*2+1]);
          }

      } 
      
      nface.push_back(face[j*3]);   nface.push_back(fine_pos[0]);  nface.push_back(fine_pos[1]);
      nface.push_back(face[j*3+1]); nface.push_back(fine_pos[0]);  nface.push_back(fine_pos[2]);
      nface.push_back(face[j*3+2]); nface.push_back(fine_pos[1]);  nface.push_back(fine_pos[2]);
      nface.push_back(fine_pos[0]); nface.push_back(fine_pos[1]);  nface.push_back(fine_pos[2]);

   }
   
   free(fine_pos);
   free(fine_cor);

   nf=nface.size()/3; 
   nv=nvertex.size()/2;
   allocMem(& refvertex, nv*2);
   allocMem(& refface, nf*3);

   for (int n=0; n<nv*2; n++)
     refvertex[n]=nvertex[n];
  
   for (int n=0; n<nf*3; n++)
     refface[n]=nface[n];
}















