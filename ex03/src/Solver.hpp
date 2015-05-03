#ifndef SOLVER_HH
#define SOLVER_HH

#include "FileReader.hpp"
#include "VTKFileWriter.hpp"
#include "Array.hpp"
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>

class Solver
{
  public:

	Solver(FileReader& parameters)
	{ 
	  omg=parameters.GetDoubleParameter("omega");
	  tstep=parameters.GetIntParameter("timesteps");
	  sx=parameters.GetIntParameter("sizex");	
	  sy=parameters.GetIntParameter("sizey");
	  vtks=parameters.GetIntParameter("vtk_step");
	  vtkname=parameters.GetStringParameter("vtk_file");
	  vtknam = vtkname.substr(0, vtkname.find("."));
      if(parameters.GetStringParameter("geometry")!="-")
	      flag = parameters.flag;
	  else {
	  	  flag=new Array(sx+2,sy+2);
          (*flag).fill(0.0);
	  }
	  init_data_struct();
	}	
	
	 //separate solve function for the case of "geometry" domain
	 bool Solve_geometry() 
   {
     int t=0,n=0;
	 std::stringstream nth;
	 moving_bound();
     //(*flag).print();
	 while(t<=tstep)
	 {  
	    correct_bound_geom();
		stream();
		comp_rho();
		comp_velocity();	
		collide();
		
		if(vtks!=0) if(t%vtks==0){
		  nth.str("");
		  nth<<n;
		  vtkPrint(std::string(vtknam)+nth.str());
		}

		t++;n++;
		//std::cout<<sum_f()<<std::endl;
	    //check_f();
	 } 
	  
      return 1;
   }
	// Solves the L-B eq
   bool Solve()
   {
	 int t=0,n=0;
	 std::stringstream nth;
	 moving_bound();
     //init_test();
	    
	 while(t<=tstep)
	 {  
	    correct_bound_geom();
		stream();
		comp_rho();
		comp_velocity();	
		collide();
		
		if(t%vtks==0){
		  nth.str("");
		  nth<<n;
		  vtkPrint(std::string(vtknam)+nth.str());
		}
		t++;n++;
	 }
	 
	return 0;	 
    }
	
    
    void comp_rho()
    { 
      for(int j=sy;j>=1;j--)
	   for(int i=sx;i>=1;i--)
	    {
	      if((*flag)(i,j)!=1)
		    (*rho)(i,j)=(*fc)(i,j)+(*fn)(i,j)+(*fne)(i,j)+(*fe)(i,j)+(*fse)(i,j)+(*fs)(i,j)+(*fsw)(i,j)+(*fw)(i,j)+(*fnw)(i,j);
	    }
    }
	
	double sum_f()
    {  
	  double sum=0;
      for(int i=1;i<=sx;i++)
	   for(int j=1;j<=sy;j++)
	    {
		  if((*flag)(i,j)!=1)
	        sum+=(*fc)(i,j)+(*fn)(i,j)+(*fne)(i,j)+(*fe)(i,j)+(*fse)(i,j)+(*fs)(i,j)+(*fsw)(i,j)+(*fw)(i,j)+(*fnw)(i,j);
	    }
	  return sum;
    }
	
	void check_f()
    {  
      for(int i=1;i<=sx;i++)
	   for(int j=1;j<=sy;j++)
	    {
	      if(((*fc)(i,j)>0.5)||((*fc)(i,j)<0)) {std::cout<<"error:: fc: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fn)(i,j)>0.5)||((*fn)(i,j)<0)) {std::cout<<"error:: fe: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fne)(i,j)>0.5)||((*fne)(i,j)<0)) {std::cout<<"error:: fne: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fe)(i,j)>0.5)||((*fe)(i,j)<0)) {std::cout<<"error:: fe: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fse)(i,j)>0.5)||((*fse)(i,j)<0)) {std::cout<<"error:: fse: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fs)(i,j)>0.5)||((*fs)(i,j)<0)) {std::cout<<"error:: fs: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fsw)(i,j)>0.5)||((*fsw)(i,j)<0)) {std::cout<<"error:: fsw: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fw)(i,j)>0.5)||((*fw)(i,j)<0)) {std::cout<<"error:: fw: "<<i<<" "<<j<<std::endl;exit(-1);}
		  if(((*fnw)(i,j)>0.5)||((*fnw)(i,j)<0)) {std::cout<<"error:: fnw: "<<i<<" "<<j<<std::endl;exit(-1);}
		  
	    }
    }
    
    void comp_velocity()
    {
      for(int j=sy;j>=1;j--)
	   for(int i=sx;i>=1;i--)
	    {
		 if((*flag)(i,j)!=1){
	      (*ux)(i,j)=((*fne)(i,j)+(*fe)(i,j)+(*fse)(i,j)-(*fsw)(i,j)-(*fw)(i,j)-(*fnw)(i,j));//(*rho)(i,j);
	      (*uy)(i,j)=((*fnw)(i,j)+(*fn)(i,j)+(*fne)(i,j)-(*fse)(i,j)-(*fs)(i,j)-(*fsw)(i,j));//(*rho)(i,j);
	      }
		}
    }
    
    void collide()
    {  
      for(int j=sy;j>=1;j--)
	   for(int i=sx;i>=1;i--)
	   { 
	    if((*flag)(i,j)!=1){
	       double vx=(*ux)(i,j);
	       double vy=(*uy)(i,j);
	       double u2 = vx*vx+vy*vy;
	       double rho_=(*rho)(i,j);
	       double feqc =(4.0/9.0 )* (rho_                                                           );
	       double feqn =(1.0/9.0 )* (rho_ +3*         vy   +4.5*vy*vy                       - 1.5*u2);
	       double feqne=(1.0/36.0)* (rho_ +3*  (vx+   vy)  +4.5*(vx+   vy)*(vx+   vy)       - 1.5*u2);
	       double feqe =(1.0/9.0 )* (rho_ +3*   vx         +4.5*vx*vx                       - 1.5*u2);
	       double feqse=(1.0/36.0)* (rho_ +3*  (vx+ (-vy)) +4.5*(vx+ (-vy))*(vx+ (-vy))     - 1.5*u2);
	       double feqs =(1.0/9.0 )* (rho_ +3*       (-vy)  +4.5*vy*vy                       - 1.5*u2);
	       double feqsw=(1.0/36.0)* (rho_ +3*((-vx)+(-vy)) +4.5*((-vx)+(-vy))*((-vx)+(-vy)) - 1.5*u2);
	       double feqw =(1.0/9.0 )* (rho_ +3* (-vx)        +4.5*vx*vx                       - 1.5*u2);
	       double feqnw=(1.0/36.0)* (rho_ +3*((-vx)+  vy)  +4.5*((-vx)+  vy)*((-vx)+  vy)   - 1.5*u2);
	       
	       (*fc)(i,j)=(*fc)(i,j)-omg*((*fc)(i,j)-feqc);
	       (*fn)(i,j)=(*fn)(i,j)-omg*((*fn)(i,j)-feqn);
	       (*fs)(i,j)=(*fs)(i,j)-omg*((*fs)(i,j)-feqs);
	       (*fe)(i,j)=(*fe)(i,j)-omg*((*fe)(i,j)-feqe);
	       (*fw)(i,j)=(*fw)(i,j)-omg*((*fw)(i,j)-feqw);
	       (*fse)(i,j)=(*fse)(i,j)-omg*((*fse)(i,j)-feqse);
	       (*fsw)(i,j)=(*fsw)(i,j)-omg*((*fsw)(i,j)-feqsw);
	       (*fne)(i,j)=(*fne)(i,j)-omg*((*fne)(i,j)-feqne);
	       (*fnw)(i,j)=(*fnw)(i,j)-omg*((*fnw)(i,j)-feqnw);
		}
	   }
    }
    
    //Stream
    void stream()
    {
        //north//
	for(int i=1;i<=sx;i++)
	  for(int j=sy;j>=1;j--)
	  {
	   if((*flag)(i,j)!=1)
		(*fn)(i,j)=(*fn)(i,j-1);
	  }
	  
	//north-east//
	for(int i=sx;i>=1;i--)
	  for(int j=sy;j>=1;j--)
	  {
	    if((*flag)(i,j)!=1) 
	     (*fne)(i,j)=(*fne)(i-1,j-1);
	  }
	  
	//north-west//
	for(int i=1;i<=sx;i++)
	  for(int j=sy;j>=1;j--)
	  { 
	    if((*flag)(i,j)!=1)
	     (*fnw)(i,j)=(*fnw)(i+1,j-1);
	  }
	  
	//south//
	for(int i=1;i<=sx;i++)
	  for(int j=1;j<=sy;j++)
	  {
	    if((*flag)(i,j)!=1)
	     (*fs)(i,j)=(*fs)(i,j+1);
	  }
	  
	//south-west//
	for(int j=1;j<=sy;j++)
	  for(int i=1;i<=sx;i++)
	  {
	    if((*flag)(i,j)!=1)
	     (*fsw)(i,j)=(*fsw)(i+1,j+1);
	  }
	  
	//south-east//
	for(int j=1;j<=sy;j++)
	  for(int i=sx;i>=1;i--)
	  {
	    if((*flag)(i,j)!=1)
		 (*fse)(i,j)=(*fse)(i-1,j+1);
	  }
	  
	//east//
	for(int j=sy;j>=1;j--)
	  for(int i=sx;i>=1;i--)
	  {
	    if((*flag)(i,j)!=1)
		 (*fe)(i,j)=(*fe)(i-1,j);
	  }
	  
	//west//
	for(int i=1;i<=sx;i++)
	  for(int j=1;j<=sy;j++)
	  {
	    if((*flag)(i,j)!=1)
	     (*fw)(i,j)=(*fw)(i+1,j);
	  }
    }

     //Correct the boundary values
    void correct_bound_geom()
	{ 
      for(int i=1;i<=sx;i++)
	    for(int j=1;j<=sy;j++)
	    {
		  if((*flag)(i,j)==1){ 
		    if((*flag)(i,j+1)!=1)    (*fn)(i,j)=(*fs)(i,j+1);
			if((*flag)(i-1,j+1)!=1)  (*fnw)(i,j)=(*fse)(i-1,j+1);
			if((*flag)(i+1,j+1)!=1)  (*fne)(i,j)=(*fsw)(i+1,j+1);
	        if((*flag)(i-1,j)!=1)    (*fw)(i,j)=(*fe)(i-1,j);
			if((*flag)(i-1,j-1)!=1)  (*fsw)(i,j)=(*fne)(i-1,j-1);
	        if((*flag)(i,j-1)!=1)    (*fs)(i,j)=(*fn)(i,j-1);
			if((*flag)(i+1,j-1)!=1)  (*fse)(i,j)=(*fnw)(i+1,j-1);
			if((*flag)(i+1,j)!=1)    (*fe)(i,j)=(*fw)(i+1,j);
		  }
		}
		
		//corners
      if((*flag)(1,1)!=1)  (*fne)(0,0)=(*fsw)(1,1);
      if((*flag)(sx,1)!=1) (*fnw)(sx+1,0)=(*fse)(sx,1);
      if((*flag)(1,sy)!=1) (*fse)(0,sy+1)=(*fnw)(1,sy)+(*ux)(0,sy+1)/6;
      if((*flag)(sx,sy)!=1)(*fsw)(sx+1,sy+1)=(*fne)(sx,sy)-(*ux)(sx+1,sy+1)/6;
	
	  //(north) 
      for(int i=1;i<=sx;i++)
      {
	   if((*flag)(i,sy)!=1)  (*fs)(i,sy+1)=(*fn)(i,sy);//-2.0*(*uy)(j,sy+1)/3;
	   if((*flag)(i+1,sy)!=1)(*fse)(i,sy+1)=(*fnw)(i+1,sy)+(*ux)(i,sy+1)/6;
	   if((*flag)(i-1,sy)!=1)(*fsw)(i,sy+1)=(*fne)(i-1,sy)-(*ux)(i,sy+1)/6;
      }
      //(south) 
      for(int i=1;i<=sx;i++)
      {
	    if((*flag)(i,1)!=1)  (*fn)(i,0)=(*fs)(i,1);
	    if((*flag)(i+1,1)!=1)(*fne)(i,0)=(*fsw)(i+1,1);
	    if((*flag)(i-1,1)!=1)(*fnw)(i,0)=(*fse)(i-1,1);
      }
      //(east) 
      for(int j=1;j<=sy;j++)
      {
	    if((*flag)(sx,j)!=1)  (*fw)(sx+1,j)=(*fe)(sx,j);
	    if((*flag)(sx,j-1)!=1)(*fsw)(sx+1,j)=(*fne)(sx,j-1);
	    if((*flag)(sx,j+1)!=1)(*fnw)(sx+1,j)=(*fse)(sx,j+1);
      }
      //(west) 
      for(int j=1;j<=sy;j++)
      {
	    if((*flag)(1,j)!=1)  (*fe)(0,j)=(*fw)(1,j);
	    if((*flag)(1,j-1)!=1)(*fse)(0,j)=(*fnw)(1,j-1);
	    if((*flag)(1,j+1)!=1)(*fne)(0,j)=(*fsw)(1,j+1);
      }
	}
     //Correct the boundary values
    void correct_bound()
    { 
      //(north) 
      for(int i=1;i<=sx;i++)
      {
	(*fs)(i,sy+1)=(*fn)(i,sy);//-2.0*(*uy)(j,sy+1)/3;
	(*fse)(i,sy+1)=(*fnw)(i+1,sy)+(*ux)(i,sy+1)/6;
	(*fsw)(i,sy+1)=(*fne)(i-1,sy)-(*ux)(i,sy+1)/6;
      }
      //(south) 
      for(int i=1;i<=sx;i++)
      {
	(*fn)(i,0)=(*fs)(i,1);
	(*fne)(i,0)=(*fsw)(i+1,1);
	(*fnw)(i,0)=(*fse)(i-1,1);
      }
      //(east) 
      for(int j=1;j<=sy;j++)
      {
	(*fw)(sx+1,j)=(*fe)(sx,j);
	(*fsw)(sx+1,j)=(*fne)(sx,j-1);
	(*fnw)(sx+1,j)=(*fse)(sx,j+1);
      }
      //(west) 
      for(int j=1;j<=sy;j++)
      {
	(*fe)(0,j)=(*fw)(1,j);
	(*fse)(0,j)=(*fnw)(1,j-1);
	(*fne)(0,j)=(*fsw)(1,j+1);
      }
      
      //corners
      (*fne)(0,0)=(*fsw)(1,1);
      (*fnw)(sx+1,0)=(*fse)(sx,1);
      (*fse)(0,sy+1)=(*fnw)(1,sy)+(*ux)(0,sy+1)/6;
      (*fsw)(sx+1,sy+1)=(*fne)(sx,sy)-(*ux)(sx+1,sy+1)/6;

    }
    
    //Set problem specific boundary conditions 
    void moving_bound()
    {  
      for (int i=0;i<=sx+1;i++)
           (*ux)(i,sy+1)=0.08;
    }


   void WriteResults() //print 
   {
	  (*fc).fprint("ResultsTXT/fc.txt");
	  (*fn).fprint("ResultsTXT/fn.txt");
	  (*fe).fprint("ResultsTXT/fe.txt");
	  (*fs).fprint("ResultsTXT/fs.txt");
	  (*fw).fprint("ResultsTXT/fw.txt");
	  (*fne).fprint("ResultsTXT/fne.txt");
	  (*fnw).fprint("ResultsTXT/fnw.txt");
	  (*fse).fprint("ResultsTXT/fse.txt");
	  (*fsw).fprint("ResultsTXT/fsw.txt");
	  (*ux).fprint("ResultsTXT/ux.txt");
	  (*uy).fprint("ResultsTXT/uy.txt");
   }
   
   void init_test()
   {
     for(int j=sy+1;j>=0;j--)
	   for(int i=sx+1;i>=0;i--)
	   { 
	     (*ux)(i,j)=0.01;
	     (*uy)(i,j)=0.01;
	     double vx=(*ux)(i,j);
	     double vy=(*uy)(i,j);
	     double u2 = vx*vx+vy*vy;
	     double rho_=1;
	     (*fc)(i,j) =(4.0/9.0 )* (rho_                                                           );
	     (*fn)(i,j) =(1.0/9.0 )* (rho_ +3*         vy   +4.5*         vy  *         vy   - 1.5*u2);
	     (*fne)(i,j)=(1.0/36.0)* (rho_ +3*  (vx+   vy)  +4.5*(  vx +  vy) *(  vx+   vy)  - 1.5*u2);
	     (*fe)(i,j) =(1.0/9.0 )* (rho_ +3*   vx         +4.5*   vx        *   vx         - 1.5*u2);
	     (*fse)(i,j)=(1.0/36.0)* (rho_ +3*  (vx+ (-vy)) +4.5*(  vx +(-vy))*(  vx +(-vy)) - 1.5*u2);
	     (*fs)(i,j) =(1.0/9.0 )* (rho_ +3*       (-vy)  +4.5*         vy  *         vy   - 1.5*u2);
	     (*fsw)(i,j)=(1.0/36.0)* (rho_ +3*((-vx)+(-vy)) +4.5*((-vx)+(-vy))*((-vx)+(-vy)) - 1.5*u2);
	     (*fw)(i,j) =(1.0/9.0 )* (rho_ +3* (-vx)        +4.5*   vx        *   vx         - 1.5*u2);
	     (*fnw)(i,j)=(1.0/36.0)* (rho_ +3*((-vx)+  vy)  +4.5*((-vx)+  vy) *((-vx)+  vy)  - 1.5*u2);
	   }
   }

   //Initialize  data structures needed
   void init_data_struct()
   {  
	  
	  rho=new Array(sx+2,sy+2);
	  (*rho).fill(1.0);
	  ux=new Array(sx+2,sy+2);
	  (*ux).fill(0.0);
	  uy=new Array(sx+2,sy+2);
	  (*uy).fill(0.0);
	  fc=new Array(sx+2,sy+2);
	  (*fc).fill(4.0/9.0);
	  fn=new Array(sx+2,sy+2);
	  (*fn).fill(1.0/9.0);
	  fne=new Array(sx+2,sy+2);
	  (*fne).fill(1.0/36.0);
	  fe=new Array(sx+2,sy+2);
	  (*fe).fill(1.0/9.0);
	  fse=new Array(sx+2,sy+2);
	  (*fse).fill(1.0/36.0);
	  fs=new Array(sx+2,sy+2);
	  (*fs).fill(1.0/9.0);
	  fsw=new Array(sx+2,sy+2);
	  (*fsw).fill(1.0/36.0);
	  fw=new Array(sx+2,sy+2);
	  (*fw).fill(1.0/9.0);
	  fnw=new Array(sx+2,sy+2);
	  (*fnw).fill(1.0/36.0);
	  
	  //boundary flag values
	  //(north) 
      for(int i=0;i<=sx+1;i++)
      {
	    (*flag)(i,sy+1)=1;
      }
      //(south) 
      for(int i=0;i<=sx+1;i++)
      {
	    (*flag)(i,0)=1;
      }
      //(east) 
      for(int j=1;j<=sy;j++)
      {
        (*flag)(0,j)=1;
      }
      //(west) 
      for(int j=1;j<=sy;j++)
      {
        (*flag)(sx+1,j)=1;
      }
   }
	
   inline void vtkPrint(const std::string &filename)
	{ std::vector < int > dim(2);
	  std::vector < float > length(2);
	  dim[0]=sx;
	  dim[1]=sy;
	  length[0]=1;
	  length[1]=1;
	
    	  VTKFileWriter vtkwrite(dim,length,filename);
	      std::vector< std::vector< double > > u(2);
          u[0].resize(sx*sy);
          u[1].resize(sx*sy);
	      std::vector<real> r(dim[0]*dim[1]);
		  std::vector<real> flags(dim[0]*dim[1]);
          for (int i=1;i<=dim[0];i++)
             for (int j=1;j<=dim[1];j++)
	        {   flags[(j-1)*dim[0]+i-1]=(*flag)(i,j);
                r[(j-1)*dim[0]+i-1]=(*rho)(i,j);
                u[0][(j-1)*dim[0]+i-1]=(*ux)(i,j);
		        u[1][(j-1)*dim[0]+i-1]=(*uy)(i,j);
                }
        vtkwrite.WriteScalar(flags,"flags");
		vtkwrite.WriteScalar(r,"density");
        vtkwrite.WriteVector(u,"velocity");  
   }

  private:
	int sx,sy,vtks,tstep;
	double omg;
	Array *flag,*rho, *ux,*uy,*fc,*fn,*fne,*fe,*fse,*fs,*fsw,*fw,*fnw;
	std::string vtkname,vtknam;
	
   // put your members here

};

#endif //SOLVER_HH
