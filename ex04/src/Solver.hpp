#ifndef SOLVER_HH
#define SOLVER_HH

#include "FileReader.hpp"
#include "VTKFileWriter.hpp"
#include "Particle.hpp"
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>




class Solver
{
  public:

   Solver(FileReader& parameters)
   {
	  vis_step=parameters.GetIntParameter("vis_space");
	  t_start=parameters.GetDoubleParameter("t_start");
	  t_end=parameters.GetDoubleParameter("t_end");
	  delta_t=parameters.GetDoubleParameter("delta_t");
	  xMin=parameters.GetDoubleParameter("x_min");
	  yMin=parameters.GetDoubleParameter("y_min");
	  zMin=parameters.GetDoubleParameter("z_min");
	  xMax=parameters.GetDoubleParameter("x_max");
	  yMax=parameters.GetDoubleParameter("y_max");
	  zMax=parameters.GetDoubleParameter("z_max");
	  r_cut=parameters.GetDoubleParameter("r_cut");
	  eps=parameters.GetDoubleParameter("epsilon");
	  sig=parameters.GetDoubleParameter("sigma");
	  name=parameters.GetStringParameter("name");

	  lineLengthx = floor ( (xMax-xMin)/r_cut );
	  lineLengthy = floor ( (yMax-yMin)/r_cut );
	  lineLengthz = floor ( (zMax-zMin)/r_cut );
	  
	  cellSizex = ( (xMax-xMin)/lineLengthx );
	  cellSizey = ( (yMax-yMin)/lineLengthy );
	  cellSizez = ( (zMax-zMin)/lineLengthz );
	  
	  particle.resize(lineLengthx);
	  forceOldx.resize(lineLengthx);
	  forceOldy.resize(lineLengthx);
	  forceOldz.resize(lineLengthx);
	  for(unsigned int i=0; i<lineLengthx; ++i)
	  {
	        forceOldx[i].resize(lineLengthy);
	        forceOldy[i].resize(lineLengthy);
	        forceOldz[i].resize(lineLengthy);
			particle[i].resize(lineLengthy);
	  }
	     
	  for(unsigned int i=0; i<lineLengthx; ++i)
	     for(unsigned int j=0; j<lineLengthy; ++j)
	     {
	        forceOldx[i][j].resize(lineLengthz);
	        forceOldy[i][j].resize(lineLengthz);
	        forceOldz[i][j].resize(lineLengthz);
	     }

	  partNum=parameters.pvec.size();
	  copyInfo(parameters.pvec);
	  
   }

   void copyInfo(std::vector<Particle>& parVec)
   {
	  int x,y,z;
	  for(unsigned int i = 0; i < partNum; ++i)
	  {
		x = (parVec[i].x-xMin)/cellSizex;	
		y = (parVec[i].y-yMin)/cellSizey;
		z = (parVec[i].z-zMin)/cellSizez;
		particle[x][y][z].push_back(parVec[i]);
	  }

   }

   void updatePartPos(double dlt)
   {
     
	double tempPos;
        Particle tempPart;
     	for(unsigned int i=0; i<lineLengthx; ++i) // Go on cells
	   for(unsigned int j=0; j<lineLengthy; ++j)
		for(unsigned int k=0; k<lineLengthz; ++k)	  
		  
			for(std::list<Particle>::iterator it = particle[i][j][k].begin(); it != particle[i][j][k].end(); ) // GO on particles in each cell
			{
				if( (*it).dlt != dlt ) // If position has not been updated yet in this time step
				{       
			   		tempPos = (*it).x + delta_t * (*it).vx + 0.5 * (*it).fx * delta_t * delta_t / (*it).mass;
					if( tempPos>xMax ) // Periodic boundary
					{	(*it).x = xMin + tempPos-xMax; }
					else if( tempPos < xMin )
						{(*it).x = xMax - (xMin- tempPos);}
					else
						{(*it).x = tempPos;		 }

			   		tempPos = (*it).y + delta_t * (*it).vy + 0.5 * (*it).fy * delta_t * delta_t / (*it).mass;
					if( tempPos>yMax )
						(*it).y = yMin + tempPos-yMax; 
					else if( tempPos < yMin )
						(*it).y = yMax - (yMin- tempPos);
					else
						(*it).y = tempPos;

			   		tempPos = (*it).z + delta_t * (*it).vz + 0.5 * (*it).fz * delta_t * delta_t / (*it).mass;
					if( tempPos>zMax )
						(*it).z = zMin + tempPos-zMax; 
					else if( tempPos < zMin )
						(*it).z = zMax - (zMin- tempPos);
					else
						(*it).z = tempPos;	
			
					tempPart = *it;  // Update paricle's cell
					it =  particle[i][j][k].erase(it);
					tempPart.dlt = dlt;
					particle[(tempPart.x-xMin)/cellSizex][(tempPart.y-yMin)/cellSizey][(tempPart.z-zMin)/cellSizez].push_back(tempPart);
				}
				else
					++it;
			
			}

   }

   void calcForce()
   {
      
     	for(unsigned int i=0; i<lineLengthx; ++i) // Go on cells
	   for(unsigned int j=0; j<lineLengthy; ++j)
		for(unsigned int k=0; k<lineLengthz; ++k)
		  
			for(std::list<Particle>::iterator it = particle[i][j][k].begin(); it != particle[i][j][k].end(); ++it) // Go on particles
			{
			        std::vector < std::vector < unsigned int > > pcell(27, std::vector<unsigned int>(3));
				// pcell will hold the coord for the cells surrounding the particle *it
				// pcell takes care of the periodic boundary too
				pcell[0][0] = i==0 ? lineLengthx-1:i-1;  pcell[0][1] = j==0 ? lineLengthy-1:j-1;  pcell[0][2] = k==0 ? lineLengthz-1:k-1;
				pcell[1][0] = i==0 ? lineLengthx-1:i-1;  pcell[1][1] = j==0 ? lineLengthy-1:j-1;  pcell[1][2] = k;
				pcell[2][0] = i==0 ? lineLengthx-1:i-1;  pcell[2][1] = j==0 ? lineLengthy-1:j-1;  pcell[2][2] = k+1<lineLengthz ? k+1:0;
				pcell[3][0] = i==0 ? lineLengthx-1:i-1;  pcell[3][1] = j;                         pcell[3][2] = k==0 ? lineLengthz-1:k-1;
				pcell[4][0] = i==0 ? lineLengthx-1:i-1;  pcell[4][1] = j;                         pcell[4][2] = k;
				pcell[5][0] = i==0 ? lineLengthx-1:i-1;  pcell[5][1] = j;                         pcell[5][2] = k+1<lineLengthz ? k+1:0;
				pcell[6][0] = i==0 ? lineLengthx-1:i-1;  pcell[6][1] = j+1<lineLengthy ? j+1:0;   pcell[6][2] = k==0 ? lineLengthz-1:k-1;
				pcell[7][0] = i==0 ? lineLengthx-1:i-1;  pcell[7][1] = j+1<lineLengthy ? j+1:0;   pcell[7][2] = k;
				pcell[8][0] = i==0 ? lineLengthx-1:i-1;  pcell[8][1] = j+1<lineLengthy ? j+1:0;   pcell[8][2] = k+1<lineLengthz ? k+1:0;
				pcell[9][0] = i==0 ? lineLengthx-1:i-1;  pcell[9][1] = j==0 ? lineLengthy-1:j-1;  pcell[9][2] = k==0 ? lineLengthz-1:k-1;
				pcell[10][0] = i;                        pcell[10][1] = j==0 ? lineLengthy-1:j-1; pcell[10][2] = k;
				pcell[11][0] = i;                        pcell[11][1] = j==0 ? lineLengthy-1:j-1; pcell[11][2] = k+1<lineLengthz ? k+1:0;
				pcell[12][0] = i;                        pcell[12][1] = j;                        pcell[12][2] = k==0 ? lineLengthz-1:k-1;
				pcell[13][0] = i;                        pcell[13][1] = j;                        pcell[13][2] = k;
				pcell[14][0] = i;                        pcell[14][1] = j;                        pcell[14][2] = k+1<lineLengthz ? k+1:0;
				pcell[15][0] = i;                        pcell[15][1] = j+1<lineLengthy ? j+1:0;  pcell[15][2] = k==0 ? lineLengthz-1:k-1;
				pcell[16][0] = i;                        pcell[16][1] = j+1<lineLengthy ? j+1:0;  pcell[16][2] = k;
				pcell[17][0] = i;                        pcell[17][1] = j+1<lineLengthy ? j+1:0;  pcell[17][2] = k+1<lineLengthz ? k+1:0;
				pcell[18][0] = i+1<lineLengthx ? i+1:0;  pcell[18][1] = j==0 ? lineLengthy-1:j-1; pcell[18][2] = k==0 ? lineLengthz-1:k-1;
				pcell[19][0] = i+1<lineLengthx ? i+1:0;  pcell[19][1] = j==0 ? lineLengthy-1:j-1; pcell[19][2] = k;
				pcell[20][0] = i+1<lineLengthx ? i+1:0;  pcell[20][1] = j==0 ? lineLengthy-1:j-1; pcell[20][2] = k+1<lineLengthz ? k+1:0;
				pcell[21][0] = i+1<lineLengthx ? i+1:0;  pcell[21][1] = j;                        pcell[21][2] = k==0 ? lineLengthz-1:k-1;
				pcell[22][0] = i+1<lineLengthx ? i+1:0;  pcell[22][1] = j;                        pcell[22][2] = k;
				pcell[23][0] = i+1<lineLengthx ? i+1:0;  pcell[23][1] = j;                        pcell[23][2] = k+1<lineLengthz ? k+1:0;
				pcell[24][0] = i+1<lineLengthx ? i+1:0;  pcell[24][1] = j+1<lineLengthy ? j+1:0;  pcell[24][2] = k==0 ? lineLengthz-1:k-1;
				pcell[25][0] = i+1<lineLengthx ? i+1:0;  pcell[25][1] = j+1<lineLengthy ? j+1:0;  pcell[25][2] = k;
				pcell[26][0] = i+1<lineLengthx ? i+1:0;  pcell[26][1] = j+1<lineLengthy ? j+1:0;  pcell[26][2] = k+1<lineLengthz ? k+1:0;

				double tempfx=0, tempfy=0, tempfz=0;
				for(unsigned int ii=0; ii<27; ii++) // Go on cells on +ve direction
					for(std::list<Particle>::iterator itq = particle[pcell[ii][0]][pcell[ii][1]][pcell[ii][2]].begin(); itq != particle[pcell[ii][0]][pcell[ii][1]][pcell[ii][2]].end(); ++itq)
					{	
						if( ((*itq).x!=(*it).x) || ((*itq).y!=(*it).y) || ((*itq).z!=(*it).z)  ) // Force was not updated yet and not same particle
						{      
							double r = sqrt( ( (*itq).x-(*it).x ) * ( (*itq).x-(*it).x ) + ( (*itq).y-(*it).y ) * ( (*itq).y-(*it).y ) + 													( (*itq).z-(*it).z ) * ( (*itq).z-(*it).z ) );
							if((r_cut-r) > 0) //within r_cut
							{
								double rSqr = r*r; 
								double val = pow(sig/r,6.0);
								double tempF = ( (1.0/rSqr) * val * ( 1.0 - 2.0 * val ) ); 
								tempfx += 24.0*eps*tempF*((*itq).x-(*it).x );
                                tempfy += 24.0*eps*tempF*((*itq).y-(*it).y ); 
                                tempfz += 24.0*eps*tempF*((*itq).z-(*it).z );
							}
						}
								
					}
				pcell.clear();
				(*it).fx=tempfx;
				(*it).fy=tempfy;
				(*it).fz=tempfz;				
			}	
   }

   void updateVel()
   {
	std::list<double>::iterator itfx, itfy, itfz;
     	for(unsigned int i=0; i<lineLengthx; ++i)
	   for(unsigned int j=0; j<lineLengthy; ++j)
		for(unsigned int k=0; k<lineLengthz; ++k)
		{	   
			itfx = forceOldx[i][j][k].begin(); 					
			itfy = forceOldy[i][j][k].begin(); 					
			itfz = forceOldz[i][j][k].begin(); 			
	
			for(std::list<Particle>::iterator it = particle[i][j][k].begin(); it != particle[i][j][k].end(); ++it)
			{
				
				(*it).vx = (*it).vx + 0.5 * ( (*it).fx + (*itfx) ) * delta_t / (*it).mass; ++itfx; 					
				(*it).vy = (*it).vy + 0.5 * ( (*it).fy + (*itfy) ) * delta_t / (*it).mass; ++itfy;					
				(*it).vz = (*it).vz + 0.5 * ( (*it).fz + (*itfz) ) * delta_t / (*it).mass; ++itfz;					
			}

		}
   }


   void copyOldForces()
   {
	
	for(unsigned int i=0; i<lineLengthx; ++i)
	    for(unsigned int j=0; j<lineLengthy; ++j)
		for(unsigned int k=0; k<lineLengthz; ++k)
		{       
		        forceOldx[i][j][k].clear();
			forceOldy[i][j][k].clear();
			forceOldz[i][j][k].clear();
			
			for(std::list<Particle>::iterator it = particle[i][j][k].begin(); it != particle[i][j][k].end(); ++it)
			{
				forceOldx[i][j][k].push_back( (*it).fx ); 					
				forceOldy[i][j][k].push_back( (*it).fy ); 					
				forceOldz[i][j][k].push_back( (*it).fz ); 					

			}
		}
	
   }

   // Solves the L-B eq
   void Solve()
   {
     
     double t=t_start;
     int n=0;
     std::stringstream nth;
     while(t<=t_end)
     {
         
	  updatePartPos(t); //std::cout<<t<<std::endl;
	  copyOldForces();
	  calcForce();
	  updateVel();
       
          if( (vis_step!=0) && (n%vis_step==0) )
          {
             nth.str("");
             nth<<n;
             vtkPrint("./ResultsVTK/"+std::string(name)+nth.str());
          }
          
          t+=delta_t;
          n++;
     }
     WriteResults(); //writes all the particles to a text file in ResultsTXT

   }
	
   void WriteResults() //print 
   {
     std::ofstream myfile;
     myfile.open("ResultsTXT/particle.txt");
     if (myfile.is_open())
     {
        for(unsigned int i=0; i<lineLengthx; ++i)
	   for(unsigned int j=0; j<lineLengthy; ++j)
	      for(unsigned int k=0; k<lineLengthz; ++k)	   
		  for(std::list<Particle>::iterator it = particle[i][j][k].begin(); it != particle[i][j][k].end(); ++it)
		  {
         	  	 myfile<<(*it).fx<<" "<<(*it).fy<<" "<<(*it).fz<<" "<<(*it).vx<<" "<<(*it).vy<<" "<<(*it).vz<<" "<<(*it).x<<" "<<(*it).y<<" "<<(*it).z<<" "<<(*it).mass<<"\n";
       		  }
     
       
    }
   }
   

  inline void vtkPrint(const std::string &filename)
   {
     
     std::vector< std::vector< double > > force(3);
     force[0].resize(partNum);
     force[1].resize(partNum);
     force[2].resize(partNum);

     std::vector< std::vector< double > > vel(3);
     vel[0].resize(partNum);
     vel[1].resize(partNum);
     vel[2].resize(partNum);
     
     std::vector< std::vector< double > > pos(3);
     pos[0].resize(partNum);
     pos[1].resize(partNum);
     pos[2].resize(partNum);
     
     std::vector<double> pmass(partNum);
     
     VTKFileWriter vtkwrite(filename);

     int t=0;
     for(unsigned int i=0; i<lineLengthx; ++i)
	for(unsigned int j=0; j<lineLengthy; ++j)
	   for(unsigned int k=0; k<lineLengthz; ++k)	   
		for(std::list<Particle>::iterator it = particle[i][j][k].begin(); it != particle[i][j][k].end(); ++it)
		{
			
         		 force[0][t]=(*it).fx;
			 force[1][t]=(*it).fy;
			 force[2][t]=(*it).fz;
			 vel[0][t]=(*it).vx;
			 vel[1][t]=(*it).vy;
			 vel[2][t]=(*it).vz;
			 pos[0][t]=(*it).x;
			 pos[1][t]=(*it).y;
			 pos[2][t]=(*it).z;
			 pmass[t]=(*it).mass;
			++t;

       		}
	   
     std::stringstream str_nr;
     str_nr.str("");
     str_nr<<partNum;
     vtkwrite.WriteVectorP(pos, str_nr.str());
     vtkwrite.WriteScalar(pmass,"mass");
     vtkwrite.WriteVector(force,"force");
     vtkwrite.WriteVector(vel,"velocity");
     
   }
   


  private:
	std::string name;
	double  cellSizex, cellSizey, cellSizez, t_start, t_end,  delta_t, xMax, yMax, zMax, xMin, yMin, zMin, vx, vy ,vz ,r_cut, sig, eps;
	int vis_step;
	unsigned int lineLengthx, lineLengthy, lineLengthz, partNum;
	std::vector< std::vector< std::vector< std::list<Particle> > > > particle;
	std::vector< std::vector< std::vector< std::list<double> > > > forceOldx, forceOldy, forceOldz;
   // put your members here

};

#endif //SOLVER_HH
