#ifndef PARTICLE_HH
#define PARTICLE_HH

#include "Types.hpp"
#include "FileReader.hpp"
#include <list>
#include <vector>


class Particle{
  
  
  
  
  public:
  double dlt, mass, x, y, z, vx, vy, vz, fx, fy, fz;
 
  Particle(): dlt(-1),mass(0.0),x(0.0),y(0.0),z(0.0),vx(0.0),vy(0.0),vz(0.0),fx(0.0),fy(0.0),fz(0.0){
    
  }
  
  Particle( double dlt, double mass, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz)
    : dlt(dlt),mass(mass),x(x),y(y),z(z),vx(vx),vy(vy),vz(vz),fx(fx),fy(fy),fz(fz){
      
    }
    
  Particle &operator=(const Particle& p) {
    
    dlt = p.dlt;
    mass=p.mass;
    x=p.x;
    y=p.y;
    z=p.z;
    vx=p.vx;
    vy=p.vy;
    vz=p.vz;
    fx=p.fx;
    fy=p.fy;
    fz=p.fz;

    return *this;
  } 

  private:
 
  
};

#endif //PARTICLE_HH
