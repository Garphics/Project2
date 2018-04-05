#include <cmath>
#include "CImg.h"
#include <limits>
#include <algorithm>
#include <iostream>
using namespace cimg_library;

struct Vector {
  double x,y,z;
  Vector(){
    x=0;
    y=0;
    z=0;
  }
  Vector(double a, double b, double c){
      x=a;
      y=b;
      z=c;
  }
  //normalize vector so length is 1
  Vector normalize(){
    double num = sqrt(x*x+y*y+z*z);
    return Vector(x/num, y/num,z/num);
  }

  //overriding opertors to make life easier
  Vector operator + (Vector v){
    return Vector(x+v.x,y+v.y,z+v.z);
  }
  Vector operator - (Vector v){
    return Vector(x-v.x,y-v.y,z-v.z);
  }
  Vector operator * (double d){
    return Vector(x*d,y*d,z*d);
  }
  Vector operator / (double d){
    return Vector(x/d,y/d,z/d);
  }
};

double dot(Vector v, Vector b){
  return (v.x*b.x+v.y*b.y+v.z*b.z);
}

struct Sphere{
  Vector origin;
  double radius;

  Sphere(Vector o, double r){
    origin =o;
    radius = r;
  }

  //use quadratic formula to solve (d*d)t^2 + 2d*(e-origin)t + (e-c)*(e-c)-R^2 = 0
  bool intersect(Vector e, Vector d, double &t1, double &t2){
    double a = dot(d,d);
    double b = 2*dot(d,(e-origin));
    double c = dot((e-origin),(e-origin)) - radius*radius;
    double discriminant = (b*b) - (4*a*c);
    if(discriminant < 0){
      return false;
    }
    else{
      double num = sqrt(discriminant);
      double t1 = ((-b)+num)/(2*a);
      double t2 = ((-b)-num)/(2*a);
      t = (t1 < t2) ? t1 :t2;
      return true;
    }
  }
  Vector normal(Vector ray){
    return (origin-ray)/radius;
  }

};
struct Object
{
    Sphere sphere = Sphere(Vector(0,0,0),0);
    Vector color;
    float ref_factor;
    float trans_factor;
    double kk,ks,ka;
    int pow;

    Object(const Sphere &s,const Vector &c,double k1, double k2,double k3,int p)
    {
      sphere = s;
      color = c;
      kk = k1;
      ks = k2;
      ka = k3;
      pow = p;
    }
};
/*
Vector/color Tracer (Vecotor ray, depth){
	max depth = INFINITE;
	if(depth > maxdepth){
		return background color;
	}
	q = isIntersect
	if(q = null){
		return backgorund color;
	}
}

 
*/
Vector Raytrace(Vector origin, Vector dir, objs, depth)
{
  Vector background(0,0,0);
  if(depth < max)
  {
    return background;
  }
  float tn = 10000000; //distance to nearest object
  Object *hitObj = NULL; //pointer to the interscted object
  Vector hitColor = (0,0,0); // color of intersected point

  //find interest of this ray with the sphere in the scene
  for(int i = 0; i < objs.size(); i++) 
  {
    float t0 = INFINITY, ti = INFINITY;
    if(objs[i].intersects(origin,dir,t0,t1)) 
    {
      if(t0 < 0){t0 = t1;}
      if(t0 < tnear){ 
        tnear = 0;
        hitObj = &objs[i];}
    }

  }
  //if there is no interestion with an object
  if (!hitObj) {return background color;};
  Vector surfaceColor = 0,0,0; //color of the ray/surface of the object intersected by the ray
  Vector hitPoint = origin + dir * tnear;
  Vector hitNorm =  objs[i].sphere.normal(hitPoint);
  bool inside = false;  

  //normalize hitnorm
  // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
  if(dir.dot(nhit) > 0) 
   {
      hitPoint = -hitNorm; inside = true;
   }
  if(obj.reflective || obj.transparent) 
  {
    float facing = -raydir.dot(hipPoint);
    float frensel =  1 * 0.1 + pow(1-facing,3) * (1 - 0.1);


     // compute reflection direction
    Vector reflectdir = raydir - hitPoint * 2 * raydir.dot(hitPoint);
    reflectdir.normalize();
    Vector refelctRay = trace(hitPoint + hitNorm, reflectdir,objs,depth+1);
    Vecotr refractRay = 0,0,0;

    // if the sphere is also transparent compute refraction ray (transmission)
    if(obj.transparent)
    {
      float  ior = 1.1;
      float eta = (inside) ? ior : 1 / ior; //are we inside or outside the surface?
      float cosi = -hitPoint.dot(raydir);
      float k = 1 - (eta*eta) * (1-cosi*cosi);
      Vecotr refractdir = rayfir * eta + hitPoint * (era * cosi - sqrt(k));
      refractdir.normalize();
      refractRay = trace(hitPoint hitNorm, refractdir,spheres, depth + 1);
    }
      // the result is a mix of reflection and refraction (if the sphere is transparent)
     surfaceColor = (refeltRay * fresnel +refractRay * (1-frensel) * obj.transparency * obj.surfaceColor);
  }

  else {
    for(int i = 0; i < objs.size; i++){
      if(objs[i].emissionColor.x > 0){
        Vector transmit = 1;
        Vector lightDir = obj[i].center - hitPoint;
        lightDirc.normalize();
        for(int j = 0; j< objs.size(); ++j){
          if(i != j){
            float t0, t1;
            if(objs[j].intersect(hitPoint+hitNorm, lightDir,t0,t1)){
              transmit = 0;
              break;
            }
          }
        }
        surfaceColor += obj.sphere.surfaceColor * transmit * max(0,hitPoint.dot(lightDir))*objs[i].emissionColor;
      }
    }
  }
 return surfaceColor + obj.sphere.emissionColor
}
int main()
{

  double intensity = .5;
  Sphere sphere1(Vector(150,150,50),50);
  Sphere sphere2(Vector(35,175,50),25);
  Vector blue(0,87,255);
  Vector green(60,167,28);
  Object objs[2] = {Object(sphere1,blue,.5,1.5,1,50),Object(sphere2,green,1,1.5,1,50)} ;//you can add spheres as objects and I could expand it to other objects
  //Compute u,v,w basis vectors
  //Creating blank 256x256 image
  CImg<unsigned char> img(256,256,1,3,0);
  double intensity = .5;
  Sphere sphere1(Vector(150,150,50),50);
  Sphere sphere2(Vector(35,175,50),25);
  Vector blue(0,87,255);
  Vector green(60,167,28);
  Object objs[2] = {Object(sphere1,blue,.5,1.5,1,50),Object(sphere2,green,1,1.5,1,50)} ;
  //for each pixel



  for (int y = 0; y<256; y++)
  {
    for(int x=0;x<256;x++)
    {
        Vector o = Vector(0, 0, 1);
        //position of pixel
        Vector dir = Vector(x, y, 0);
        Vector color = Raytrace(o,dir,objs,0)
    
  }
  img.display();
  return 0;
}