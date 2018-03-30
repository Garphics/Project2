#include <cmath>
#include "CImg.h"
#include <limits>
#include <algorithm>
using namespace cimg_library;

struct Vector {
  double x,y,z;
  Vector(){
    x=0;
    y-0;
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
  bool intersect(Vector e, Vector d, double &t){
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
      t=(t1 < t2) ? t1 :t2;
      return true;
    }
  }
  Vector normal(Vector ray){
    return (origin-ray)/radius;
  }

};

int main(){

  double intensity = .5;
  Sphere sphere1(Vector(150,150,50),50);
  Sphere sphere2(Vector(35,175,50),25);
  //Compute u,v,w basis vectors
  //Creating blank 256x256 image
  CImg<unsigned char> img(256,256,1,3,0);

  //for each pixel
  for (int y = 0; y<256; y++){
    for(int x=0;x<256;x++){
        //---compute viewing ray---//
        //distance from eye to point on half-line set purposely to a very high number
        double t = std::numeric_limits<double>::infinity();
        //position of eye
        Vector d = Vector(0,0,1);
        //position of pixel
        Vector e = Vector(x,y,0);
        //find first object hit by ray and its surface normal n
        //Vector d = s-e;
        //ray e+t*d
        Vector viewRay;
        if(sphere1.intersect(e,d,t) && t>0){
          //set pixel color
          //light source top middle
          viewRay = e + d*t;
          Vector light = Vector(0,1,1);
          Vector normal = sphere1.normal(viewRay);
          Vector r = normal*(2*dot(normal,light)) - light;
          // double lightDis = sqrt((light.x-sphere1.origin.x)*(light.x-sphere1.origin.x) + (light.y-sphere1.origin.y)*(light.y-sphere1.origin.y) + (light.z-sphere1.origin.z)*(light.z-sphere1.origin.z));
          double lk= (intensity)*std::max(0.0,dot(normal.normalize(),light.normalize()));
          double ls= 1.5*(intensity)*pow(std::max(0.0,dot(normal.normalize(),r.normalize())),50);
          double la = (intensity);
          if(0*la + 255*ls + 0*lk > 255){
            img(x,y,0) = 255;
          }else{
            img(x,y,0) = 0*la + 255*ls + 0*lk;
          }
          if(87*la + 255*ls + 87*lk > 255){
            img(x,y,1) = 255;
          }else{
            img(x,y,1) = 87*la + 255*ls + 87*lk;
          }
          if(255*la + 255*ls + 255*lk > 255){
            img(x,y,2) = 255;
          }else{
            img(x,y,2) = 255*la + 255*ls + 255*lk;
          }
        }
        else if(sphere2.intersect(e,d,t) && t>0){
          //set pixel color
          //light source top middle
          viewRay = e + d*t;
          Vector light = Vector(1,3,1);
          Vector normal = sphere2.normal(viewRay);
          Vector r = normal*(2*dot(normal,light)) - light;
          // double lightDis = sqrt((light.x-sphere2.origin.x)*(light.x-sphere2.origin.x) + (light.y-sphere2.origin.y)*(light.y-sphere2.origin.y) + (light.z-sphere2.origin.z)*(light.z-sphere2.origin.z));
          double lk= (intensity)*std::max(0.0,dot(normal.normalize(),light.normalize()));
          double ls= 1.5*(intensity)*std::max(0.0,dot(normal.normalize(),r.normalize()));
          double la = (intensity);
          if(60*la + 255*ls + 60*lk > 255){
            img(x,y,0) = 255;
          }else{
            img(x,y,0) = 60*la + 255*ls + 60*lk;
          }
          if(167*la + 255*ls + 167*lk > 255){
            img(x,y,1) = 255;
          }else{
            img(x,y,1) = 167*la + 255*ls + 167*lk;
          }
          if(28*la + 255*ls + 28*lk > 255){
            img(x,y,2) = 255;
          }else{
            img(x,y,2) = 28*la + 255*ls + 28*lk;
          }
        }
    }
  }
  img.display();
  return 0;
}
