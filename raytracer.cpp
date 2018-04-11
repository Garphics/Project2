#include <cmath>
#include "CImg.h"
#include <limits>
#include <algorithm>
#include <iostream>
using namespace cimg_library;
using namespace std;
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
struct Object
{
    Sphere sphere = Sphere(Vector(0,0,0),0);
    Vector color;
   // Vecotr emission_color;
    double kk,ks,ka,flect,fract;
    int pow;

    Object(const Sphere &s,const Vector &c, double k1, double k2,double k3, int p)
    {
      sphere = s;
      color = c;
    //  emission_color = ec;
      kk = k1;
      ks = k2;
      ka = k3;
     // fl = flect;
    //  fr = fract;
      pow = p;
    }
};


Vector rayTrace(Vector e, Vector d, Object objs[], int size )
{
  // Vector rgb;
  // Vector none = Vector(0,0,0);
  int step = 0;
  Object hitObj = Object(Sphere(Vector(0,0,0),0),Vector(0,0,0),0,0,0,0);
  double intensity = 0.5;
  double t = std::numeric_limits<double>::infinity();
   Vector hitColor(0,0,0);
 //  Vector rgb; 
  // Sphere hitObj = NULL;
  for(int i = 0;i<size;i++) 
  {
   // double t = std::numeric_limits<double>::infinity();
   // double t1 = std::numeric_limits<double>::infinity();
        if (objs[i].sphere.intersect(e, d, t) && t >= 0) 
        {
           hitObj = objs[i];

        }
            
   }
        Vector viewRay = e + d * t;
        Vector light = Vector(1, 0, 1);
        Vector normal = hitObj.sphere.normal(viewRay);
        Vector r = normal * (2 * dot(normal, light)) - light;  
        double lk = hitObj.kk *(intensity) * std::max(0.0, dot(normal.normalize(), light.normalize()));
        double ls = hitObj.ks * (intensity) * pow(std::max(0.0, dot(normal.normalize(), r.normalize())), hitObj.pow);
        double la = hitObj.ka * (intensity);        
        return Vector(
          (hitObj.color.x * la + 255 * ls + hitObj.color.x * lk),
          (hitObj.color.y * la + 255 * ls + hitObj.color.y * lk),
          (hitObj.color.z * la + 255 * ls + hitObj.color.z * lk)); 
    // return rgb;
     
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
  //for each pixel


  for (int y = 0; y<256; y++)
  {
    for(int x=0;x<256;x++)
    {
        Vector dir = Vector(0, 0, 1);
        //position of pixel
        Vector o = Vector(x, y, 0);
        Vector mycolor = rayTrace(o,dir,objs, sizeof(objs)/sizeof(objs[0]));
      //  std::cout << mycolor << std::endl;
        if (mycolor.x > 255) 
          {
            img(x, y, 0) = 255;
          } else 
          {
          }

          if (mycolor.y > 255) 
          {
            img(x, y, 1) = 255;
          } else 
          {
            img(x, y, 1) = mycolor.y;
          }

          if (mycolor.z > 255) 
          {
            img(x, y, 2) = 255;
          } else 
          {
            img(x, y, 2) = mycolor.z;
          }
    
  } 
}
img.display();
  return 0;
}
