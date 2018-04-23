#include <cmath>
#include "CImg.h"
#include <limits>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <random>

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

struct Box {

	Vector color;
	Vector normal;
	double kk, ks, ka;
	int pow;

	Box(const Vector &vMin, const Vector &vMax, Vector n, Vector c, double k1, double k2, double k3, int p) {

		bounds[0] = vMin;
		bounds[1] = vMax;

		color = c;
		normal = Vector(n.x, n.y, n.z);
		kk = k1;
		ks = k2;
		ka = k3;
		pow = p;

	}

	bool intersect(Vector e, Vector d, double &t) {

		float tMin;
		float tMax;
		float tYMin;
		float tYMax;
		float tZMin;
		float tZMax;

		tMin = (bounds[0].x - e.x) / d.x;
		tMax = (bounds[1].x - e.x) / d.x;

		if (tMin > tMax) {

			swap(tMin, tMax);
		}

		tYMin = (bounds[0].y - e.y) / d.y;
		tYMax = (bounds[1].y - e.y) / d.y;

		if(tYMin > tYMax) {

			swap(tYMin, tYMax);
		}

		if((tMin > tYMax) || (tYMin > tMax)) {

			return false;
		}

		if(tYMin > tMin) {

			tMin = tYMin;
		}

		if(tYMax < tMax) {

			tMax = tYMax;
		}

		tZMin = (bounds[2].z - e.z) / d.z;
		tZMax = (bounds[2].z - e.z) / d.z;

		if(tZMin > tZMax) {

			swap(tZMin, tZMax);
		}

		if((tMin > tZMax) || (tZMin > tMax)) {

			return false;
		}

		if(tZMin > tMin) {

			tMin = tZMin;
		}

		if(tZMax < tMax) {

			tMax = tZMax;
		}

		t = tMin;

		if(t < 0) {

			t = tMax;

			if(t < 0) {

				return false;
			}
		}

		return true;

	}

	Vector bounds[2];

};

struct plane
{
    Vector point;
    Vector normal;
    Vector color;
    double kk,ks,ka;
    int pow;
    plane(Vector po, Vector n,Vector c, double k1, double k2, double k3, int p)
    {
        point = po;
        normal = Vector(n.x,n.y,n.z);
        color = c;
        kk = k1;
        ks = k2;
        ka = k3;
        pow = p;

    }

    bool intersect(Vector RayOri,  Vector RayDir, double &t)
    {
        double ln = dot(normal,RayDir);
        if(abs(ln) > 1e-6)
        {
            Vector PROri = point - RayOri;
            t = dot(PROri,normal)/ln;
            return(t >= 0);
        }
        return false;
    }
};

struct Object
{
    Sphere sphere = Sphere(Vector(0,0,0),0);
    Vector color;
   // Vecotr emission_color;
    double kk,ks,ka,flect,fract;
    int pow;
    int refr;
    Object(const Sphere &s,const Vector &c, double k1, double k2,double k3, int p, int r, double fl)
    {
      sphere = s;
      color = c;
      refr = r;
    //  emission_color = ec;
      kk = k1;
      ks = k2;
      ka = k3;
      flect = fl;
    //  fr = fract;
      pow = p;
    }
};

bool isShadow(Object objs[],int size,Vector viewRay, Vector light,Vector normal, int currentObj){
  Vector l = Vector(-light.x, -light.y, -light.z);
  for(int i = 0;i<size;i++){
    double distance = std::numeric_limits<double>::infinity();
    bool test = objs[i].sphere.intersect(viewRay, l , distance);
    //to make sure shadow is rendered on correct side of object
    if(i == currentObj){
      // if((normal.x*l.x > 0) || (normal.y*l.y > 0) || (normal.x*l.x > 0)){
      //   break;
      // }
    }
    //check if it a sphere or a plane
    //else if(objs[i].isSphere){
      //remove the else when we add the above if statement about isSphere
      else if (objs[i].sphere.intersect(viewRay, l, distance) && distance > 0) {
        return true;
      }
    //}
    //else if(objs.[i].isPlane){
      //check for intersection and return true if intersected
    //}
  }
  return false;
}

Vector rayTrace(Vector e, Vector d, Object objs[], int size, int y, int step )
{
  // Vector rgb;
  // Vector none = Vector(0,0,0);
  Vector returnColor(0,0,0);
  int indexTemp;
  Object hitObj = Object(Sphere(Vector(0,0,0),0),Vector(0,0,0),0,0,0,0,0,0);
  double intensity = 1.5;
  double t = std::numeric_limits<double>::infinity();
  double minT = t;
   Vector hitColor(0,0,0);
 //  Vector rgb;
  // Sphere hitObj = NULL;

  Vector red(200,30,30);

  Vector blue(0, 0, 255);

  Box b(Vector(50,50,50), Vector(200, 200, 90), Vector(0,-2,-1), blue, 0.5, 1.5, 1, 50);

  plane p(Vector(0,100,300),Vector(0,-2,-1),red,.5,1.5,1,50);

  Vector light = Vector(2, 1, 2);


  //Box intersection
  if (b.intersect(e, d, t) && t > 0) {
      if (t < minT){
        minT = t;
        //set pixel color
        //light source top middle
        Vector viewRay = e + d * t;
        Vector plight = Vector(light.x, light.y, light.z);
        Vector normal = b.normal;
        indexTemp = -1;
        double la;
        Vector r = normal * (2 * dot(normal, plight)) - plight;
        if(isShadow(objs,size,viewRay,plight,normal,indexTemp)){
          la = .5*b.ka * (intensity);
          la = (la * y/192);
          //lk =0;
          //ls =0;
        }
        else {
          la = b.ka * (intensity);
          la = (la * y/192);
        }
        if (b.color.x * la > 255) {
          returnColor.x = 255;
        } else {
          returnColor.x = b.color.x * la;
        }
        if (b.color.y * la > 255) {
          returnColor.y = 255;
        } else {
          returnColor.y = b.color.y * la;
        }
        if (b.color.z * la > 255) {
          returnColor.z = 255;
        } else {
          returnColor.z = b.color.z * la;
        }
      }
  }

  if (p.intersect(e, d, t) && t > 0) {
      if (t < minT){
        minT = t;
        //set pixel color
        //light source top middle
        Vector viewRay = e + d * t;
        Vector plight = Vector(light.x, light.y, light.z);
        Vector normal = p.normal;
        indexTemp = -1;
        double la;
        Vector r = normal * (2 * dot(normal, plight)) - plight;
        if(isShadow(objs,size,viewRay,plight,normal,indexTemp)){
          la = .5*p.ka * (intensity);
          la = (la * y/192);
          //lk =0;
          //ls =0;
        }
        else {
          la = p.ka * (intensity);
          la = (la * y/192);
        }
        if (p.color.x * la > 255) {
          returnColor.x = 255;
        } else {
          returnColor.x = p.color.x * la;
        }
        if (p.color.y * la > 255) {
          returnColor.y = 255;
        } else {
          returnColor.y = p.color.y * la;
        }
        if (p.color.z * la > 255) {
          returnColor.z = 255;
        } else {
          returnColor.z = p.color.z * la;
        }
      }
  }

  for(int i = 0;i<size;i++)
  {
   // double t = std::numeric_limits<double>::infinity();
   // double t1 = std::numeric_limits<double>::infinity();
        if (objs[i].sphere.intersect(e, d, t) && t >= 0)
        {
          if (t < minT){
            minT = t;
            hitObj = objs[i];
            indexTemp =i;
            Vector viewRay = e + d * t;
            Vector slight = Vector(light.x, light.y, light.z);
            Vector normal = hitObj.sphere.normal(viewRay);
            Vector r = normal * (2 * dot(normal, slight)) - slight;
            Vector reflColor;
            double lk;
            double ls;
            double la;
            double nconst = 1.5;
            double c1 = dot(normal.normalize(),viewRay);
            double c2 = sqrt(1-nconst*nconst*(1-c1*c1));
            Vector T = viewRay * nconst + normal.normalize()*(nconst*c1 - c2);
            Vector refractColor =  Vector(0,0,0);
            //Vector refractColor = rayTrace(viewRay, T, objs,size,y,100);
            if(isShadow(objs,size,viewRay,slight,normal,indexTemp)){
              la = objs[indexTemp].ka * (intensity);
              lk =0;
              ls =0;
            }
            else{
              lk = hitObj.kk *(intensity) * std::max(0.0, dot(normal.normalize(), slight.normalize()));
              ls = hitObj.ks * (intensity) * pow(std::max(0.0, dot(normal.normalize(), r.normalize())), hitObj.pow);
              la = hitObj.ka * (intensity);
            }
            if(step < 5){
              reflColor = rayTrace(viewRay, d - normal.normalize() * 2 * dot(d,normal.normalize()),objs,size,y,step+1);
              returnColor = Vector(
                (hitObj.color.x * la + 255 * ls + hitObj.color.x * lk) + hitObj.flect*reflColor.x,
                (hitObj.color.y * la + 255 * ls + hitObj.color.y * lk) + hitObj.flect*reflColor.y,
                (hitObj.color.z * la + 255 * ls + hitObj.color.z * lk) + hitObj.flect*reflColor.z);
              if(hitObj.refr = 1){
                returnColor = Vector(returnColor.x + refractColor.x, returnColor.y + refractColor.y, returnColor.z + refractColor.z );
              }
            }
            else{
              returnColor = Vector(
                (hitObj.color.x * la + 255 * ls + hitObj.color.x * lk),
                (hitObj.color.y * la + 255 * ls + hitObj.color.y * lk),
                (hitObj.color.z * la + 255 * ls + hitObj.color.z * lk));
              if(hitObj.refr = 1){
                returnColor = Vector(returnColor.x + refractColor.x, returnColor.y + refractColor.y, returnColor.z + refractColor.z );
              }
                // return rgb;
            }
          }
        }
   }

   return returnColor;
  }

int main()
{

  double intensity = 1.5;
  Sphere sphere1(Vector(165,130,140),50);
  Sphere sphere2(Vector(85,170,100),30);
  Vector blue(15,15,75);
  Vector green(15,75,15);

  Object objs[2] = {Object(sphere1,blue,1.5,1.5,.5,50,0,.3),Object(sphere2,green,1.5,1.5,.5,50,1,.3)} ;//you can add spheres as objects and I could expand it to other objects
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
        Vector mycolor = rayTrace(o,dir,objs, sizeof(objs)/sizeof(objs[0]),y,0);
      //  std::cout << mycolor << std::endl;
        if (mycolor.x > 255)
          {
            img(x, y, 0) = 255;
          } else
          {
            img(x, y, 0) = mycolor.x;
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
