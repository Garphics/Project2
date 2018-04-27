#include <cmath>
#include "CImg.h"
#include <limits>
#include <algorithm>
#include <iostream>

const double PI = 3.141592653589793;

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

struct Vector4 {
  double a,b,c,d;
  Vector4(){
    a=0;
    b=0;
    c=0;
    d=0;
  }
  Vector4(double aa, double bb, double cc,double dd){
      a=aa;
      b=bb;
      c=cc;
      d=dd;
  }
  //normalize vector so length is 1
  Vector4 normalize(){
    double num = sqrt(a*a+b*b+c*c+d*d);
    return Vector4(a/num, b/num,c/num, d/num);
  }

  //overriding opertors to make life easier
  Vector4 operator + (Vector4 v){
    return Vector4(a+v.a,b+v.b,c+v.c,d+v.d);
  }
  Vector4 operator - (Vector4 v){
    return Vector4(a-v.a,b-v.b,c-v.c,d-v.d);
  }
  Vector4 operator * (double num){
    return Vector4 (a*num,b*num,c*num,d*num);
  }
  Vector4 operator / (double num){
    return Vector4(a/num,b/num,c/num,d/num);
  }
};


Vector cross (Vector v, Vector b){
  return Vector(v.y*b.z - v.z*b.y, v.z*b.x - v.x*b.z, v.x*b.y - v.y*b.x);
}

double dot4(Vector4 v1, Vector4 v2){
  return (v1.a*v2.a+v1.b*v2.b+v1.c*v2.c,v1.d*v2.d);
}

struct matrix4 {
  Vector4 v1;
  Vector4 v2;
  Vector4 v3;
  Vector4 v4;

  matrix4(){
    v1 = Vector4(0,0,0,0);
    v2 = Vector4(0,0,0,0);
    v3 = Vector4(0,0,0,0);
    v4 = Vector4(0,0,0,0);
  }

  matrix4(Vector4 vec1,Vector4 vec2,Vector4 vec3,Vector4 vec4){
    v1 = vec1;
    v2 = vec2;
    v3 = vec3;
    v4 = vec4;
  }

  matrix4 operator * (matrix4 M){
    matrix4 finalM = matrix4();
    finalM.v1.a = v1.a*M.v1.a + v1.b*M.v2.a + v1.c*M.v3.a + v1.d*M.v4.a;
    finalM.v1.b = v1.a*M.v1.b + v1.b*M.v2.b + v1.c*M.v3.b + v1.d*M.v4.b;
    finalM.v1.c = v1.a*M.v1.c + v1.b*M.v2.c + v1.c*M.v3.c + v1.d*M.v4.c;
    finalM.v1.d = v1.a*M.v1.d + v1.b*M.v2.d + v1.c*M.v3.d + v1.d*M.v4.d;

    finalM.v2.a = v2.a*M.v1.a + v2.b*M.v2.a + v2.c*M.v3.a + v2.d*M.v4.a;
    finalM.v2.b = v2.a*M.v1.b + v2.b*M.v2.b + v2.c*M.v3.b + v2.d*M.v4.b;
    finalM.v2.c = v2.a*M.v1.c + v2.b*M.v2.c + v2.c*M.v3.c + v2.d*M.v4.c;
    finalM.v2.d = v2.a*M.v1.d + v2.b*M.v2.d + v2.c*M.v3.d + v2.d*M.v4.d;

    finalM.v3.a = v3.a*M.v1.a + v3.b*M.v2.a + v3.c*M.v3.a + v3.d*M.v4.a;
    finalM.v3.b = v3.a*M.v1.b + v3.b*M.v2.b + v3.c*M.v3.b + v3.d*M.v4.b;
    finalM.v3.c = v3.a*M.v1.c + v3.b*M.v2.c + v3.c*M.v3.c + v3.d*M.v4.c;
    finalM.v3.d = v3.a*M.v1.d + v3.b*M.v2.d + v3.c*M.v3.d + v3.d*M.v4.d;

    finalM.v4.a = v4.a*M.v1.a + v4.b*M.v2.a + v4.c*M.v3.a + v4.d*M.v4.a;
    finalM.v4.b = v4.a*M.v1.b + v4.b*M.v2.b + v4.c*M.v3.b + v4.d*M.v4.b;
    finalM.v4.c = v4.a*M.v1.c + v4.b*M.v2.c + v4.c*M.v3.c + v4.d*M.v4.c;
    finalM.v4.d = v4.a*M.v1.d + v4.b*M.v2.d + v4.c*M.v3.d + v4.d*M.v4.d;
    return finalM;
  }
};

matrix4 cameraToWorld(Vector e, Vector d, Vector up){
  Vector zaxis = e - d;
  zaxis = zaxis.normalize();
  Vector xaxis = cross(up,zaxis);
  xaxis = xaxis.normalize();
  Vector yaxis = cross(zaxis,xaxis);

  matrix4 orientation = matrix4();
  orientation.v1 = Vector4(xaxis.x,yaxis.x,zaxis.x,0);
  orientation.v2 = Vector4(xaxis.y,yaxis.y,zaxis.y,0);
  orientation.v3 = Vector4(xaxis.z,yaxis.z,zaxis.z,0);
  orientation.v4 = Vector4(0,0,0,1);

  matrix4 translation = matrix4();
  translation.v1 = Vector4(1,0,0,0);
  translation.v2 = Vector4(0,1,0,0);
  translation.v3 = Vector4(0,0,1,0);
  translation.v4 = Vector4(-e.x,-e.y,-e.z,1);

  return (orientation * translation);


}

// struct Sphere{
//   Vector origin;
//   double radius;
//
//   Sphere(Vector o, double r){
//     origin =o;
//     radius = r;
//   }
//
//   //use quadratic formula to solve (d*d)t^2 + 2d*(e-origin)t + (e-c)*(e-c)-R^2 = 0
//     bool intersect(Vector e, Vector d, double &t){
//     double a = dot(d,d);
//     double b = 2*dot(d,(e-origin));
//     double c = dot((e-origin),(e-origin)) - radius*radius;
//     double discriminant = (b*b) - (4*a*c);
//     if(discriminant < 0){
//       return false;
//     }
//     else{
//       double num = sqrt(discriminant);
//       double t1 = ((-b)+num)/(2*a);
//       double t2 = ((-b)-num)/(2*a);
//       t=(t1 < t2) ? t1 :t2;
//       return true;
//     }
//   }
//   Vector normal(Vector ray){
//     return (origin-ray)/radius;
//   }
//
// };
//
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
//
// struct Object
// {
//     Sphere sphere = Sphere(Vector(0,0,0),0);
//     Vector color;
//    // Vecotr emission_color;
//     double kk,ks,ka,flect,fract,ior;
//     int pow;
//     Object(const Sphere &s,const Vector &c, double k1, double k2,double k3, int p, double r, double fl,double indexOfReflection)
//     {
//       sphere = s;
//       color = c;
//     //  emission_color = ec;
//       kk = k1;
//       ks = k2;
//       ka = k3;
//       flect = fl;
//       fract = r;
//       pow = p;
//       ior = indexOfReflection;
//     }
// };

struct Object{
  double kk,ks,ka,flect,fract,ior;
  int pow;
  Vector color;
  virtual bool intersect(Vector e, Vector d, double &t){};
  virtual Vector normal(Vector ray){};
};

struct Sphere : Object {

  double radius;
  Vector origin;

  Sphere(double k1,double k2,double k3, double fl, double r,double iOFr,int p, double rad, const Vector &c, Vector o){
    kk = k1;
    ks = k2;
    ka = k3;
    flect = fl;
    fract = r;
    ior = iOFr;
    pow = p;
    color = c;
    origin = o;
    radius = rad;
  }

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

double isShadow(Object** objs,int size,Vector viewRay, Vector light,Vector normal, int currentObj){
  Vector l = Vector(-light.x, -light.y, -light.z);
  for(int i = 0;i<size;i++){
    double distance = std::numeric_limits<double>::infinity();
    bool test = objs[i]->intersect(viewRay, l , distance);
    //to make sure shadow is rendered on correct side of object
    if(i == currentObj){
      // if((normal.x*l.x > 0) || (normal.y*l.y > 0) || (normal.x*l.x > 0)){
      //   break;
      // }
    }
    //check if it a sphere or a plane
    //else if(objs[i].isSphere){
      //remove the else when we add the above if statement about isSphere
      else if (objs[i]->intersect(viewRay, l, distance) && distance > 0 && objs[i]->fract > 0 ) {
        return objs[i]->fract;
      }
      else if (objs[i]->intersect(viewRay, l, distance) && distance > 0) {
        return 1;
      }
    //}
    //else if(objs.[i].isPlane){
      //check for intersection and return true if intersected
    //}
  }
  return -1;
}

Vector rayTrace(Vector e, Vector d, Object** objs, int size, int y, int step )
{
  Vector returnColor(0,0,0);
  int indexTemp;
  Object* hitObj;
  double intensity = 1.5;
  double t = std::numeric_limits<double>::infinity();
  double minT = t;

  Vector red(100,10,10);

  plane p(Vector(0,0,-17),Vector(0,1,0),red,.5,1.5,1,50);

  Vector light = Vector(2, -1, -2);

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
        double shadowVal = isShadow(objs,size,viewRay,plight,normal,indexTemp);
        if( shadowVal >= 0){
          la = shadowVal*p.ka*intensity ;
          //la = (la * y/192);
          //lk =0;
          //ls =0;
        }
        else {
          // lk = p.kk *(intensity) * std::max(0.0, dot(normal.normalize(), slight.normalize()));
          // ls = p.ks * (intensity) * pow(std::max(0.0, dot(normal.normalize(), r.normalize())), p.pow);
          la = p.ka * (intensity);
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
    //cout << "size: " << size << endl;
   //double t = std::numeric_limits<double>::infinity();
   // double t1 = std::numeric_limits<double>::infinity();
        if (objs[i]->intersect(e, d, t) && t >= 0)
        {
          //cout << "intersect"<<endl;
          if (t < minT){
            //cout << "t"<< endl;
            minT = t;
            hitObj = objs[i];
            indexTemp =i;
            Vector viewRay = e + d * t;
            Vector slight = Vector(light.x, light.y, light.z);
            Vector normal = hitObj->normal(viewRay);
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
            double shadowVal = isShadow(objs,size,viewRay,slight,normal,indexTemp);
            if(shadowVal > 0){
              //cout << "shadow"<<endl;
              la = (shadowVal)*objs[indexTemp]->ka;
              lk =0;
              ls =0;
            }
            else{
              lk = hitObj->kk *(intensity) * std::max(0.0, dot(normal.normalize(), slight.normalize()));
              ls = hitObj->ks * (intensity) * pow(std::max(0.0, dot(normal.normalize(), r.normalize())), hitObj->pow);
              la = hitObj->ka * (intensity);
            }
            Vector refractNormal = normal;
            Vector frenselNormal = normal;
            double kr;
            double fcosi = dot(frenselNormal,d);
            if(fcosi > 1){
              fcosi = 1;
            }
            if (fcosi < -1){
              fcosi = -1;
            }
            double fetai = 1;
            double fetat = 1.5;

            if(fcosi > 0){
              double ftemp = fetai;
              fetai = fetat;
              fetat = ftemp;
            }
            double m = 1 - fcosi*fcosi;
            if(m < 0){
              m=0;
            }
            double sint = (1/(fetai/fetat)) * sqrt(m);

            if(sint>=1){
              kr =1;
            }
            else{
              double n = 1-sint*sint;
              if(n < 0){
                n =0;
              }
              double cost = sqrt(n);
              fcosi = abs(fcosi);
              double Rs = ((fetat * fcosi) - (fetat*cost))/((fetat * fcosi) + (fetai*cost));
              double Rp = ((fetai * fcosi) - (fetat*cost))/((fetai * fcosi) + (fetat*cost));
              kr = (Rs * Rs + Rp * Rp)/2;
            }
            if(step < 10){
              //cout << "refract"<<endl;
              reflColor = rayTrace(viewRay, d - normal.normalize() * 2 * dot(d,normal.normalize()),objs,size,y,step+1);
              if(hitObj->fract > 0 && kr < 1){
                returnColor = Vector(
                  (hitObj->color.x * la + 255 * ls + hitObj->color.x * lk) + (kr)*reflColor.x,
                  (hitObj->color.y * la + 255 * ls + hitObj->color.y * lk) + (kr)*reflColor.y,
                  (hitObj->color.z * la + 255 * ls + hitObj->color.z * lk) + (kr)*reflColor.z);
              }
              else{
                returnColor = Vector(
                  (hitObj->color.x * la + 255 * ls + hitObj->color.x * lk) + hitObj->flect*reflColor.x,
                  (hitObj->color.y * la + 255 * ls + hitObj->color.y * lk) + hitObj->flect*reflColor.y,
                  (hitObj->color.z * la + 255 * ls + hitObj->color.z * lk) + hitObj->flect*reflColor.z);
              }
              if(hitObj->fract > 0 && kr < 1){
                Vector I = d;
                double cosi = dot(I,refractNormal);
                if(cosi > 1){
                  cosi =1;
                }
                if(cosi < -1){
                  cosi = -1;
                }
                double etai = 1;
                double etat = 1.5;
                if(cosi < 0){
                  cosi = -cosi;
                }
                else{
                  double temp = etai;
                  etai = etat;
                  etat = temp;
                  refractNormal = Vector(0,0,0) - refractNormal;
                }
                double eta = 1/(etai/etat);
                double k = 1 - eta*eta *(1-cosi*cosi);
                Vector refdir;
                if(k<0){
                  refdir = Vector(0,0,0);
                }
                else{
                  refdir = I * eta + refractNormal * (eta * cosi - sqrt(k));
                }
                //Vector normRefdir = refdir.normalize();
                refractColor = rayTrace(viewRay - refractNormal*1e-4, refdir.normalize(),objs,size,y,step+1);
                returnColor = Vector(returnColor.x + (1-kr)*hitObj->fract*refractColor.x, returnColor.y + (1-kr)*hitObj->fract*refractColor.y, returnColor.z + (1-kr)*hitObj->fract*refractColor.z );
              }
            }
            else{
              //cout << "color"<<endl;
              returnColor = Vector(
                (hitObj->color.x * la + 255 * ls + hitObj->color.x * lk),
                (hitObj->color.y * la + 255 * ls + hitObj->color.y * lk),
                (hitObj->color.z * la + 255 * ls + hitObj->color.z * lk));
                // return rgb;
            }
          }
         }
   }
   //cout << returnColor.x <<endl;
   return returnColor;
  }

int main()
{

  int imageWidth = 1000;
  int imageHeight = 1000;
  double intensity = 1;
  Vector black(0,0,0);
  Vector green(15,75,15);

  Object * test1;
  Sphere s1 = Sphere(1.5,1.5,.5,0,1,0,100,1,black,Vector(0,1,-6));
  test1 = &s1;

  Object * test2;
  Sphere s2 = Sphere(1.5,1.5,.5,.5,0,0,100,1,green,Vector(-.3,1,-12));
  test2 = &s2;

  Object * objs[2]={ test1, test2 };
  //Compute u,v,w basis vectors
  //Creating blank 256x256 image
  CImg<unsigned char> img(imageWidth,imageHeight,1,3,0);
  //for each pixel
  for (int y = 0; y<imageHeight; y++)
  {
    for(int x=0;x<imageWidth;x++)
    {
        //pixelNDC is the normalized pixel position
        //NDC is Normalized Device Coordinates
        //shifted 0.5 pixels because we want to be in the middle of the pixels
        //range of 0-1
        double PixelNDCx= (x + 0.5)/imageWidth;
        double PixelNDCy=( y + 0.5)/imageHeight;

        //maps the coordinates so that the center of the screen is origin
        double PixelScreenx = 2*PixelNDCx -1;
        double PixelScreeny = 1-2*PixelNDCy;

        //adjusting based on aspect ratio and also field of view angle
        double FOVangle = 45;
        FOVangle = (FOVangle*PI)/180;
        float ImageAspectRatio = imageWidth/(float)imageHeight;

        //We are now in Camera Space :)
        double PixelCamerax = (PixelScreenx)*ImageAspectRatio*tan(FOVangle/2);
        double PixelCameray = (PixelScreeny)*tan(FOVangle/2);

        //final coordinate of the pixel on the image plane is
        //(PixelCamerax, PixelCameray);
        //ray origin
        Vector o = Vector(0, 2, 0);
        double angleY = .2;
        //ray direction
        Vector dir = Vector(PixelCamerax, PixelCameray - angleY, -1) - o;
        dir = dir.normalize();

        matrix4 cameraToWorldM = cameraToWorld(o,dir,Vector(0,1,0));

        Vector4 dir4 = Vector4(PixelCamerax,PixelCameray- angleY,-1,0);
        matrix4 dirM = matrix4(dir4,Vector4(),Vector4(),Vector4());

        matrix4 newDirM = cameraToWorldM*dirM;

        dir = Vector(newDirM.v1.a,newDirM.v1.b, newDirM.v1.c);
        dir = dir.normalize();

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
