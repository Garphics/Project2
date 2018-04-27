#include <iostream>
using namespace std;

struct Object{
  double kk,ks,ka,flect,fract,ior;
  int pow;
  // Vector color;
  // Vector position;
  virtual int itWorks(){};
};

struct Sphere : Object {
  double radius;

  Sphere(double k1,double k2,double k3, double fl, double r,double iOFr,int p, double rad){
    kk = k1;
    ks = k2;
    ka = k3;
    flect = fl;
    fract = r;
    ior = iOFr;
    pow = p;
    // color = c;
    // position = pos;
    radius = rad;
  }

  int itWorks(){
    return 30;
  }
};

bool testArray(Object** objs) {

cout << objs[1]->kk << endl;

return true;

}
int main(){

  Object * test;
  Sphere s1 = Sphere(1,2,3,4,5,6,7,8);
  test = &s1;

  Object * test2;
  Sphere s2 = Sphere(9,10,11,12,13,14,15,16);
  test2 = &s2;

  Object * objs[2] =  {test,test2};

  testArray(objs);

  /*for(int i = 0; i < 2; i++){
    cout << objs[i]->kk << endl;

    if(objs[i]->itWorks() == 30){
      cout << "it works" << endl;
    }
  }*/

}
