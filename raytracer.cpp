#include <cmath>
#include "CImg.h"
#include <limits>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>

const double PI = 3.141592653589793;

using namespace cimg_library;
using namespace std;

struct Vector {
    double x, y, z;

    Vector() {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector(double a, double b, double c) {
        x = a;
        y = b;
        z = c;
    }

    //normalize vector so length is 1
    Vector normalize() {
        double num = sqrt(x * x + y * y + z * z);
        return Vector(x / num, y / num, z / num);
    }

    //overriding opertors to make life easier
    Vector operator+(Vector v) {
        return Vector(x + v.x, y + v.y, z + v.z);
    }

    Vector operator-(Vector v) {
        return Vector(x - v.x, y - v.y, z - v.z);
    }

    Vector operator*(double d) {
        return Vector(x * d, y * d, z * d);
    }

    Vector operator/(double d) {
        return Vector(x / d, y / d, z / d);
    }
};

double dot(Vector v, Vector b) {
    return (v.x * b.x + v.y * b.y + v.z * b.z);
}

struct Vector4 {
    double a, b, c, d;

    Vector4() {
        a = 0;
        b = 0;
        c = 0;
        d = 0;
    }

    Vector4(double aa, double bb, double cc, double dd) {
        a = aa;
        b = bb;
        c = cc;
        d = dd;
    }

    //normalize vector so length is 1
    Vector4 normalize() {
        double num = sqrt(a * a + b * b + c * c + d * d);
        return Vector4(a / num, b / num, c / num, d / num);
    }

    //overriding opertors to make life easier
    Vector4 operator+(Vector4 v) {
        return Vector4(a + v.a, b + v.b, c + v.c, d + v.d);
    }

    Vector4 operator-(Vector4 v) {
        return Vector4(a - v.a, b - v.b, c - v.c, d - v.d);
    }

    Vector4 operator*(double num) {
        return Vector4(a * num, b * num, c * num, d * num);
    }

    Vector4 operator/(double num) {
        return Vector4(a / num, b / num, c / num, d / num);
    }
};


Vector cross(Vector v, Vector b) {
    return Vector(v.y * b.z - v.z * b.y, v.z * b.x - v.x * b.z, v.x * b.y - v.y * b.x);
}

double dot4(Vector4 v1, Vector4 v2) {
    return (v1.a * v2.a + v1.b * v2.b + v1.c * v2.c, v1.d * v2.d);
}

struct matrix4 {
    Vector4 v1;
    Vector4 v2;
    Vector4 v3;
    Vector4 v4;

    matrix4() {
        v1 = Vector4(0, 0, 0, 0);
        v2 = Vector4(0, 0, 0, 0);
        v3 = Vector4(0, 0, 0, 0);
        v4 = Vector4(0, 0, 0, 0);
    }

    matrix4(Vector4 vec1, Vector4 vec2, Vector4 vec3, Vector4 vec4) {
        v1 = vec1;
        v2 = vec2;
        v3 = vec3;
        v4 = vec4;
    }

    matrix4 operator*(matrix4 M) {
        matrix4 finalM = matrix4();
        finalM.v1.a = v1.a * M.v1.a + v1.b * M.v2.a + v1.c * M.v3.a + v1.d * M.v4.a;
        finalM.v1.b = v1.a * M.v1.b + v1.b * M.v2.b + v1.c * M.v3.b + v1.d * M.v4.b;
        finalM.v1.c = v1.a * M.v1.c + v1.b * M.v2.c + v1.c * M.v3.c + v1.d * M.v4.c;
        finalM.v1.d = v1.a * M.v1.d + v1.b * M.v2.d + v1.c * M.v3.d + v1.d * M.v4.d;

        finalM.v2.a = v2.a * M.v1.a + v2.b * M.v2.a + v2.c * M.v3.a + v2.d * M.v4.a;
        finalM.v2.b = v2.a * M.v1.b + v2.b * M.v2.b + v2.c * M.v3.b + v2.d * M.v4.b;
        finalM.v2.c = v2.a * M.v1.c + v2.b * M.v2.c + v2.c * M.v3.c + v2.d * M.v4.c;
        finalM.v2.d = v2.a * M.v1.d + v2.b * M.v2.d + v2.c * M.v3.d + v2.d * M.v4.d;

        finalM.v3.a = v3.a * M.v1.a + v3.b * M.v2.a + v3.c * M.v3.a + v3.d * M.v4.a;
        finalM.v3.b = v3.a * M.v1.b + v3.b * M.v2.b + v3.c * M.v3.b + v3.d * M.v4.b;
        finalM.v3.c = v3.a * M.v1.c + v3.b * M.v2.c + v3.c * M.v3.c + v3.d * M.v4.c;
        finalM.v3.d = v3.a * M.v1.d + v3.b * M.v2.d + v3.c * M.v3.d + v3.d * M.v4.d;

        finalM.v4.a = v4.a * M.v1.a + v4.b * M.v2.a + v4.c * M.v3.a + v4.d * M.v4.a;
        finalM.v4.b = v4.a * M.v1.b + v4.b * M.v2.b + v4.c * M.v3.b + v4.d * M.v4.b;
        finalM.v4.c = v4.a * M.v1.c + v4.b * M.v2.c + v4.c * M.v3.c + v4.d * M.v4.c;
        finalM.v4.d = v4.a * M.v1.d + v4.b * M.v2.d + v4.c * M.v3.d + v4.d * M.v4.d;
        return finalM;
    }
};

matrix4 cameraToWorld(Vector e, Vector d, Vector up) {
    Vector zaxis = e - d;
    zaxis = zaxis.normalize();
    Vector xaxis = cross(up, zaxis);
    xaxis = xaxis.normalize();
    Vector yaxis = cross(zaxis, xaxis);

    matrix4 orientation = matrix4();
    orientation.v1 = Vector4(xaxis.x, yaxis.x, zaxis.x, 0);
    orientation.v2 = Vector4(xaxis.y, yaxis.y, zaxis.y, 0);
    orientation.v3 = Vector4(xaxis.z, yaxis.z, zaxis.z, 0);
    orientation.v4 = Vector4(0, 0, 0, 1);

    matrix4 translation = matrix4();
    translation.v1 = Vector4(1, 0, 0, 0);
    translation.v2 = Vector4(0, 1, 0, 0);
    translation.v3 = Vector4(0, 0, 1, 0);
    translation.v4 = Vector4(-e.x, -e.y, -e.z, 1);

    return (orientation * translation);


}

struct Object {
    double kk, ks, ka, flect, fract, ior;
    int pow;
    Vector color;

    virtual bool intersect(Vector e, Vector d, double &t) {};

    virtual Vector normal(Vector ray) {};
};

struct Sphere : Object {

    double radius;
    Vector origin;

    Sphere(double k1, double k2, double k3, double fl, double r, double iOFr, int p, double rad, const Vector &c,
           Vector o) {
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

    bool intersect(Vector e, Vector d, double &t) {

        double a = dot(d, d);
        double b = 2 * dot(d, (e - origin));
        double c = dot((e - origin), (e - origin)) - radius * radius;
        double discriminant = (b * b) - (4 * a * c);
        if (discriminant < 0) {
            return false;
        } else {
            double num = sqrt(discriminant);
            double t1 = ((-b) + num) / (2 * a);
            double t2 = ((-b) - num) / (2 * a);
            t = (t1 < t2) ? t1 : t2;
            return true;
        }
    }

    Vector normal(Vector ray) {
        return (origin - ray) / radius;
    }

};

struct Plane : Object {
    Vector point;
    Vector n;

    Plane(double k1, double k2, double k3, double fl, double r, double iOFr, int p, const Vector &c, Vector origin,
          Vector norm) {
        kk = k1;
        ks = k2;
        ka = k3;
        flect = 0;
        fract = 0;
        ior = iOFr;
        pow = p;
        color = c;
        point = origin;
        n = norm;
    }

    bool intersect(Vector RayOri, Vector RayDir, double &t) {
        double ln = dot(n, RayDir);
        if (abs(ln) > 1e-6) {
            Vector PROri = point - RayOri;
            t = dot(PROri, n) / ln;
            return (t >= 0);
        }
        return false;
    }

    Vector normal(Vector ray) {
        return n;
    }
};

struct Ray {

    Ray(const Vector &orig, const Vector &dir) {

        ori = orig;
        direc = dir;
        invdir = Vector(1 / direc.x, 1 / direc.y, 1 / direc.z);
        sign[0] = (invdir.x < 0);
        sign[1] = (invdir.y < 0);
        sign[2] = (invdir.z < 0);
    }

    Vector ori;
    Vector direc;
    Vector invdir;
    int sign[3];
};

struct Box : Object {

    Vector origin;
    double size;

    Box(double k1, double k2, double k3, double fl, double r, double iOFr, int p, const Vector &c, const Vector &o,
        double s) {
        kk = k1;
        ks = k2;
        ka = k3;
        flect = fl;
        fract = r;
        ior = iOFr;
        pow = p;
        color = c;
        size = s;
        origin = o;
        // const Vector &vMin, const Vector &vMax
        // bounds[0] = vMin;
        // bounds[1] = vMax;
    }

    bool intersect(Vector e, Vector d, double &t) {

        Vector vMin = Vector(origin.x - size / 2, origin.y - size / 2, origin.z - size / 2);
        Vector vMax = Vector(origin.x + size / 2, origin.y + size / 2, origin.z + size / 2);
        Vector bounds[2] = {vMin, vMax};

        Ray r(e, d);
        float tMin;
        float tMax;
        float tYMin;
        float tYMax;
        float tZMin;
        float tZMax;

        tMin = (bounds[r.sign[0]].x - r.ori.x) * r.invdir.x;
        tMax = (bounds[1 - r.sign[0]].x - r.ori.x) * r.invdir.x;
        tYMin = (bounds[r.sign[1]].y - r.ori.y) * r.invdir.y;
        tYMax = (bounds[1 - r.sign[1]].y - r.ori.y) * r.invdir.y;

        if ((tMin > tYMax) || (tYMin > tMax)) {
            return false;
        }

        if (tYMin > tMin) {
            tMin = tYMin;
        }

        if (tYMax < tMax) {
            tMax = tYMax;
        }

        tZMin = (bounds[r.sign[2]].z - r.ori.z) * r.invdir.z;
        tZMax = (bounds[1 - r.sign[2]].z - r.ori.z) * r.invdir.z;

        if ((tMin > tZMax) || (tZMin > tMax)) {
            return false;
        }

        if (tZMin > tMin) {
            tMin = tZMin;
        }

        if (tZMax < tMax) {
            tMax = tZMax;
        }

        t = tMin;

        if (t < 0) {

            t = tMax;

            if (t < 0) {

                return false;
            }
        }

        return true;
    }

    Vector normal(Vector ray) {
        Vector norm;
        Vector localPoint = ray - origin;
        double min = std::numeric_limits<double>::max();
        double dis = std::abs(size - std::abs(localPoint.x));
        if (dis < min) {
            min = dis;
            norm = Vector(-1, 0, 0);
            if (localPoint.x < 0) {
                norm = Vector(0, 0, 0) - norm;
            }
        }
        dis = std::abs(size - std::abs(localPoint.y));
        if (dis < min) {
            min = dis;
            norm = Vector(0, -1, 0);
            if (localPoint.y < 0) {
                norm = Vector(0, 0, 0) - norm;
            }
        }
        dis = std::abs(size - std::abs(localPoint.z));
        if (dis < min) {
            min = dis;
            norm = Vector(0, 0, -1);
            if (localPoint.z < 0) {
                norm = Vector(0, 0, 0) - norm;
            }
        }
        return norm;
    }

};

struct Cylinder : Object {

    Vector n;
    double radius;
    Vector origin;
    double height;

    Cylinder(double k1, double k2, double k3, double fl, double r, double iOFr, int p, const Vector &c, const Vector &o,
             double rad, double h) {
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
        height = h;
    }

    bool intersect(Vector e, Vector d, double &t) {
        Vector translate = e - origin;

        double a = d.x * d.x + d.z * d.z;
        double b = 2 * translate.x * d.x + 2 * translate.z * d.z;
        double c = translate.x * translate.x + translate.z * translate.z - radius;

        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            return false;
        } else {
            double num = sqrt(discriminant);
            double t1 = ((-b) + num) / (2 * a);
            double t2 = ((-b) - num) / (2 * a);

            if (t1 > t2) {
                double tmp = t1;
                t1 = t2;
                t2 = tmp;
            }

            double y1 = translate.y + t1 * d.y;
            double y2 = translate.y + t2 * d.y;

            double low = -height / 2;
            double high = height / 2;

            if (y1 < low) {
                if (y2 < low) {
                    return false;
                } else {
                    double th = t1 + (t2 - t1) * (y1 + high) / (y1 - y2);
                    if (th <= 0) {
                        return false;
                    }
                    t = th;
                    n = Vector(0, 1, 0);
                    return true;
                }
            } else if (y1 >= low && y1 <= high) {
                if (t1 <= 0) {
                    return false;
                }
                t = t1;
                Vector hitPos = translate + d * t;
                n = Vector(-hitPos.x, 0, -hitPos.z);
                n = n.normalize();
                return true;
            } else if (y1 > high) {
                if (y2 > high) {
                    return false;
                } else {
                    double th = t1 + (t2 - t1) * (y1 + low) / (y1 - y2);
                    if (th <= 0) {
                        return false;
                    }

                    t = th;
                    n = Vector(0, -1, 0);
                    return true;
                }
            }
            return false;
        }
    }

    Vector normal(Vector ray) {
        return n;
    }

};

double isShadow(Object **objs, int size, Vector viewRay, Vector light, Vector normal, int currentObj) {
    Vector l = Vector(-light.x, -light.y, -light.z);
    for (int i = 0; i < size; i++) {
        double distance = std::numeric_limits<double>::infinity();
        bool test = objs[i]->intersect(viewRay, l, distance);
        //to make sure shadow is rendered on correct side of object
        if (i == currentObj) {
        } else if (objs[i]->intersect(viewRay, l, distance) && distance > 0 && objs[i]->fract > 0) {
            return objs[i]->fract;
        } else if (objs[i]->intersect(viewRay, l, distance) && distance > 0) {
            return 0;
        }
    }
    return -1;
}

Vector light = Vector(2, -1, -2);

Vector rayTrace(Vector e, Vector d, Object **objs, int size, int y, int step) {
    Vector returnColor(80, 80, 80);
    int indexTemp;
    Object *hitObj;
    double intensity = 1.5;
    double t = std::numeric_limits<double>::infinity();
    double minT = t;


    for (int i = 0; i < size; i++) {
        //cout << "size: " << size << endl;
        //double t = std::numeric_limits<double>::infinity();
        // double t1 = std::numeric_limits<double>::infinity();
        if (objs[i]->intersect(e, d, t) && t >= 0) {
            //cout << "intersect"<<endl;
            if (t < minT) {
                //cout << "t"<< endl;
                minT = t;
                hitObj = objs[i];
                indexTemp = i;
                Vector viewRay = e + d * t;
                Vector slight = Vector(light.x, light.y, light.z);
                Vector normal = hitObj->normal(viewRay);
                Vector r = normal * (2 * dot(normal, slight)) - slight;
                Vector reflColor;
                double lk;
                double ls;
                double la;
                double nconst = 1.5;
                double c1 = dot(normal.normalize(), viewRay);
                double c2 = sqrt(1 - nconst * nconst * (1 - c1 * c1));
                Vector T = viewRay * nconst + normal.normalize() * (nconst * c1 - c2);
                Vector refractColor = Vector(0, 0, 0);
                double shadowVal = isShadow(objs, size, viewRay, slight, normal, indexTemp);
                if (shadowVal >= 0) {
                    lk = (shadowVal * 0.3 + 0.5) * hitObj->kk * (intensity) *
                         std::max(0.0, dot(normal.normalize(), slight.normalize()));
                    ls = (shadowVal * 0.3 + 0.5) * hitObj->ks * (intensity) *
                         pow(std::max(0.0, dot(normal.normalize(), r.normalize())), hitObj->pow);
                    la = (shadowVal * 0.3 + 0.5) * hitObj->ka * (intensity);
                } else {
                    lk = hitObj->kk * (intensity) * std::max(0.0, dot(normal.normalize(), slight.normalize()));
                    ls = hitObj->ks * (intensity) *
                         pow(std::max(0.0, dot(normal.normalize(), r.normalize())), hitObj->pow);
                    la = hitObj->ka * (intensity);
                }
                Vector refractNormal = normal;
                Vector frenselNormal = normal;
                double kr;
                double fcosi = dot(frenselNormal, d);
                if (fcosi > 1) {
                    fcosi = 1;
                }
                if (fcosi < -1) {
                    fcosi = -1;
                }
                double fetai = 1;
                double fetat = hitObj->ior;

                if (fcosi > 0) {
                    double ftemp = fetai;
                    fetai = fetat;
                    fetat = ftemp;
                }
                double m = 1 - fcosi * fcosi;
                if (m < 0) {
                    m = 0;
                }
                double sint = (1 / (fetai / fetat)) * sqrt(m);

                if (sint >= 1) {
                    kr = 1;
                } else {
                    double n = 1 - sint * sint;
                    if (n < 0) {
                        n = 0;
                    }
                    double cost = sqrt(n);
                    fcosi = abs(fcosi);
                    double Rs = ((fetat * fcosi) - (fetat * cost)) / ((fetat * fcosi) + (fetai * cost));
                    double Rp = ((fetai * fcosi) - (fetat * cost)) / ((fetai * fcosi) + (fetat * cost));
                    kr = (Rs * Rs + Rp * Rp) / 2;
                }
                if (step < 50) {
                    //cout << "refract"<<endl;
                    reflColor = rayTrace(viewRay - normal * 1e-4,
                                         d - normal.normalize() * 2 * dot(d, normal.normalize()), objs, size, y,
                                         step + 1);
                    if (hitObj->fract > 0 && kr < 1) {
                        returnColor = Vector(
                                (hitObj->color.x * la + 255 * ls + hitObj->color.x * lk) + (kr) * reflColor.x,
                                (hitObj->color.y * la + 255 * ls + hitObj->color.y * lk) + (kr) * reflColor.y,
                                (hitObj->color.z * la + 255 * ls + hitObj->color.z * lk) + (kr) * reflColor.z);
                    } else {
                        returnColor = Vector(
                                (hitObj->color.x * la + 255 * ls + hitObj->color.x * lk) + hitObj->flect * reflColor.x,
                                (hitObj->color.y * la + 255 * ls + hitObj->color.y * lk) + hitObj->flect * reflColor.y,
                                (hitObj->color.z * la + 255 * ls + hitObj->color.z * lk) + hitObj->flect * reflColor.z);
                    }
                    if (hitObj->fract > 0 && kr < 1) {
                        Vector I = d;
                        double cosi = dot(I, refractNormal);
                        if (cosi > 1) {
                            cosi = 1;
                        }
                        if (cosi < -1) {
                            cosi = -1;
                        }
                        double etai = 1;
                        double etat = hitObj->ior;
                        if (cosi < 0) {
                            cosi = -cosi;
                        } else {
                            double temp = etai;
                            etai = etat;
                            etat = temp;
                            refractNormal = Vector(0, 0, 0) - refractNormal;
                        }
                        double eta = 1 / (etai / etat);
                        double k = 1 - eta * eta * (1 - cosi * cosi);
                        Vector refdir;
                        if (k < 0) {
                            refdir = Vector(0, 0, 0);
                        } else {
                            refdir = I * eta + refractNormal * (eta * cosi - sqrt(k));
                        }
                        //Vector normRefdir = refdir.normalize();
                        refractColor = rayTrace(viewRay - refractNormal * 1e-4, refdir.normalize(), objs, size, y,
                                                step + 1);
                        returnColor = Vector(returnColor.x + (1 - kr) * hitObj->fract * refractColor.x,
                                             returnColor.y + (1 - kr) * hitObj->fract * refractColor.y,
                                             returnColor.z + (1 - kr) * hitObj->fract * refractColor.z);
                    }
                } else {
                    returnColor = Vector(
                            (hitObj->color.x * la + 255 * ls + hitObj->color.x * lk),
                            (hitObj->color.y * la + 255 * ls + hitObj->color.y * lk),
                            (hitObj->color.z * la + 255 * ls + hitObj->color.z * lk));
                    // return rgb;
                }
            }
        }
    }
    return returnColor;
}

int imageWidth = 1000;
int imageHeight = 1000;
double intensity = 1;
Vector black(0, 0, 0);
Vector green(15, 75, 15);
Vector red(120, 10, 10);
Vector gray(120, 120, 120);
Vector orange(209, 107, 0);
Vector blue(16, 59, 147);

string pColor = "gray";

Sphere s1 = Sphere(1, 1, .5, 0, 0, 1.5, 100, .5, orange, Vector(-2, .5, -9));
Object *test1 = &s1;

Sphere s2 = Sphere(1, 0, .5, 0, 0, 1.5, 5, 2, blue, Vector(5, 2, -20));
Object *test2 = &s2;

Plane p3 = Plane(1, 1, 1, 1, 0, 0, 100, gray, Vector(0, 0, -17), Vector(0, 1, 0));
Object *test3 = &p3;

Box b4 = Box(1, 1, .5, .8, 0, 1, 100, black, Vector(-2, .751, -17), 1.5);
Object *test4 = &b4;

Cylinder c5 = Cylinder(1, 1, 1, 0, .8, 1.3, 100, black, Vector(-1.2, .751, -5), .25, 1.5);
Object *test5 = &c5;

Sphere s6 = Sphere(1, 1, .5, .8, 0, 1.5, 100, 1, black, Vector(1.5, 1, -16));
Object *test6 = &s6;

Cylinder c7 = Cylinder(1, 1, .5, 0, 0, 1.3, 100, orange, Vector(1.4, .501, -7), .5, 1);
Object *test7 = &c7;

Sphere s8 = Sphere(1, 1, .5, .2, 0, 1.5, 100, .75, blue, Vector(0, .75, -10));
Object *test8 = &s8;

Sphere s9 = Sphere(1, 1, .5, 0, 1, 1.5, 100, .25, black, Vector(.75, 1, -3));
Object *test9 = &s9;

Box b10 = Box(1, 1, .5, 0, 0, 1.5, 100, orange, Vector(-1.5, 5, -13), 1);
Object *test10 = &b10;

Object *objs[10] = {test1, test2, test3, test4, test5, test6, test7, test8, test9, test10};

void *change_color(int i) {
    unsigned char textcolor[] = {255, 255, 255};
    int numItems = 6;
    string it;
    ostringstream convert;
    convert << "Current Color: " << pColor << "\nblack" << "\ngreen" << "\nred" << "\ngray" << "\norange" << "\nblue";
    it = convert.str();

    const char *items = it.c_str();
    int height = numItems * 18 + 30;
    CImg<unsigned char> menu(300, height, 1, 3);
    menu.draw_text(1, 1, items, textcolor, 0, 1, 18);
    CImgDisplay disp(menu, "Menu");
    int s_number = 0;


    while (!disp.is_closed() && !disp.is_keyESC()) {

        if (disp.button()) {

            s_number = 1 + (disp.mouse_y()) / 18;


            disp.set_button();
            if (s_number == 2) {
                objs[i]->color = black;
                pColor = "black";
            } else if (s_number == 3) {
                objs[i]->color = green;
                pColor = "green";
            } else if (s_number == 4) {
                objs[i]->color = red;
                pColor = "red";
            } else if (s_number == 5) {
                objs[i]->color = gray;
                pColor = "gray";
            } else if (s_number == 6) {
                objs[i]->color = orange;
                pColor = "orange";
            } else if (s_number == 7) {
                objs[i]->color = blue;
                pColor = "blue";
            }
        }
        menu.assign(300, height, 1, 3, 0);
        convert.str(std::string());
        convert << "Current Color: " << pColor << "\nblack" << "\ngreen" << "\nred" << "\ngray" << "\norange"
                << "\nblue";
        it = convert.str();

        items = it.c_str();
        menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

        disp.display(menu);
    }

    return 0;
}

void *render() {

    //Compute u,v,w basis vectors
    //Creating blank 256x256 image
    CImg<unsigned char> img(imageWidth, imageHeight, 1, 3, 0);
    //for each pixel
    for (int y = 0; y < imageHeight; y++) {
        for (int x = 0; x < imageWidth; x++) {
            //pixelNDC is the normalized pixel position
            //NDC is Normalized Device Coordinates
            //shifted 0.5 pixels because we want to be in the middle of the pixels
            //range of 0-1
            double PixelNDCx = (x + 0.5) / imageWidth;
            double PixelNDCy = (y + 0.5) / imageHeight;

            //maps the coordinates so that the center of the screen is origin
            double PixelScreenx = 2 * PixelNDCx - 1;
            double PixelScreeny = 1 - 2 * PixelNDCy;

            //adjusting based on aspect ratio and also field of view angle
            double FOVangle = 45;
            FOVangle = (FOVangle * PI) / 180;
            float ImageAspectRatio = imageWidth / (float) imageHeight;

            //We are now in Camera Space :)
            double PixelCamerax = (PixelScreenx) * ImageAspectRatio * tan(FOVangle / 2);
            double PixelCameray = (PixelScreeny) * tan(FOVangle / 2);

            //final coordinate of the pixel on the image plane is
            //(PixelCamerax, PixelCameray);
            //ray origin
            Vector o = Vector(0, 2, 0);
            double angleY = .2;
            //ray direction
            Vector dir = Vector(PixelCamerax, PixelCameray - angleY, -1) - o;
            dir = dir.normalize();

            matrix4 cameraToWorldM = cameraToWorld(o, dir, Vector(0, 1, 0));

            Vector4 dir4 = Vector4(PixelCamerax, PixelCameray - angleY, -1, 0);
            matrix4 dirM = matrix4(dir4, Vector4(), Vector4(), Vector4());

            matrix4 newDirM = cameraToWorldM * dirM;

            dir = Vector(newDirM.v1.a, newDirM.v1.b, newDirM.v1.c);
            dir = dir.normalize();

            Vector mycolor = rayTrace(o, dir, objs, sizeof(objs) / sizeof(objs[0]), y, 0);
            //  std::cout << mycolor << std::endl;
            if (mycolor.x > 255) {
                img(x, y, 0) = 255;
            } else {
                img(x, y, 0) = mycolor.x;
            }

            if (mycolor.y > 255) {
                img(x, y, 1) = 255;
            } else {
                img(x, y, 1) = mycolor.y;
            }

            if (mycolor.z > 255) {
                img(x, y, 2) = 255;
            } else {
                img(x, y, 2) = mycolor.z;
            }

        }
    }
    img.display();
    return 0;
}

void *change_light() {

    unsigned char textcolor[] = {255, 255, 255};
    int numItems = 3;
    string it;
    ostringstream convert;
    convert << "X Value:" << light.x << "\nY Value:" << light.y << "\nZ Value:" << light.z;
    it = convert.str();

    const char *items = it.c_str();
    int height = numItems * 18 + 30;
    CImg<unsigned char> menu(300, height, 1, 3);
    menu.draw_text(1, 1, items, textcolor, 0, 1, 18);
    CImgDisplay disp(menu, "Menu");
    int s_number = 0;
    bool key0, key1, key2, key3, key4, key5, key6, key7, key8, key9, change, neg, negchange;
    neg = true;
    negchange = false;

    while (!disp.is_closed() && !disp.is_keyESC()) {

        if (disp.button()) {

            s_number = 1 + (disp.mouse_y()) / 18;
            disp.set_button();
        }
        if (s_number == 1) {
            int num = light.x;
            if (num < 0) {
                neg = true;
            } else {
                neg = false;
            }
            while (!disp.is_closed() && (s_number == 1)) {
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    if (s_number == 4) {
                        negchange = true;
                        neg = !neg;
                        s_number = 1;
                    }
                    disp.set_button();
                }


                if (disp.is_key0() && !key0) {
                    key0 = true;
                    num = 0;
                } else if (!disp.is_key0() && key0) { key0 = false; }

                if (disp.is_key1() && !key1) {
                    key1 = true;
                    num = 1;
                } else if (!disp.is_key1() && key1) { key1 = false; }

                if (disp.is_key2() && !key2) {
                    key2 = true;
                    num = 2;
                } else if (!disp.is_key2() && key2) { key2 = false; }

                if (disp.is_key3() && !key3) {
                    key3 = true;
                    num = 3;
                } else if (!disp.is_key3() && key3) { key3 = false; }

                if (disp.is_key4() && !key4) {
                    key4 = true;
                    num = 4;
                } else if (!disp.is_key4() && key4) { key4 = false; }

                if (disp.is_key5() && !key5) {
                    key5 = true;
                    num = 5;
                } else if (!disp.is_key5() && key5) { key5 = false; }

                if (disp.is_key6() && !key6) {
                    key6 = true;
                    num = 6;
                } else if (!disp.is_key6() && key6) { key6 = false; }

                if (disp.is_key7() && !key7) {
                    key7 = true;
                    num = 7;
                } else if (!disp.is_key7() && key7) { key7 = false; }

                if (disp.is_key8() && !key8) {
                    key8 = true;
                    num = 8;
                } else if (!disp.is_key8() && key8) { key8 = false; }

                if (disp.is_key9() && !key9) {
                    key9 = true;
                    num = 9;
                } else if (!disp.is_key9() && key9) { key9 = false; }

                if (light.x != num || negchange) {
                    change = true;
                } else {
                    change = false;
                }

                menu.assign(300, height, 1, 3, 0);
                if (change && neg) {
                    light.x = -num;

                } else if (change && !neg) {
                    light.x = num;
                }
                change = false;
                negchange = false;


                convert.str(std::string());
                convert << "X Value:" << light.x << "\nY Value:" << light.y << "\nZ Value:" << light.z
                        << "\nToggle Current Variable Negative";
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);

            }
        } else if (s_number == 2) {
            int num = light.y;
            if (num < 0) {
                neg = true;
            } else {
                neg = false;
            }
            while (!disp.is_closed() && (s_number == 2)) {
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    if (s_number == 4) {
                        negchange = true;
                        neg = !neg;
                        s_number = 2;
                    }
                    disp.set_button();
                }


                if (disp.is_key0() && !key0) {
                    key0 = true;
                    num = 0;
                } else if (!disp.is_key0() && key0) { key0 = false; }

                if (disp.is_key1() && !key1) {
                    key1 = true;
                    num = 1;
                } else if (!disp.is_key1() && key1) { key1 = false; }

                if (disp.is_key2() && !key2) {
                    key2 = true;
                    num = 2;
                } else if (!disp.is_key2() && key2) { key2 = false; }

                if (disp.is_key3() && !key3) {
                    key3 = true;
                    num = 3;
                } else if (!disp.is_key3() && key3) { key3 = false; }

                if (disp.is_key4() && !key4) {
                    key4 = true;
                    num = 4;
                } else if (!disp.is_key4() && key4) { key4 = false; }

                if (disp.is_key5() && !key5) {
                    key5 = true;
                    num = 5;
                } else if (!disp.is_key5() && key5) { key5 = false; }

                if (disp.is_key6() && !key6) {
                    key6 = true;
                    num = 6;
                } else if (!disp.is_key6() && key6) { key6 = false; }

                if (disp.is_key7() && !key7) {
                    key7 = true;
                    num = 7;
                } else if (!disp.is_key7() && key7) { key7 = false; }

                if (disp.is_key8() && !key8) {
                    key8 = true;
                    num = 8;
                } else if (!disp.is_key8() && key8) { key8 = false; }

                if (disp.is_key9() && !key9) {
                    key9 = true;
                    num = 9;
                } else if (!disp.is_key9() && key9) { key9 = false; }

                if (light.y != num || negchange) {
                    change = true;
                } else {
                    change = false;
                }

                menu.assign(300, height, 1, 3, 0);
                if (change && neg) {
                    light.y = -num;

                } else if (change && !neg) {
                    light.y = num;
                }
                change = false;
                negchange = false;


                convert.str(std::string());
                convert << "X Value:" << light.x << "\nY Value:" << light.y << "\nZ Value:" << light.z
                        << "\nToggle Current Variable Negative";
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);

            }
        } else if (s_number == 3) {
            int num = light.z;
            if (num < 0) {
                neg = true;
            } else {
                neg = false;
            }
            while (!disp.is_closed() && (s_number == 3)) {
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    if (s_number == 4) {
                        negchange = true;
                        neg = !neg;
                        s_number = 3;
                    }
                    disp.set_button();
                }


                if (disp.is_key0() && !key0) {
                    key0 = true;
                    num = 0;
                } else if (!disp.is_key0() && key0) { key0 = false; }

                if (disp.is_key1() && !key1) {
                    key1 = true;
                    num = 1;
                } else if (!disp.is_key1() && key1) { key1 = false; }

                if (disp.is_key2() && !key2) {
                    key2 = true;
                    num = 2;
                } else if (!disp.is_key2() && key2) { key2 = false; }

                if (disp.is_key3() && !key3) {
                    key3 = true;
                    num = 3;
                } else if (!disp.is_key3() && key3) { key3 = false; }

                if (disp.is_key4() && !key4) {
                    key4 = true;
                    num = 4;
                } else if (!disp.is_key4() && key4) { key4 = false; }

                if (disp.is_key5() && !key5) {
                    key5 = true;
                    num = 5;
                } else if (!disp.is_key5() && key5) { key5 = false; }

                if (disp.is_key6() && !key6) {
                    key6 = true;
                    num = 6;
                } else if (!disp.is_key6() && key6) { key6 = false; }

                if (disp.is_key7() && !key7) {
                    key7 = true;
                    num = 7;
                } else if (!disp.is_key7() && key7) { key7 = false; }

                if (disp.is_key8() && !key8) {
                    key8 = true;
                    num = 8;
                } else if (!disp.is_key8() && key8) { key8 = false; }

                if (disp.is_key9() && !key9) {
                    key9 = true;
                    num = 9;
                } else if (!disp.is_key9() && key9) { key9 = false; }

                if (light.z != num || negchange) {
                    change = true;
                } else {
                    change = false;
                }

                menu.assign(300, height, 1, 3, 0);
                if (change && neg) {
                    light.z = -num;

                } else if (change && !neg) {
                    light.z = num;
                }
                change = false;
                negchange = false;


                convert.str(std::string());
                convert << "X Value:" << light.x << "\nY Value:" << light.y << "\nZ Value:" << light.z
                        << "\nToggle Current Variable Negative";
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);

            }
        }

        key0 = false;
        key1 = false;
        key2 = false;
        key3 = false;
        key4 = false;
        key5 = false;
        key6 = false;
        key7 = false;
        key8 = false;
        key9 = false;
        neg = false;
    }

    return 0;


}

void *change_obj(int i) {
    unsigned char textcolor[] = {255, 255, 255};
    int numItems = 8;
    string it;
    ostringstream convert;


    convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks << "\nSpecular Power Value: "
            << objs[i]->pow << "\nAmbient Value: "
            << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
            << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
            << objs[i]->ior << "\nChange Color: " << pColor;
    it = convert.str();

    const char *items = it.c_str();
    int height = numItems * 18 + 30;
    CImg<unsigned char> menu(250, height, 1, 3);
    menu.draw_text(1, 1, items, textcolor, 0, 1, 18);
    CImgDisplay disp(menu, "Menu");
    int s_number = 0;
    bool left, right, change;


    while (!disp.is_closed() && !disp.is_keyESC()) {

        if (disp.button()) {

            s_number = 1 + (disp.mouse_y()) / 18;
            disp.set_button();
        }


        if (s_number == 1) {
            double num = objs[i]->kk;

            while (!disp.is_closed() && (s_number == 1)) {

                if (disp.is_keyARROWRIGHT() && !right) {
                    right = true;
                    if (num < 1) {
                        num += .1;
                    }
                } else if (!disp.is_keyARROWRIGHT() && right) { right = false; }

                if (disp.is_keyARROWLEFT() && !left) {
                    left = true;
                    if (num > .1) {
                        num -= .1;
                    }
                } else if (!disp.is_keyARROWLEFT() && left) { left = false; }


                menu.assign(250, height, 1, 3, 0);

                objs[i]->kk = num;

                convert.str(std::string());
                convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                        << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                        << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                        << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                        << objs[i]->ior << "\nChange Color: " << pColor;
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    disp.set_button();
                }

            }
        } else if (s_number == 2) {
            double num = objs[i]->ks;

            while (!disp.is_closed() && (s_number == 2)) {

                if (disp.is_keyARROWRIGHT() && !right) {
                    right = true;
                    if (num < 1) {
                        num += .1;
                    }
                } else if (!disp.is_keyARROWRIGHT() && right) { right = false; }

                if (disp.is_keyARROWLEFT() && !left) {
                    left = true;
                    if (num > .1) {
                        num -= .1;
                    }
                } else if (!disp.is_keyARROWLEFT() && left) { left = false; }


                menu.assign(250, height, 1, 3, 0);

                objs[i]->ks = num;

                convert.str(std::string());
                convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                        << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                        << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                        << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                        << objs[i]->ior << "\nChange Color: " << pColor;
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    disp.set_button();
                }

            }
        } else if (s_number == 3) {
            double num = objs[i]->pow;

            while (!disp.is_closed() && (s_number == 3)) {

                if (disp.is_keyARROWRIGHT() && !right) {
                    right = true;
                    if (num < 100) {
                        num += 10;
                    }
                } else if (!disp.is_keyARROWRIGHT() && right) { right = false; }

                if (disp.is_keyARROWLEFT() && !left) {
                    left = true;
                    if (num > 0) {
                        num -= 10;
                    }
                } else if (!disp.is_keyARROWLEFT() && left) { left = false; }


                menu.assign(250, height, 1, 3, 0);

                objs[i]->pow = num;

                convert.str(std::string());
                convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                        << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                        << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                        << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                        << objs[i]->ior << "\nChange Color: " << pColor;
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    disp.set_button();
                }

            }
        } else if (s_number == 4) {
            double num = objs[i]->ka;

            while (!disp.is_closed() && (s_number == 4)) {

                if (disp.is_keyARROWRIGHT() && !right) {
                    right = true;
                    if (num < .9) {
                        num += .1;
                    }
                } else if (!disp.is_keyARROWRIGHT() && right) { right = false; }

                if (disp.is_keyARROWLEFT() && !left) {
                    left = true;
                    if (num > .1) {
                        num -= .1;
                    }
                } else if (!disp.is_keyARROWLEFT() && left) { left = false; }


                menu.assign(250, height, 1, 3, 0);

                objs[i]->ka = num;

                convert.str(std::string());
                convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                        << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                        << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                        << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                        << objs[i]->ior << "\nChange Color: " << pColor;
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    disp.set_button();
                }

            }
        } else if (s_number == 5) {
            double num = objs[i]->flect;

            while (!disp.is_closed() && (s_number == 5)) {

                if (disp.is_keyARROWRIGHT() && !right) {
                    right = true;
                    if (num < .9) {
                        num += .1;
                    }
                } else if (!disp.is_keyARROWRIGHT() && right) { right = false; }

                if (disp.is_keyARROWLEFT() && !left) {
                    left = true;
                    if (num > .1) {
                        num -= .1;
                    }
                } else if (!disp.is_keyARROWLEFT() && left) { left = false; }


                menu.assign(250, height, 1, 3, 0);

                objs[i]->flect = num;

                convert.str(std::string());
                convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                        << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                        << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                        << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                        << objs[i]->ior << "\nChange Color: " << pColor;
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    disp.set_button();
                }

            }
        } else if (s_number == 6) {
            double num = objs[i]->fract;

            while (!disp.is_closed() && (s_number == 6)) {

                if (disp.is_keyARROWRIGHT() && !right) {
                    right = true;
                    if (num < .9) {
                        num += .1;
                    }
                } else if (!disp.is_keyARROWRIGHT() && right) { right = false; }

                if (disp.is_keyARROWLEFT() && !left) {
                    left = true;
                    if (num > .1) {
                        num -= .1;
                    }
                } else if (!disp.is_keyARROWLEFT() && left) { left = false; }


                menu.assign(250, height, 1, 3, 0);

                objs[i]->fract = num;

                convert.str(std::string());
                convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                        << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                        << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                        << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                        << objs[i]->ior << "\nChange Color: " << pColor;
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    disp.set_button();
                }

            }
        } else if (s_number == 7) {
            double num = objs[i]->ior;

            while (!disp.is_closed() && (s_number == 7)) {

                if (disp.is_keyARROWRIGHT() && !right) {
                    right = true;
                    if (num < 2) {
                        num += .1;
                    }
                } else if (!disp.is_keyARROWRIGHT() && right) { right = false; }

                if (disp.is_keyARROWLEFT() && !left) {
                    left = true;
                    if (num > 1) {
                        num -= .1;
                    }
                } else if (!disp.is_keyARROWLEFT() && left) { left = false; }


                menu.assign(250, height, 1, 3, 0);

                objs[i]->ior = num;

                convert.str(std::string());
                convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                        << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                        << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                        << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                        << objs[i]->ior << "\nChange Color: " << pColor;
                it = convert.str();

                items = it.c_str();
                menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

                disp.display(menu);
                if (disp.button()) {
                    s_number = 1 + (disp.mouse_y()) / 18;
                    disp.set_button();
                }

            }
        } else if (s_number == 8) {
            change_color(i);
            s_number = 0;
            menu.assign(250, height, 1, 3, 0);
            convert.str(std::string());
            convert << "Diffuse Value: " << objs[i]->kk << "\nSpecular Value: " << objs[i]->ks
                    << "\nSpecular Power Value: " << objs[i]->pow << "\nAmbient Value: "
                    << objs[i]->ka << "\nReflective Value: " << objs[i]->flect
                    << "\nTransparency Value: " << objs[i]->fract << "\nIndex of Refraction Value: "
                    << objs[i]->ior << "\nChange Color: " << pColor;
            it = convert.str();

            items = it.c_str();
            menu.draw_text(1, 1, items, textcolor, 0, 1, 18);

            disp.display(menu);
        }
        left = false;
        right = false;
    }

    return 0;
}

void menu_item(const unsigned int demo_number) {
    switch (demo_number) {
        case 1:
            render();
            break;
        case 2:
            change_light();
            break;
        case 3:
            change_color(2);
            break;
        case 4:
            change_obj(0);
            break;
        case 5:
            change_obj(1);
            break;
        case 6:
            change_obj(3);
            break;
        case 7:
            change_obj(4);
            break;
        case 8:
            change_obj(5);
            break;
        case 9:
            change_obj(6);
            break;
        case 10:
            change_obj(7);
            break;
        case 11:
            change_obj(8);
            break;
        case 12:
            change_obj(9);
            break;
        default:
            break;
    }
}

int main(int argc, char **argv) {
    unsigned int demo_number = cimg_option("-run", 0, 0);
    if (demo_number) menu_item(demo_number);
    else {
        //render menu

        unsigned char textcolor[] = {255, 255, 255};
        int numItems = 11;
        const char *const items = "Render\n"
                "Change Light Direction\n"
                "Change Plane Color\n"
                "Sphere1\n"
                "Cube1\n"
                "Cylinder1\n"
                "Sphere2\n"
                "Cylinder2\n"
                "Sphere3\n"
                "Sphere4\n"
                "Cube2\n";
        CImg<unsigned char> menu(200, numItems * 18 + 30, 1, 3);
        menu.draw_text(1, 1, items, textcolor, 0, 1, 18);
        CImgDisplay disp(menu, "Menu");


        while (!disp.is_closed() && !disp.is_keyESC()) {
            if (disp.button()) {
                demo_number = 1 + (disp.mouse_y()) / 18;
                disp.set_button();
            }
            menu_item(demo_number);
            demo_number = 0;
        }


    }

    return 0;
}
