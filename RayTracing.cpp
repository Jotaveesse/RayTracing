#include <iostream>
#include <cmath>
#include <list>
using namespace std;

class Vector{
    public:
        float x, y, z;

        Vector(float xx, float yy, float zz){
            x = xx;
            y = yy;
            z = zz;
        }
        Vector(float xx){
            x = xx;
            y = xx;
            z = xx;
        }
        Vector(){
            x = 0;
            y = 0;
            z = 0;
        }

        float length(){
            return sqrt(x * x + y * y + z * z); 
        }

        float sqrdLength(){
            return (x * x + y * y + z * z); 
        }

        Vector normalize(){
            float len = length();
            if(len>0){
                float invLen = 1/len;

                x *= invLen;
                y *= invLen;
                z *= invLen;
            }

            return *this;
        }

        Vector normalized(){
            float len = length();
            if(len>0){
                float invLen = 1/len;
                Vector vec(x * invLen, y * invLen, z * invLen);
                return vec;
            }

            return Vector();
        }

        float dot(Vector vec){
            return x * vec.x + y * vec.y + z * vec.z; 
        }

        Vector cross(Vector vec){
            return Vector(
            y * vec.z - z * vec.y, 
            z * vec.x - x * vec.z, 
            x * vec.y - y * vec.x);
        }

        Vector operator + (const Vector &v) const
        {return Vector(x + v.x, y + v.y, z + v.z);}

        Vector operator - (const Vector &v) const
        {return Vector(x - v.x, y - v.y, z - v.z);}

        Vector operator * (const float &n) const
        {return Vector(x * n, y * n, z * n);}

        const float operator [] (uint8_t i) const { return (&x)[i]; }
        float& operator [] (uint8_t i) { return (&x)[i]; }

        friend std::ostream& operator << (std::ostream &s, const Vector &v)
        {
            return s << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
        }
};

class Point{
    public:
        float x, y, z;

        Point(float xx, float yy, float zz){
            x = xx;
            y = yy;
            z = zz;
        }
        Point(float xx){
            x = xx;
            y = xx;
            z = xx;
        }
        Point(){
            x = 0;
            y = 0;
            z = 0;
        }

        Vector operator + (const Point &v) const
        {return Vector(x + v.x, y + v.y, z + v.z);}

        Vector operator - (const Point &v) const
        {return Vector(x - v.x, y - v.y, z - v.z);}

        Vector operator * (const float &n) const
        {return Vector(x * n, y * n, z * n);}

        const float operator [] (uint8_t i) const { return (&x)[i]; }
        float& operator [] (uint8_t i) { return (&x)[i]; }

        friend std::ostream& operator << (std::ostream &s, const Point &v)
        {
            return s << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
        }
};



class Camera{
    
    public:
        Point center;
        Point target;
        Vector up;

        Vector orthoW;
        Vector orthoV;
        Vector orthoU;

        float distScreen;
        int height;
        int width;

        Camera(int inHeight, int inWidth, float inDistScreen, Vector inUp, Point inCenter, Point inTarget){
            height = inHeight;
            width = inWidth;
            distScreen = inDistScreen;
            up = inUp;
            center = inCenter;
            target = inTarget;

            orthoU = (target - center).normalize();
            orthoV = orthoU.cross(up).normalize();
            orthoW = orthoU.cross(orthoV);
        }
};

class Color{
    public:
        int R, G, B;
        
        Color(int ColorR, int ColorG, int ColorB){
            R = ColorR;
            G = ColorG;
            B = ColorB;
            cout << "uhh";
        }

        Color(){
            R = 0;
            G = 0;
            B = 0;
        }

        Color normalized(){
            float invLen = 1/255;

            return Color(R * invLen, G * invLen, B * invLen);
        }

        const int operator [] (uint8_t i) const {
                if (i <= 2)
                    return (&R)[i];
                else
                    __throw_out_of_range;
            }
        int& operator [] (uint8_t i) { 
            if (i <= 2)
                return (&R)[i];
            else
                __throw_out_of_range;
            return (&R)[0];
        }
};

class Light{
    public:
        Point center;
        Color color;

        Light(Point center, Color color){
            color = color;
            center = center;
        }
};

class Scene{
    public:
        Color ambient;
        list<Light> lights = {};
        Scene(Color ambient, list<Light> lights){
            ambient = ambient;
            lights= lights;
        }
};

class Sphere{
    public:
        Point center;
        float radius;
        Color color;

        Sphere(Point inCenter, float inRadius, Color inColor){
            center = inCenter;
            radius = inRadius;
            color = inColor.normalized();
        }
};

class Plane{
    public:
        Point point;
        Vector vector;
        Color color;

        Plane(Point inPoint, Vector inVector, Color inColor){
            point = inPoint;
            vector = inVector;
            color = inColor.normalized();
        }
};

class Mesh{
    public:
        int triCount;
        int vertCount;
        list<Point> vertices;
        list<tuple<int, int, int>> triangles;
        list<Vector> triNormals;
        list<Vector> vertNormals;

        Color color;
        Mesh(){

        }
};

int main() {
    Camera cam(10,10,10,Vector(0,2,0),Point(0,0,0),Point(10,0,0));
    Color c(123,4,12);
    cout << c[2];
    return 0;
};