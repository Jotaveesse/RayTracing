#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <sstream>  

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
        float R, G, B;
        
        Color(float ColorR, float ColorG, float ColorB){
            R = ColorR;
            G = ColorG;
            B = ColorB;
        }

        Color(){
            R = 0;
            G = 0;
            B = 0;
        }

        void normalize(){
            float invLen = 1.0f/255.0f;

            R *= invLen;
            G *= invLen;
            B *= invLen;
        }

        const float operator [] (uint8_t i) const {
                if (i <= 2)
                    return (&R)[i];
                else
                    __throw_out_of_range;
            }
        float& operator [] (uint8_t i) { 
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
            color = inColor;
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
            color = inColor;
        }
};

class Mesh{
    public:
        int triCount;
        int vertCount;
        list<Point> vertices;
        list<int> triangles;
        list<Vector> triNormals;
        list<Vector> vertNormals;

        Color color;
        Mesh(){

        }
};

int main() {
    Camera cam(10,10,10,Vector(0,2,0),Point(0,0,0),Point(10,0,0));
    Color c(123,4,12);

    string data;
    vector<string> dataList;
    string line;
    ifstream inputFile("input.txt");

    while (getline(inputFile, line)) {
        dataList.clear();
        stringstream streamLine(line);

        while (getline(streamLine, data, ' ')) {
            //cout << data << endl;
            dataList.push_back(data);
        }

        if(dataList[0] == "s"){
            Point center(stoi(dataList[1]),
                stoi(dataList[2]),
                stoi(dataList[3]));
            
            float radius = stoi(dataList[4]);

            Color col(stoi(dataList[5]),
                stoi(dataList[6]),
                stoi(dataList[7]));
            
            col.normalize();
            Sphere sph(center, radius, col);
            cout << sph.color.R << endl;
        }

            
    }

    inputFile.close();

    return 0;
};