#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <sstream>
#include <tuple>

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

        Camera(int height, int width, float distScreen, Vector up,
        Point center, Point target){
            this->height = height;
            this->width = width;
            this->distScreen = distScreen;
            this->up = up;
            this->center = center;
            this->target = target;

            orthoU = (target - center).normalize();
            orthoV = orthoU.cross(up).normalize();
            orthoW = orthoU.cross(orthoV);
        }
};

class Color{
    public:
        float R, G, B;
        
        Color(float R, float G, float B){
            this->R = R;
            this->G = G;
            this->B = B;
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
                    throw out_of_range("Out of range");
            }
        float& operator [] (uint8_t i) { 
            if (i <= 2)
                return (&R)[i];
            else
                throw out_of_range("Out of range");
            return (&R)[0];
        }
};

class Light{
    public:
        Point center;
        Color color;

        Light(Point center, Color color){
            this->color = color;
            this->center = center;
        }
};

class Scene{
    public:
        Color ambient;
        vector<Light> lights = {};
        Scene(Color ambient, vector<Light> lights){
            this->ambient = ambient;
            this->lights = lights;
        }
};

class Sphere{
    public:
        Point center;
        float radius;
        Color color;
        float difCo;
        float espCo;
        float ambCo;
        float refCo;
        float tranCo;
        float rugCo;

        Sphere(Point center, float radius, Color color, float difCo,
        float espCo, float ambCo, float refCo, float tranCo, float rugCo){
            this->center = center;
            this->radius = radius;
            this->color = color;
            this->difCo = difCo;
            this->espCo = espCo;
            this->ambCo = ambCo;
            this->refCo = refCo;
            this->tranCo = tranCo;
            this->rugCo = rugCo;
        }
};

class Plane{
    public:
        Point point;
        Vector vector;
        Color color;
        float difCo;
        float espCo;
        float ambCo;
        float refCo;
        float tranCo;
        float rugCo;

        Plane(Point point, Vector vector, Color color, float difCo,
        float espCo, float ambCo, float refCo, float tranCo, float rugCo){
            this->point = point;
            this->vector = vector;
            this->color = color;
            this->difCo = difCo;
            this->espCo = espCo;
            this->ambCo = ambCo;
            this->refCo = refCo;
            this->tranCo = tranCo;
            this->rugCo = rugCo;
        }
};

class Mesh{
    public:
        int triCount;
        int vertCount;
        vector<Point> vertices;
        vector<tuple<int, int, int>> triangles;
        vector<Vector> triNormals;
        vector<Vector> vertNormals;
        Color color;
        float difCo;
        float espCo;
        float ambCo;
        float refCo;
        float tranCo;
        float rugCo;

        Mesh(int triCount, int vertCount, vector<Point> vertices,
        vector<tuple<int, int, int>> triangles, Color color, float difCo, float espCo,
        float ambCo, float refCo, float tranCo, float rugCo){
            this->triCount = triCount;
            this->vertCount = vertCount;
            this->vertices = vertices;
            this->triangles = triangles;
            this->color = color;
            this->difCo = difCo;
            this->espCo = espCo;
            this->ambCo = ambCo;
            this->refCo = refCo;
            this->tranCo = tranCo;
            this->rugCo = rugCo;
        }
};

int main() {
    vector<Sphere> sphereList;
    vector<Plane> planeList;
    vector<Mesh> meshList;
    vector<Light> lightList;

    Camera *globalCam;
    Scene *globalScene;

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
            Point center(stof(dataList[1]),
                stof(dataList[2]),
                stof(dataList[3]));
            
            float radius = stof(dataList[4]);

            Color col(stof(dataList[5]),
                stof(dataList[6]),
                stof(dataList[7]));
            col.normalize();

            float difCo = stof(dataList[8]);
            float espCo = stof(dataList[9]);
            float ambCo = stof(dataList[10]);
            float refCo = stof(dataList[11]);
            float tranCo = stof(dataList[12]);
            float rugCo = stof(dataList[13]);

            Sphere sph(center, radius, col, difCo, espCo, ambCo,
                        refCo, tranCo, rugCo);

            sphereList.push_back(sph);
        }
        else if(dataList[0] == "p"){
            Point p(stof(dataList[1]),
                stof(dataList[2]),
                stof(dataList[3]));
            
            Vector normal(stof(dataList[4]),
                stof(dataList[5]),
                stof(dataList[6]));

            Color col(stof(dataList[7]),
                stof(dataList[8]),
                stof(dataList[9]));
            col.normalize();

            float difCo = stof(dataList[10]);
            float espCo = stof(dataList[11]);
            float ambCo = stof(dataList[12]);
            float refCo = stof(dataList[13]);
            float tranCo = stof(dataList[14]);
            float rugCo = stof(dataList[15]);

            Plane pln(p, normal, col, difCo, espCo, ambCo,
                        refCo, tranCo, rugCo);

            planeList.push_back(pln);
        }
        else if(dataList[0] == "t"){
            int triCount = stoi(dataList[1]);
            int vertCount = stoi(dataList[2]);
            vector<Point> vertices;
            vector<tuple<int, int, int>> triangles;

            for (int i=0;i<vertCount;i++){
                getline(inputFile, line);
                dataList.clear();
                stringstream streamLine(line);

                while (getline(streamLine, data, ' ')) {
                    dataList.push_back(data);
                }

                Point vertPoint(stof(dataList[0]),
                                stof(dataList[1]),
                                stof(dataList[2]));

                vertices.push_back(vertPoint);
            }

            for (int j=0;j<triCount;j++){
                getline(inputFile, line);
                dataList.clear();
                stringstream streamLine(line);

                while (getline(streamLine, data, ' ')) {
                    dataList.push_back(data);
                }

                tuple<int, int, int> vertIndex{stof(dataList[0]),
                                stof(dataList[1]),
                                stof(dataList[2])};

                triangles.push_back(vertIndex);
            }

            getline(inputFile, line);
            dataList.clear();
            stringstream streamLine(line);

            while (getline(streamLine, data, ' ')) {
                dataList.push_back(data);
            }

            Color col(stof(dataList[0]),
                stof(dataList[1]),
                stof(dataList[2]));
            col.normalize();
            float difCo = stof(dataList[3]);
            float espCo = stof(dataList[4]);
            float ambCo = stof(dataList[5]);
            float refCo = stof(dataList[6]);
            float tranCo = stof(dataList[7]);
            float rugCo = stof(dataList[8]);

            Mesh mesh(triCount, vertCount, vertices, triangles, col,
            difCo, espCo, ambCo, refCo, tranCo, rugCo);

            meshList.push_back(mesh);
        }
        else if(dataList[0] == "c"){
            
            int height = stoi(dataList[1]);
            int width = stoi(dataList[2]);
            float distScreen = stof(dataList[3]);
            Vector up(stof(dataList[4]),
                stof(dataList[5]),
                stof(dataList[6])); 
            
            Point center(stof(dataList[7]),
                stof(dataList[8]),
                stof(dataList[9]));
            
            Point focus(stof(dataList[10]),
                stof(dataList[11]),
                stof(dataList[12]));

            globalCam = new Camera(height, width, distScreen, up, center, focus);
        }
        else if(dataList[0] == "l"){
            Point center(stof(dataList[1]),
                stof(dataList[2]),
                stof(dataList[3]));
            
            Color intensity(stof(dataList[4]),
                stof(dataList[5]),
                stof(dataList[6]));

            Light light(center, intensity);
            lightList.push_back(light);
        }
        else if(dataList[0] == "a"){
            Color col(stof(dataList[1]),
                stof(dataList[2]),
                stof(dataList[3]));

            globalScene = new Scene(col, lightList);
        }

            
    }

    inputFile.close();

    return 0;
};