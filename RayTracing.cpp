#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <sstream>
#include <tuple>

#define PI 3.14159265
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

        Vector& operator += (const Vector &v)
        {
            this->x += v.x;
            this->y += v.y;
            this->z += v.z;
            return *this;
        }

        Vector operator - (const Vector &v) const
        {return Vector(x - v.x, y - v.y, z - v.z);}

        Vector operator * (const float &n) const
        {return Vector(x * n, y * n, z * n);}

        Vector operator / (const float &n) const
        {return Vector(x / n, y / n, z / n);}

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

        Point operator + (const Vector &v) const
        {return Point(x + v.x, y + v.y, z + v.z);}

        Point operator - (const Vector &v) const
        {return Point(x - v.x, y - v.y, z - v.z);}

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

        tuple<Point, Vector, float, Color> intersect(Point origin, Vector dir){
            Vector L = center - origin;
            float lengthL = L.length();
            dir.normalize();
            float tc = L.dot(dir);

            if (tc >= 0){
                float d = sqrt(lengthL*lengthL - tc*tc);

                if (d <= radius){
                    float tlc = sqrt(radius*radius - d*d);
                    float t1 = tc - tlc;
                    float t2 = tc + tlc;
                    
                    if(t1 > 0 || t2 > 0){
                        if (t1 < t2){
                            Point inters = origin + (dir * t1);
                            Vector normal = (inters - center).normalize();
                            
                            return tuple<Point, Vector, float, Color>{inters, normal, t1, color};
                        }
                        else{
                            Point inters = origin + (dir * t2);
                            Vector normal = (inters - center).normalize();
                            
                            return tuple<Point, Vector, float, Color>{inters, normal, t2, color};
                        }
                    }

                }

            }

            return tuple<Point, Vector, float, Color>{Point(), Vector(), -1, Color()};
        }
};

class Plane{
    public:
        Point point;
        Vector normal;
        Color color;
        float difCo;
        float espCo;
        float ambCo;
        float refCo;
        float tranCo;
        float rugCo;

        Plane(Point point, Vector normal, Color color, float difCo,
        float espCo, float ambCo, float refCo, float tranCo, float rugCo){
            this->point = point;
            this->normal = normal;
            this->color = color;
            this->difCo = difCo;
            this->espCo = espCo;
            this->ambCo = ambCo;
            this->refCo = refCo;
            this->tranCo = tranCo;
            this->rugCo = rugCo;
        }

        tuple<Point, Vector, float, Color> intersect(Point origin, Vector dir){
            dir.normalize();
            float denom = normal.dot(dir);

            if(abs(denom) > 0.0001f){
                Vector v = point - origin;
                float t = v.dot(normal) / denom;
                if(t >= 0){
                    Point inters = origin + (dir * t);
                    return tuple<Point, Vector, float, Color>{inters, normal, t, color};
                }
            }
            return tuple<Point, Vector, float, Color>{Point(), Vector(), -1, Color()};
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

            getTriNormals();
        }

        void getTriNormals(){
            for(int i = 0; i<vertCount;i++){
                vertNormals.push_back(Vector());
            }
            
            for (tuple<int, int, int> tri : triangles){
                Point triVerts[3]={vertices[get<0>(tri)],
                                vertices[get<1>(tri)],
                                vertices[get<2>(tri)]};

                Vector vec1 = triVerts[1] - triVerts[0];
                Vector vec2 = triVerts[2] - triVerts[0];

                Vector normal = vec1.cross(vec2).normalize();

                vertNormals[get<0>(tri)] += normal;
                vertNormals[get<1>(tri)] += normal;
                vertNormals[get<2>(tri)] += normal;

                triNormals.push_back(normal);
            }

            for(int i = 0; i<vertCount;i++){
                vertNormals[i].normalize();
            }
        }
};

void trace(Camera cam, Scene scn, vector<Sphere> spheres, vector<Plane> planes){
    ofstream imagePpm("image.ppm");

    imagePpm << "P3\n"<< cam.width << " " << cam.height << "\n255\n";

    Vector t = cam.target - cam.center;
    Vector b = cam.up.cross(t);
    Vector v = t.cross(b);
    float fov = 90;
    
    t.normalize();
    b.normalize();
    v.normalize();

    float gx = cam.distScreen * tan(fov * PI / 180.0 / 2.0);
    float gy = gx * cam.height / cam.width;

    Vector pixWidth = (b * 2 * gx)/(cam.width-1);
    Vector pixHeight = (v * 2 * gy)/(cam.height-1);

    Vector firstPix = t * cam.distScreen - b * gx + v * gy;

    for(int i=0; i<cam.height;i++){
        for(int j=0; j<cam.width;j++){
            tuple<Point, Vector, float, Color> closestInter;
            float closestDist = numeric_limits<float>::infinity();
            Color paintColor = scn.ambient;

            for(Sphere sph : spheres){
                Vector pixVector = firstPix + pixWidth * (j-1) - pixHeight * (i-1);
                tuple<Point, Vector, float, Color> inter = sph.intersect(cam.center, pixVector);
                
                if(get<2>(inter) >= pixVector.length() && get<2>(inter) >= 0 && get<2>(inter) < closestDist){
                    closestDist = get<2>(inter);
                    closestInter = inter;
                    paintColor = get<3>(inter);
                }
            }

            for(Plane pln : planes){
                Vector pixVector = firstPix + pixWidth * (j-1) - pixHeight * (i-1);
                tuple<Point, Vector, float, Color> inter = pln.intersect(cam.center, pixVector);
                
                if(get<2>(inter) >= pixVector.length() && get<2>(inter) >= 0 && get<2>(inter) < closestDist){
                    closestDist = get<2>(inter);
                    closestInter = inter;
                    paintColor = get<3>(inter);
                }
            }

            imagePpm << (int)(paintColor.R*255) << " ";
            imagePpm << (int)(paintColor.G*255) << " ";
            imagePpm << (int)(paintColor.B*255) << " ";
        }
        imagePpm << "\n";
    }

    
    imagePpm.close();

}

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

                tuple<int, int, int> vertIndex{stof(dataList[0])-1,
                                stof(dataList[1])-1,
                                stof(dataList[2])-1};

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
            int width = stoi(dataList[1]);
            int height = stoi(dataList[2]);
            
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
            col.normalize();
            globalScene = new Scene(col, lightList);
        }

        
    }

    inputFile.close();

    trace(*globalCam, *globalScene, sphereList, planeList);

    return 0;
};