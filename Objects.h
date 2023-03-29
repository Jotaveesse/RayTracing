#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>

#define PI 3.14159265
using namespace std;

class Vector{
    public:
        double x, y, z;

        Vector(double xx, double yy, double zz){
            x = xx;
            y = yy;
            z = zz;
        }
        Vector(double xx){
            x = xx;
            y = xx;
            z = xx;
        }
        Vector(){
            x = 0;
            y = 0;
            z = 0;
        }

        double length(){
            return sqrt(x * x + y * y + z * z); 
        }

        double sqrdLength(){
            return (x * x + y * y + z * z); 
        }

        Vector normalize(){
            double len = length();
            if(len>0){
                double invLen = 1/len;
                
                x *= invLen;
                y *= invLen;
                z *= invLen;
            }

            return *this;
        }

        Vector normalized(){
            double len = length();
            if(len>0){
                double invLen = 1/len;
                Vector vec(x * invLen, y * invLen, z * invLen);
                return vec;
            }

            return Vector();
        }

        double dot(Vector vec){
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

        Vector operator * (const double &n) const
        {return Vector(x * n, y * n, z * n);}

        Vector operator / (const double &n) const
        {return Vector(x / n, y / n, z / n);}

        const double operator [] (uint8_t i) const { return (&x)[i]; }
        double& operator [] (uint8_t i) { return (&x)[i]; }

        friend std::ostream& operator << (std::ostream &s, const Vector &v)
        {
            return s << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
        }
};

class Point{
    public:
        double x, y, z;

        Point(double xx, double yy, double zz){
            x = xx;
            y = yy;
            z = zz;
        }
        Point(double xx){
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

        Vector operator * (const double &n) const
        {return Vector(x * n, y * n, z * n);}

        const double operator [] (uint8_t i) const { return (&x)[i]; }
        double& operator [] (uint8_t i) { return (&x)[i]; }

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

        double distScreen;
        int height;
        int width;

        Camera(int height, int width, double distScreen, Vector up,
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
        double R, G, B;
        
        Color(double R, double G, double B){
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
            double invLen = 1.0f/255.0f;

            R *= invLen;
            G *= invLen;
            B *= invLen;
        }

        void clamp(){
            if(R<0)
                R = 0;
            if(G<0)
                G = 0;
            if(B<0)
                B = 0;
            if(R>1)
                R = 1;
            if(G>1)
                G = 1;
            if(B>1)
                B = 1;
        }

        const double operator [] (uint8_t i) const {
                if (i <= 2)
                    return (&R)[i];
                else
                    throw out_of_range("Out of range");
            }
        double& operator [] (uint8_t i) { 
            if (i <= 2)
                return (&R)[i];
            else
                throw out_of_range("Out of range");
            return (&R)[0];
        }

        Color operator * (const double &n) const
        {return Color(R * n,G * n, B * n);}

        Color operator + (const Color &c) const
        {return Color(R + c.R, G + c.G, B + c.B);}

        Color operator * (const Color &c) const
        {return Color(R * c.R, G * c.G, B * c.B);}
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

class Object{
    public:
        Color color;
        double difCo;
        double espCo;
        double ambCo;
        double refCo;
        double tranCo;
        double rugCo;

        explicit Object(Color color, double difCo,
        double espCo, double ambCo, double refCo, double tranCo, double rugCo){
            this->color = color;
            this->difCo = difCo;
            this->espCo = espCo;
            this->ambCo = ambCo;
            this->refCo = refCo;
            this->tranCo = tranCo;
            this->rugCo = rugCo;
        }

        virtual tuple<Point, Vector, double> intersect(Point origin, Vector dir){
            return tuple<Point, Vector, double>{Point(), Vector(), -1};
        }
};

class Sphere: public virtual Object{
    public:
        Point center;
        double radius;
        
        Sphere(Point center, double radius, Color color, double difCo,
        double espCo, double ambCo, double refCo, double tranCo, double rugCo):
        Object(color, difCo, espCo, ambCo, refCo, tranCo, rugCo){
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

        tuple<Point, Vector, double> intersect(Point origin, Vector dir){
            Vector L = center - origin;
            double lengthL = L.length();
            dir.normalize();
            double tc = L.dot(dir);

            if (tc >= 0){
                double d = lengthL*lengthL - tc*tc;

                if (d <= radius*radius){
                    double tlc = sqrt(radius*radius - d);
                    double t1 = tc - tlc;
                    double t2 = tc + tlc;
                    
                    if(t1 > 0 || t2 > 0){
                        if (t1 < t2){
                            Point inters = origin + (dir * t1);
                            Vector normal = (inters - center).normalize();
                            
                            return tuple<Point, Vector, double>{inters, normal, t1};
                        }
                        else{
                            Point inters = origin + (dir * t2);
                            Vector normal = (inters - center).normalize();
                            
                            return tuple<Point, Vector, double>{inters, normal, t2};
                        }
                    }

                }

            }

            return tuple<Point, Vector, double>{Point(), Vector(), -1};
        }
};

class Plane: public virtual Object{
    public:
        Point point;
        Vector normal;

        Plane(Point point, Vector normal, Color color, double difCo,
        double espCo, double ambCo, double refCo, double tranCo, double rugCo):
        Object(color, difCo, espCo, ambCo, refCo, tranCo, rugCo){
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

        tuple<Point, Vector, double> intersect(Point origin, Vector dir){
            dir.normalize();
            double denom = normal.dot(dir);

            if(abs(denom) > 0.0001f){
                Vector v = point - origin;
                double t = v.dot(normal) / denom;
                if(t >= 0){
                    Point inters = origin + (dir * t);
                    return tuple<Point, Vector, double>{inters, normal, t};
                }
            }
            return tuple<Point, Vector, double>{Point(), Vector(), -1};
        }
};

class Mesh: public Object{
    public:
        int triCount;
        int vertCount;
        vector<Point> vertices;
        vector<tuple<int, int, int>> triangles;
        vector<Vector> triNormals;
        vector<Vector> vertNormals;

        Mesh(int triCount, int vertCount, vector<Point> vertices,
        vector<tuple<int, int, int>> triangles, Color color, double difCo, double espCo,
        double ambCo, double refCo, double tranCo, double rugCo):
        Object(color, difCo, espCo, ambCo, refCo, tranCo, rugCo){
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