#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>

#define PI 3.14159265
#define kEpsilon 0.001f
using namespace std;

class Point;

class Vector{
    public:
        float x, y, z;

        Vector(float x, float y, float z){
            this->x = x;
            this->y = y;
            this->z = z;
        }
        Vector(float x){
            this->x = x;
            this->y = x;
            this->z = x;
        }
        Vector(){
            this->x = 0;
            this->y = 0;
            this->z = 0;
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

        Point(float x, float y, float z){
            this->x = x;
            this->y = y;
            this->z = z;
        }
        Point(float x){
            this->x = x;
            this->y = x;
            this->z = x;
        }
        Point(){
            this->x = 0;
            this->y = 0;
            this->z = 0;
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
            orthoV = up.cross(orthoU).normalize();
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

        Color(float C){
            this->R = C;
            this->G = C;
            this->B = C;
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

        Color operator * (const float &n) const
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
        float difCo;
        float espCo;
        float ambCo;
        float refCo;
        float tranCo;
        float rugCo;

        explicit Object(Color color, float difCo,
        float espCo, float ambCo, float refCo, float tranCo, float rugCo){
            this->color = color;
            this->difCo = difCo;
            this->espCo = espCo;
            this->ambCo = ambCo;
            this->refCo = refCo;
            this->tranCo = tranCo;
            this->rugCo = rugCo;
        }

        virtual tuple<Point, Vector, float> intersect(Point origin, Vector dir, bool backCulling = true){
            return tuple<Point, Vector, float>{Point(), Vector(), -1};
        }
};

class Sphere: public virtual Object{
    public:
        Point center;
        float radius;
        
        Sphere(Point center, float radius, Color color, float difCo,
        float espCo, float ambCo, float refCo, float tranCo, float rugCo):
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

        tuple<Point, Vector, float> intersect(Point origin, Vector dir, bool backCulling = true){
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
                            
                            return tuple<Point, Vector, float>{inters, normal, t1};
                        }
                        else{
                            Point inters = origin + (dir * t2);
                            Vector normal = (inters - center).normalize();
                            
                            return tuple<Point, Vector, float>{inters, normal, t2};
                        }
                    }

                }

            }

            return tuple<Point, Vector, float>{Point(), Vector(), -1};
        }
};

class Plane: public virtual Object{
    public:
        Point point;
        Vector normal;

        Plane(Point point, Vector normal, Color color, float difCo,
        float espCo, float ambCo, float refCo, float tranCo, float rugCo):
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

        tuple<Point, Vector, float> intersect(Point origin, Vector dir, bool backCulling = true){
            dir.normalize();
            float denom = normal.dot(dir);

            if(abs(denom) > kEpsilon){
                Vector v = point - origin;
                float t = v.dot(normal) / denom;
                if(t >= 0){
                    Point inters = origin + (dir * t);
                    return tuple<Point, Vector, float>{inters, normal, t};
                }
            }
            return tuple<Point, Vector, float>{Point(), Vector(), -1};
        }
};

class Mesh: public Object{
    public:
        int triCount;
        int vertCount;
        vector<Point> vertices;
        vector<tuple<int, int, int>> triangles;
        vector<Vector> triNormals;
        vector<Vector> fullTriNormals;
        vector<Vector> vertNormals;

        Mesh(int triCount, int vertCount, vector<Point> vertices,
        vector<tuple<int, int, int>> triangles, Color color, float difCo, float espCo,
        float ambCo, float refCo, float tranCo, float rugCo):
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

            getNormals();
        }

        void getNormals(){
            for(int i = 0; i<vertCount;i++){
                vertNormals.push_back(Vector());
            }
            
            for (tuple<int, int, int> tri : triangles){
                Point triVerts[3]={vertices[get<0>(tri)],
                                vertices[get<1>(tri)],
                                vertices[get<2>(tri)]};

                Vector vec1 = triVerts[1] - triVerts[0];
                Vector vec2 = triVerts[2] - triVerts[0];

                Vector normal = vec1.cross(vec2);
                triNormals.push_back(normal);

                normal.normalize();

                vertNormals[get<0>(tri)] += normal;
                vertNormals[get<1>(tri)] += normal;
                vertNormals[get<2>(tri)] += normal;
            }

            for(int i = 0; i<vertCount;i++){
                vertNormals[i].normalize();
            }
        }

        tuple<Point, Vector, float> intersect(Point origin, Vector dir, bool backCulling = true){
            float closestDist = numeric_limits<float>::infinity();
            tuple<Point, Vector, float> closestInter = {Point(), Vector(), -1};

            for (int i = 0; i < triangles.size(); i++){
                Vector normal = triNormals[i];
                tuple<int, int, int> triangle = triangles[i];
                Point vert0 = vertices[get<0>(triangle)];
                Point vert1 = vertices[get<1>(triangle)];
                Point vert2 = vertices[get<2>(triangle)];
                float NdotRayDirection = normal.dot(dir);
                float denom = normal.dot(normal);
                float u;
                float v;
                
                if (fabs(NdotRayDirection) > kEpsilon){
                    float d = -(normal.x * vert0.x + normal.y * vert0.y + normal.z * vert0.z);
                    
                    float t = -((normal.x * origin.x + normal.y * origin.y + normal.z * origin.z) + d) / NdotRayDirection;
                    
                    //se t<0 o triangulo está atras do raio
                    if (t < 0 || t >= closestDist)
                        break;

                    Point interPoint = origin + dir * t;

                    //checa se a face do triangulo está na direção da camera
                    if(backCulling){
                        Vector origInter = (origin - interPoint).normalized();
                        float dot = origInter.dot(normal.normalized());
                        if(dot<0) 
                            continue;
                    }
                    
                    //vetor perpendicular ao plano do triangulo
                    Vector interPerp;

                    Vector edge0 = vert1 - vert0; 
                    Vector vp0 = interPoint - vert0;
                    interPerp = edge0.cross(vp0);

                    //checa interseção estado do lado certo desta aresta
                    if (normal.dot(interPerp) < 0)
                        continue;

                    Vector edge1 = vert2 - vert1; 
                    Vector vp1 = interPoint - vert1;
                    interPerp = edge1.cross(vp1);

                    //checa interseção estado do lado certo desta aresta
                    if ((u = normal.dot(interPerp)) < 0)
                        continue;

                    Vector edge2 = vert0 - vert2; 
                    Vector vp2 = interPoint - vert2;
                    interPerp = edge2.cross(vp2);

                    //checa interseção estado do lado certo desta aresta
                    if ((v = normal.dot(interPerp)) < 0)
                        continue;

                    Vector vNormal0 = vertNormals[get<0>(triangle)];
                    Vector vNormal1 = vertNormals[get<1>(triangle)];
                    Vector vNormal2 = vertNormals[get<2>(triangle)];

                    u /= denom;
                    v /= denom;

                    //interpola as normais dos vetores
                    Vector hitNormal = vNormal0 * u + vNormal1 * v +vNormal2 * (1 - u - v); 
                    
                    closestInter = {interPoint, hitNormal.normalized(), t};
                    closestDist = t;
                }
                    
            }
            return closestInter;
        }
};