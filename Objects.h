#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>

#define PI 3.14159265
using namespace std;

class Point;

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

        virtual tuple<Point, Vector, float> intersect(Point origin, Vector dir){
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

        tuple<Point, Vector, float> intersect(Point origin, Vector dir){
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

        tuple<Point, Vector, float> intersect(Point origin, Vector dir){
            dir.normalize();
            float denom = normal.dot(dir);

            if(abs(denom) > 0.0001f){
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
                fullTriNormals.push_back(normal);

                normal.normalize();

                vertNormals[get<0>(tri)] += normal;
                vertNormals[get<1>(tri)] += normal;
                vertNormals[get<2>(tri)] += normal;

                triNormals.push_back(normal);
            }

            for(int i = 0; i<vertCount;i++){
                vertNormals[i].normalize();
            }
        }

        tuple<Point, Vector, float> intersect(Point origin, Vector dir){
            float closestDist = numeric_limits<float>::infinity();
            tuple<Point, Vector, float> closestInter = {Point(), Vector(), -1};

            for (int i = 0; i < triangles.size(); i++){
                Vector normal = fullTriNormals[i];
                tuple<int, int, int> triangle = triangles[i];
                Point vert0 = vertices[get<0>(triangle)];
                Point vert1 = vertices[get<1>(triangle)];
                Point vert2 = vertices[get<2>(triangle)];
                float NdotRayDirection = normal.dot(dir);
                float denom = normal.dot(normal);
                float u;
                float v;

                if (fabs(NdotRayDirection) > 0.001){
                     // compute d parameter using equation 2
                    float d = -(normal.x * vert0.x + normal.y * vert0.y + normal.z * vert0.z);
                    
                    // compute t (equation 3)
                    float t = -((normal.x * origin.x + normal.y * origin.y + normal.z * origin.z) + d) / NdotRayDirection;
                    // check if the triangle is behind the ray
                    if (t >= 0 && t < closestDist){
                        
                
                        // compute the intersection point using equation 1
                        Point P = origin + dir * t;
                    
                        // Step 2: inside-outside test
                        Vector C; // vector perpendicular to triangle's plane
                    
                        // edge 0
                        Vector edge0 = vert1 - vert0; 
                        Vector vp0 = P - vert0;
                        C = edge0.cross(vp0);
                        if (normal.dot(C) >= 0){
                            Vector edge1 = vert2 - vert1; 
                            Vector vp1 = P - vert1;
                            C = edge1.cross(vp1);
                            if ((u = normal.dot(C)) >= 0){
                                Vector edge2 = vert0 - vert2; 
                                Vector vp2 = P - vert2;
                                C = edge2.cross(vp2);
                                if ((v = normal.dot(C)) >= 0){
                                    closestDist = t;

                                    Vector vNormal0 = vertNormals[get<0>(triangle)];
                                    Vector vNormal1 = vertNormals[get<1>(triangle)];
                                    Vector vNormal2 = vertNormals[get<2>(triangle)];

                                    u /= denom;
                                    v /= denom;

                                    Vector hitNormal = vNormal2 * (1 - u - v) + vNormal1 * v + vNormal0 * u; 

                                    closestInter = {P, hitNormal.normalized(), t};
                                    //closestInter = {P, normal, t};
                                }
                            }
                        }
                    }
                }
                    
            }
            return closestInter;
        }
};

void swap(float* row1, float* row2){
    for(int i = 0; i < 4; i++){
        float tmp = row1[i];
        row1[i] = row2[i];
        row2[i] = tmp;
    }
}

void partialPivot(float (*m)[4], float (*inverse)[4], int i) {
    double maxElement = abs(m[i][i]);
    int maxRow = i;
    // Encontra a linha com o maior elemento na coluna i
    for (int j = i + 1; j < 4; j++) {
        if (abs(m[j][i]) > maxElement) {
            maxElement = abs(m[j][i]);
            maxRow = j;
        }
    }
    // Troca a linha atual com a linha com o maior elemento
    if (maxRow != i) {
        swap(m[i], m[maxRow]);
        swap(inverse[i], inverse[maxRow]);
    }
}

class Transform{
    public: 
        float matrix[4][4] =
            {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1}
            };

        Point apply(Point& p){
            Point result = Point();
            float res[4] = {0, 0, 0, 0};

            for(int i = 0; i < 4; i++){
                res[i] =
                    this->matrix[i][0]*p.x +
                    this->matrix[i][1]*p.y +
                    this->matrix[i][2]*p.z +
                    this->matrix[i][0]*1;
            }

            result.x = res[0];
            result.y = res[1];
            result.z = res[2];
            return result;
        }

        Vector apply(Vector& v){
            Vector result;
            float res[3] = {0, 0, 0};

            for(int i = 0; i < 3; i++){
                res[i] =
                    this->matrix[i][0]*v.x +
                    this->matrix[i][1]*v.y +
                    this->matrix[i][2]*v.z;
            }

            result.x = res[0];
            result.y = res[1];
            result.z = res[2];
            return result;
        }
        

        Transform copy(){
            Transform result;

            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    result.matrix[i][j] = this->matrix[i][j];
                }
            }

            return result;
        }

        Transform inverse(){
            Transform result;
            Transform copy = this->copy();

            float (*m)[4] = copy.matrix;
            float (*inverse)[4] = result.matrix;

            // Aplica a eliminação gaussiana com pivoteamento parcial
            for (int i = 0; i < 4; i++) {
                partialPivot(m, inverse, i);
                double pivot = m[i][i];
                // Divide a linha atual pelo pivô
                for (int j = 0; j < 4; j++) {
                    m[i][j] /= pivot;
                    inverse[i][j] /= pivot;
                }
                // Subtrai a linha atual das linhas abaixo dela
                for (int j = i + 1; j < 4; j++) {
                    double factor = m[j][i];
                    for (int k = 0; k < 4; k++) {
                        m[j][k] -= factor * m[i][k];
                        inverse[j][k] -= factor * inverse[i][k];
                    }
                }
            }
            // Aplica a eliminação gaussiana com pivoteamento parcial inverso
            for (int i = 4 - 1; i >= 0; i--) {
                for (int j = i - 1; j >= 0; j--) {
                    double factor = m[j][i];
                    for (int k = 0; k < 4; k++) {
                        m[j][k] -= factor * m[i][k];
                        inverse[j][k] -= factor * inverse[i][k];
                    }
                }
            }

            return result;
        }

        Transform operator * (const Transform &other) const{
            Transform result;

            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    result.matrix[i][j] = 
                        this->matrix[i][0]*other.matrix[0][j] +
                        this->matrix[i][1]*other.matrix[1][j] +
                        this->matrix[i][2]*other.matrix[2][j] +
                        this->matrix[i][3]*other.matrix[3][j];
                }
            }

            return result;
        }


        void print(){
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    printf("%f ", this->matrix[i][j]);
                }
                printf("\n");
            }
        }
};


class RotationTransform : public Transform{
    public:
        RotationTransform(float angle, char axis){
            switch (axis)
            {
                case 'x':
                    this->matrix[0][0] =  1;
                    this->matrix[1][1] =  cos(angle);
                    this->matrix[2][1] =  sin(angle);
                    this->matrix[1][2] = -sin(angle);
                    this->matrix[2][2] =  cos(angle);
                    break;
                case 'y':
                    this->matrix[0][0] =  cos(angle);
                    this->matrix[1][1] =  1;
                    this->matrix[0][2] =  sin(angle);
                    this->matrix[2][0] = -sin(angle);
                    this->matrix[2][2] =  cos(angle);
                    break;
                case 'z':
                    this->matrix[0][0] =  cos(angle);
                    this->matrix[1][0] =  sin(angle);
                    this->matrix[0][1] = -sin(angle);
                    this->matrix[1][1] =  cos(angle);
                    this->matrix[2][2] =  1;
                    break;
                default:
                    printf("Eixo de rotação inválido passado, rotação será substituida por transformação identidade\n");
                    break;
                }
        }        
};

