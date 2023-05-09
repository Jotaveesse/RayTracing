#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>

#define PI 3.14159265f
#define EPSILON 0.001f
using namespace std;

class Point;

void swap(float* row1, float* row2);
void partialPivot(float (*m)[4], float (*inverse)[4], int i);

class Vector{
    public:
        float x, y, z;

        Vector(float in_x, float in_y, float in_z){
            x = in_x;
            y = in_y;
            z = in_z;
        }
        Vector(float in_x){
            x = in_x;
            y = in_x;
            z = in_x;
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
                x /= len;
                y /= len;
                z /= len;
            }

            return *this;
        }

        Vector normalized(){
            float len = length();
            if(len>0){
                Vector vec(x / len, y / len, z / len);
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

        Vector reflect(Vector normal){
            return ((normal * 2 * normal.dot(*this)) - *this).normalize();
        }

        Vector refract(Vector normal, float relRefrIndex){
            float cos = normal.dot(*this);
            Vector projOnNormal = normal*cos;
            Vector sinVector = *this - projOnNormal;
            float sinRay = sinVector.length();
            Vector normalizedSinVector = sinVector.normalized();

            float sinRef;

            if(cos>0)   //entrando no objeto
                sinRef = sinRay / relRefrIndex;
            else   //saindo saindo do objeto
                sinRef = relRefrIndex * sinRay;
            
            if(abs(sinRef) > 1)
                sinRef = sinRef > 0 ? 1 : -1;

            float cosRef = sqrt(1 - (sinRef*sinRef));
            
            Vector refraction = (normal*cosRef + normalizedSinVector*sinRef) * -1;

            return refraction;
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

        friend std::ostream& operator << (std::ostream &s, const Vector &v)
        {
            return s << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
        }
};

class Point{
    public:
        float x, y, z;

        Point(float in_x, float in_y, float in_z){
            x = in_x;
            y = in_y;
            z = in_z;
        }
        Point(float in_x){
            x = in_x;
            y = in_x;
            z = in_x;
        }
        Point(){
            x = 0;
            y = 0;
            z = 0;
        }

        float dist(Point p){
            return sqrt(pow(x - p.x, 2) + pow(y - p.y, 2) + pow(z - p.z, 2));
        }

        float sqrdDist(Point p){
            return pow(x - p.x, 2) + pow(y - p.y, 2) + pow(z - p.z, 2);
        }

        Point operator + (const Vector &v) const
        {return Point(x + v.x, y + v.y, z + v.z);}

        Point operator - (const Vector &v) const
        {return Point(x - v.x, y - v.y, z - v.z);}

        Vector operator + (const Point &v) const
        {return Vector(x + v.x, y + v.y, z + v.z);}

        Vector operator - (const Point &v) const
        {return Vector(x - v.x, y - v.y, z - v.z);}

        Point operator * (const float &n) const
        {return Point(x * n, y * n, z * n);}

        Point operator / (const float &n) const
        {return Point(x / n, y / n, z / n);}

        float operator [] (uint8_t i) const { return (&x)[i]; }
        float& operator [] (uint8_t i) { return (&x)[i]; }

        friend std::ostream& operator << (std::ostream &s, const Point &v)
        {
            return s << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
        }
};


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
                    this->matrix[i][3]*1;
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
                float pivot = m[i][i];
                // Divide a linha atual pelo pivô
                for (int j = 0; j < 4; j++) {
                    m[i][j] /= pivot;
                    inverse[i][j] /= pivot;
                }
                // Subtrai a linha atual das linhas abaixo dela
                for (int j = i + 1; j < 4; j++) {
                    float factor = m[j][i];
                    for (int k = 0; k < 4; k++) {
                        m[j][k] -= factor * m[i][k];
                        inverse[j][k] -= factor * inverse[i][k];
                    }
                }
            }
            // Aplica a eliminação gaussiana com pivoteamento parcial inverso
            for (int i = 4 - 1; i >= 0; i--) {
                for (int j = i - 1; j >= 0; j--) {
                    float factor = m[j][i];
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

        Camera(int in_height, int in_width, float in_distScreen, Vector in_up,
        Point in_center, Point in_target){
            height = in_height;
            width = in_width;
            distScreen = in_distScreen;
            up = in_up;
            center = in_center;
            target = in_target;

            orthoU = (target - center).normalize();
            orthoV = orthoU.cross(up).normalize();
            orthoW = orthoV.cross(orthoU);
        }

        void apply(Transform& t){
            center = t.apply(center);
            target = t.apply(target);
            up = t.apply(up);
            orthoU = t.apply(orthoU);
            orthoV = t.apply(orthoV);
            orthoW = t.apply(orthoW);
        }
};

class Color{
    public:
        float R, G, B;
        
        Color(float in_R, float in_G, float in_B){
            R = in_R;
            G = in_G;
            B = in_B;
        }

        Color(float C){
            R = C;
            G = C;
            B = C;
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
            else if(R>1)
                R = 1;
            if(G<0)
                G = 0;
            else if(G>1)
                G = 1;
            if(B<0)
                B = 0;
            else if(B>1)
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
        {return Color(R * n, G * n, B * n);}

        Color operator + (const Color &c) const
        {return Color(R + c.R, G + c.G, B + c.B);}

        Color operator * (const Color &c) const
        {return Color(R * c.R, G * c.G, B * c.B);}
};

class Light{
    public:
        Point center;
        Color color;

        Light(Point in_center, Color in_color){
            color = in_color;
            center = in_center;
        }
};

class Scene{
    public:
        Color ambient;
        vector<Light> lights = {};
        Scene(Color in_ambient, vector<Light> in_lights){
            ambient = in_ambient;
            lights = in_lights;
        }
};

class Object{
    public:
        Color color;
        float difCo;
        float espCo;
        float ambCo;
        float reflCo;
        float tranCo;
        float rugCo;
        float refrInd;

        explicit Object(Color in_color, float in_difCo,
        float in_espCo, float in_ambCo, float in_reflCo, float in_tranCo, float in_rugCo, float in_refrInd){
            color = in_color;
            difCo = in_difCo;
            espCo = in_espCo;
            ambCo = in_ambCo;
            reflCo = in_reflCo;
            tranCo = in_tranCo;
            rugCo = in_rugCo;
            refrInd = in_refrInd;
        }

        explicit Object(){
        }

        virtual void apply(Transform& t){
            
        }

        virtual tuple<Point, Vector, float> intersect(Point origin, Vector dir){
            return tuple<Point, Vector, float>{Point(), Vector(), -1};
        }
};

class Sphere: public virtual Object{
    public:
        Point center;
        float radius;
        
        Sphere(Point in_center, float in_radius, Color in_color, float in_difCo,
        float in_espCo, float in_ambCo, float in_reflCo, float in_tranCo, float in_rugCo, float in_refrInd):
        Object(in_color, in_difCo, in_espCo, in_ambCo, in_reflCo, in_tranCo, in_rugCo, in_refrInd){
            center = in_center;
            radius = in_radius;
            color = in_color;
            difCo = in_difCo;
            espCo = in_espCo;
            ambCo = in_ambCo;
            reflCo = in_reflCo;
            tranCo = in_tranCo;
            rugCo = in_rugCo;
            refrInd = in_refrInd;
        }

        Sphere():
        Object(){
        }

        void apply(Transform& t){
            center = t.apply(center);
        }

        tuple<Point, Vector, float> intersect(Point origin, Vector dir){
            Vector L = center - origin;
            float lengthL = L.sqrdLength();
            dir.normalize();
            float tc = L.dot(dir);

            if (tc >= 0){
                float d = sqrt(lengthL - tc*tc);

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

        Plane(Point in_point, Vector in_normal, Color in_color, float in_difCo,
        float in_espCo, float in_ambCo, float in_reflCo, float in_tranCo, float in_rugCo, float in_refrInd):
        Object(in_color, in_difCo, in_espCo, in_ambCo, in_reflCo, in_tranCo, in_rugCo, in_refrInd){
            point = in_point;
            normal = in_normal;
            color = in_color;
            difCo = in_difCo;
            espCo = in_espCo;
            ambCo = in_ambCo;
            reflCo = in_reflCo;
            tranCo = in_tranCo;
            rugCo = in_rugCo;
            refrInd = in_refrInd;
        }

        void apply(Transform& t){
            point = t.apply(point);
            normal = t.apply(normal);
        }

        tuple<Point, Vector, float> intersect(Point origin, Vector dir){
            dir.normalize();
            float denom = normal.dot(dir);

            if(abs(denom) > EPSILON){
                Vector v = point - origin;
                float t = v.dot(normal) / denom;
                if(t >= 0){
                    Point inters = origin + (dir * t);
                    
                    if(normal.dot(dir)<0)
                        return tuple<Point, Vector, float>{inters, normal, t};
                    else
                        return tuple<Point, Vector, float>{inters, normal * -1, t};
                }
            }
            return tuple<Point, Vector, float>{Point(), Vector(), -1};
        }
};

class Paraboloid: public virtual Object{
    public:
        Point focus;
        Point planePoint;
        Vector planeNormal;

        Paraboloid(Point in_focus, Point in_planePoint, Vector in_planeNormal, Color in_color, float in_difCo,
        float in_espCo, float in_ambCo, float in_reflCo, float in_tranCo, float in_rugCo, float in_refrInd):
        Object(in_color, in_difCo, in_espCo, in_ambCo, in_reflCo, in_tranCo, in_rugCo, in_refrInd){
            focus = in_focus;
            planePoint = in_planePoint;
            planeNormal = in_planeNormal.normalize();
            color = in_color;
            difCo = in_difCo;
            espCo = in_espCo;
            ambCo = in_ambCo;
            reflCo = in_reflCo;
            tranCo = in_tranCo;
            rugCo = in_rugCo;
            refrInd = in_refrInd;
        }

        void apply(Transform& t){
            focus = t.apply(focus);
            planePoint = t.apply(planePoint);
            planeNormal = t.apply(planeNormal);
        }

        tuple<Point, Vector, float> intersect(Point origin, Vector dir){
            dir.normalize();

            //definições das variaveis para o calculo
            float d = - (planeNormal.x * planePoint.x + planeNormal.y * planePoint.y + planeNormal.z * planePoint.z);
            float r = planeNormal.sqrdLength();
            float s = planeNormal.x * origin.x + planeNormal.y * origin.y + planeNormal.z * origin.z + d;
            float t = planeNormal.x * dir.x + planeNormal.y * dir.y + planeNormal.z * dir.z;
            float u = origin.x - focus.x;
            float v = origin.y - focus.y;
            float w = origin.z - focus.z;

            //a b e c  da formula de bhaskara
            float bhaskaraA = pow(t , 2) - pow(dir.x , 2) * r - pow(dir.y , 2) * r - pow(dir.z , 2) * r;
            float bhaskaraB = 2 * t * s - 2 * dir.x * u - 2 * dir.y * v - 2 * dir.z * w;
            float bhaskaraC = pow(s , 2) - pow(u , 2) * r - pow(v , 2) * r - pow(w , 2) * r;

            float delta = pow(bhaskaraB , 2) - 4 * bhaskaraA * bhaskaraC;

            //se axiste alguma interseção
            if(delta >= 0){
                float dist1 = (-bhaskaraB + sqrt(delta)) / (2 * bhaskaraA);
                float dist2 = (-bhaskaraB - sqrt(delta)) / (2 * bhaskaraA);
                
                //se alguma das interseções está na frente da camera
                if(dist1 > 0 || dist2 > 0){
                    float closestDist;

                    if (dist1 < dist2)
                        closestDist = dist1;
                    else
                        closestDist = dist2;
                    
                    Point inters = origin + (dir * closestDist);
                    
                    Vector focusInters = inters - focus;
                    Vector intersPlane = planeNormal.normalize() * -focusInters.length();
                    Vector normal = (focusInters + intersPlane).normalize();
                    
                    return tuple<Point, Vector, float>{inters, normal, closestDist};
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
        vector<Vector> vertNormals;

        Sphere boundingSphere;

        Mesh(int in_triCount, int in_vertCount, vector<Point> in_vertices,
        vector<tuple<int, int, int>> in_triangles, Color in_color, float in_difCo, float in_espCo,
        float in_ambCo, float in_reflCo, float in_tranCo, float in_rugCo, float in_refrInd):
        Object(in_color, in_difCo, in_espCo, in_ambCo, in_reflCo, in_tranCo, in_rugCo, in_refrInd){
            triCount = in_triCount;
            vertCount = in_vertCount;
            vertices = in_vertices;
            triangles = in_triangles;
            color = in_color;
            difCo = in_difCo;
            espCo = in_espCo;
            ambCo = in_ambCo;
            reflCo = in_reflCo;
            tranCo = in_tranCo;
            rugCo = in_rugCo;
            refrInd = in_refrInd;

            getNormals();
            getBoundingSphere();
        }

        void apply(Transform& t){
            for(unsigned int i = 0; i < this->vertices.size(); i++){
                this->vertices[i] = t.apply(this->vertices[i]);
            }

            for(unsigned int i = 0; i < this->triNormals.size(); i++){
                this->triNormals[i] = t.apply(this->triNormals[i]);
            }

            for(unsigned int i = 0; i < this->vertNormals.size(); i++){
                this->vertNormals[i] = t.apply(this->vertNormals[i]);
            }
        }

        void getBoundingSphere(){
            Point avrgVert;
            for(int i = 0; i<vertCount;i++){
                avrgVert.x = avrgVert.x + vertices[i].x;
                avrgVert.y = avrgVert.y + vertices[i].y;
                avrgVert.z = avrgVert.z + vertices[i].z;
            }
            boundingSphere.center = avrgVert / vertCount;

            float furthestDist = 0;
            for(int i = 0; i<vertCount;i++){
                float dist = boundingSphere.center.sqrdDist(vertices[i]);
                if(dist>furthestDist)
                    furthestDist = dist;
            }
            boundingSphere.radius = sqrt(furthestDist);
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

        tuple<Point, Vector, float> intersect(Point origin, Vector dir){
            tuple<Point, Vector, float> closestInter = {Point(), Vector(), -1};

            //checa se raio colide com boundingSphere
            float boundInter = get<2>(boundingSphere.intersect(origin, dir));
            float distOrigCenter = origin.dist(boundingSphere.center);
            if(boundInter<0 && distOrigCenter>boundingSphere.radius)
                return closestInter;

            Vector normal;
            tuple<int, int, int> triangle;
            Point vert0;
            Point vert1;
            Point vert2;
            float closestDist = numeric_limits<float>::infinity();
            float NdotRayDir ;
            float denom;
            float u;
            float v;

            for (unsigned int i = 0; i < triangles.size(); i++){
                normal = triNormals[i];
                triangle = triangles[i];
                vert0 = vertices[get<0>(triangle)];
                vert1 = vertices[get<1>(triangle)];
                vert2 = vertices[get<2>(triangle)];
                NdotRayDir = normal.dot(dir);
                denom = normal.dot(normal);

                //checa se o raio é paralelo ao triangulo
                if (fabs(NdotRayDir) > EPSILON){
                    float d = -(normal.x * vert0.x + normal.y * vert0.y + normal.z * vert0.z);
                    
                    float t = -((normal.x * origin.x + normal.y * origin.y + normal.z * origin.z) + d) / NdotRayDir;
                    
                    //se t < -EPSILON o triangulo está atras do raio
                    if (t < -EPSILON*5 || t >= closestDist)
                        continue;

                    Point interPoint = origin + dir * t;
                    
                    //vetor perpendicular ao plano do triangulo
                    Vector interPerp;

                    Vector edge0 = vert1 - vert0; 
                    Vector vp0 = interPoint - vert0;
                    interPerp = edge0.cross(vp0);

                    //checa se interseção esta do lado certo desta aresta
                    if (normal.dot(interPerp) < 0)
                        continue;

                    Vector edge1 = vert2 - vert1; 
                    Vector vp1 = interPoint - vert1;
                    interPerp = edge1.cross(vp1);

                    //checa se interseção esta do lado certo desta aresta
                    if ((u = normal.dot(interPerp)) < 0)
                        continue;

                    Vector edge2 = vert0 - vert2; 
                    Vector vp2 = interPoint - vert2;
                    interPerp = edge2.cross(vp2);

                    //checa se interseção esta do lado certo desta aresta
                    if ((v = normal.dot(interPerp)) < 0)
                        continue;

                    Vector vNormal0 = vertNormals[get<0>(triangle)];
                    Vector vNormal1 = vertNormals[get<1>(triangle)];
                    Vector vNormal2 = vertNormals[get<2>(triangle)];

                    u /= denom;
                    v /= denom;

                    //interpola as normais dos vetores
                    Vector hitNormal = (vNormal0 * u + vNormal1 * v + vNormal2 * (1 - u - v)).normalize(); 

                    float cosNormalDir = dir.dot(normal.normalized());

                    if(cosNormalDir>0)
                        closestInter = {interPoint, hitNormal * -1, t};
                    else
                        closestInter = {interPoint, hitNormal, t};
                    
                    closestDist = t;
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
    float maxElement = abs(m[i][i]);
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


class RotationTransform : public Transform{
    float _angle;
    char _axis;
    
    public:
        RotationTransform(float angle, char axis){
            this->_angle = angle * (PI/180);
            this->_axis = axis;

            switch (axis)
            {
                case 'x':
                    this->matrix[0][0] =  1;
                    this->matrix[1][1] =  cos(_angle);
                    this->matrix[2][1] =  sin(_angle);
                    this->matrix[1][2] = -sin(_angle);
                    this->matrix[2][2] =  cos(_angle);
                    break;
                case 'y':
                    this->matrix[0][0] =  cos(_angle);
                    this->matrix[1][1] =  1;
                    this->matrix[0][2] =  sin(_angle);
                    this->matrix[2][0] = -sin(_angle);
                    this->matrix[2][2] =  cos(_angle);
                    break;
                case 'z':
                    this->matrix[0][0] =  cos(_angle);
                    this->matrix[1][0] =  sin(_angle);
                    this->matrix[0][1] = -sin(_angle);
                    this->matrix[1][1] =  cos(_angle);
                    this->matrix[2][2] =  1;
                    break;
                default:
                    printf("Eixo de rotação inválido passado, rotação será substituida por transformação identidade\n");
                    break;
            }
        }

        RotationTransform inverse(){
            return RotationTransform(-this->_angle, this->_axis);
        }
};

class TranslationTransform : public Transform {
    public:

        /*
            Este construtor inicializa uma
            transformacao na qual os objetos afetados
            serao trasladados na direcao do vetor
            passado como argumento
        */
        TranslationTransform(Vector& translation){
            this->matrix[0][3] = translation.x;
            this->matrix[1][3] = translation.y;
            this->matrix[2][3] = translation.z;
        }

        
        /*
            Este construtor inicializa uma
            transformacao na qual os objetos afetados
            serao trasladados de forma que o ponto
            passado como argumento seja como a "nova
            origem" do sistema de coordenadas
        */
        TranslationTransform(Point& origin){
            this->matrix[0][3] = -origin.x;
            this->matrix[1][3] = -origin.y;
            this->matrix[2][3] = -origin.z;
        }


        TranslationTransform inverse() {
            Vector newTranslation = Vector(
                -this->matrix[0][3],
                -this->matrix[1][3],
                -this->matrix[2][3]
            );

            return TranslationTransform(newTranslation);
        }
};