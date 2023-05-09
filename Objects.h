#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include "Geometry.h"

#define EPSILON 0.001f
using namespace std;

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
            boundingSphere.apply(t);

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
            //acha media dos vertices
            for(int i = 0; i<vertCount;i++){
                avrgVert.x = avrgVert.x + vertices[i].x;
                avrgVert.y = avrgVert.y + vertices[i].y;
                avrgVert.z = avrgVert.z + vertices[i].z;
            }
            boundingSphere.center = avrgVert / vertCount;

            float furthestDist = 0;
            //acha vertice mais distance da media
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