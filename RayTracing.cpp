#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <sstream>
#include <tuple>
#include "Objects.h"

#define PI 3.14159265

#define kEpsilon 0.001f
#define MAXBOUNCE 5
#define ENV_REFINDEX 1
#define OBJ_REFINDEX 1.5

using namespace std;

Color intersectRay(Scene scn, vector<Object*> objects, Vector dir, Point origin, int count, bool insideObj);

Color phong( vector<Object*> objects, Scene scn, Object& obj, Point interPoint, Point specPoint, Vector normal, int count, bool insideObj){
    Color finalColor = scn.ambient * obj.ambCo;
    Vector V = (specPoint - interPoint).normalized();
    Color refColor;
    Color tranColor;

    count -= 1;
    
    //recursão de reflexão
    if(obj.refCo != 0 && count>=0){
        Vector reflection = ((normal * 2 * V.dot(normal)) - V).normalized();
        Color Ir = intersectRay(scn, objects, reflection, interPoint, count, insideObj);
        refColor = Ir * obj.refCo;
        refColor.clamp();
    }

    //recursão de refração
    if(obj.tranCo != 0 && count>=0){
        float cos = V.dot(normal);
        Vector projVOnNormal = normal*cos;
        Vector sinVector = V - projVOnNormal;
        float sinRay = sinVector.length();
        Vector normalizedSinVector = sinVector.normalized();

        float sinRef = ENV_REFINDEX*sinRay/OBJ_REFINDEX;
        
        if(insideObj)
            sinRef = OBJ_REFINDEX*sinRay/ENV_REFINDEX;
        
        if(abs(sinRef) > 1)
            sinRef = sinRef > 0 ? 1 : -1;

        float cosRef = 1 - (sinRef*sinRef);
        
        Vector refraction = Vector(0, 0, 0) - (normal*cosRef) - (normalizedSinVector*sinRef);

        Color It = intersectRay(scn, objects, refraction, interPoint, count, 
            (obj.hasInterior() && !insideObj)
        );
        tranColor = It * obj.tranCo;
        tranColor.clamp();
    }

    for(Light light : scn.lights){
        Color diffColor;
        Color specColor;
        Vector Li = (light.center - interPoint);
        float distLight = Li.length();
        Li.normalize();

        Vector Ri = (normal * 2 * Li.dot(normal)) - Li;
        
        bool blocked = false;

        //checa se a luz esta bloqueada por algum objeto
        vector<Object*>::iterator iter;
        for(iter = objects.begin(); iter != objects.end(); iter++) 
        {
            break;
            //ponto de interseção
            tuple<Point, Vector, float> inter = (*iter)->intersect(interPoint, Li, false);
            
            float dist = get<2>(inter);

            //se dist > distLight ponto esta atras da luz
            if(dist >= kEpsilon && dist + kEpsilon< distLight){
                blocked = true;
                break;
            }
        }

        if(!blocked){
            float RiDotV = Ri.dot(V);
            //impede que RiDotV seja negativo
            if(RiDotV < 0)
                RiDotV = 0;

            //phong
            diffColor = light.color * obj.color * obj.difCo * normal.dot(Li);
            diffColor.clamp();

            specColor = light.color * obj.espCo * pow(RiDotV, obj.rugCo);
            specColor.clamp();
            
            finalColor = finalColor + diffColor + specColor;
            finalColor.clamp();
        }
    }

    finalColor = finalColor + refColor + tranColor;
    finalColor.clamp();
    
    return finalColor;
}

Color intersectRay(Scene scn, vector<Object*> objects, Vector dir, Point origin, int count, bool insideObj){
    float closestDist = numeric_limits<float>::infinity();
    tuple<Point, Vector, float> closestInter;
    Object* closestObj = NULL;

    //itera sobre todos os objetos
    vector<Object*>::iterator iter;
    for(iter = objects.begin(); iter != objects.end(); iter++) 
    {
        //ponto de interseção
        tuple<Point, Vector, float> inter = (*iter)->intersect(origin, dir);
        
        float dist = get<2>(inter);

        //se dist < que tamanho do pixVector ponto está entre tela e foco
        if(dist >= kEpsilon && dist < closestDist){
            closestDist = dist;
            closestObj = *iter;
            closestInter = inter;
        }
        
    }
    
    Color finalColor;

    if(closestObj != NULL)
        finalColor = phong(objects, scn, *closestObj, get<0>(closestInter), origin, get<1>(closestInter), count, insideObj);

    return finalColor;
}


void trace(Camera cam, Scene scn, vector<Object*> objects, string fileName){
    fileName = fileName + ".ppm";

    ofstream imagePpm(fileName);

    imagePpm << "P3\n"<< cam.width << " " << cam.height << "\n255\n";

    float fov = 90;

    float gx = cam.distScreen * tan(fov * PI / 180.0 / 2.0);
    float gy = gx * cam.height / cam.width;

    Vector pixWidth = (cam.orthoV * 2 * gx)/(cam.width-1);
    Vector pixHeight = (cam.orthoW * 2 * gy)/(cam.height-1);

    //pixel no canto superior esquerdo
    Vector firstPix = cam.orthoU * cam.distScreen - cam.orthoV * gx + cam.orthoW * gy;

    for(int i=0; i<cam.height;i++){
        for(int j=0; j<cam.width;j++){
            
            //vector que vai do foco pro pixel
            Vector pixVector = firstPix + pixWidth * (j-1) - pixHeight * (i-1);

            Color finalColor = intersectRay(scn, objects, pixVector, cam.center, MAXBOUNCE, false);

            imagePpm << (int)(finalColor.R*255) << " ";
            imagePpm << (int)(finalColor.G*255) << " ";
            imagePpm << (int)(finalColor.B*255) << " ";
        }
        imagePpm << "\n";
    }

    imagePpm.close();

}

Sphere* extractSphere(vector<string> valueArr){
    Point center(stof(valueArr[1]),
        stof(valueArr[2]),
        stof(valueArr[3]));
    
    float radius = stof(valueArr[4]);

    Color col(stof(valueArr[5]),
        stof(valueArr[6]),
        stof(valueArr[7]));
    col.normalize();

    float difCo = stof(valueArr[8]);
    float espCo = stof(valueArr[9]);
    float ambCo = stof(valueArr[10]);
    float refCo = stof(valueArr[11]);
    float tranCo = stof(valueArr[12]);
    float rugCo = stof(valueArr[13]);

    Sphere* sph= new Sphere(center, radius, col,
        difCo, espCo, ambCo, refCo, tranCo, rugCo);

    return sph;
}

Plane* extractPlane(vector<string> valueArr){
    Point p(
        stof(valueArr[1]),
        stof(valueArr[2]),
        stof(valueArr[3]));
    
    Vector normal(
        stof(valueArr[4]),
        stof(valueArr[5]),
        stof(valueArr[6]));

    Color col(
        stof(valueArr[7]),
        stof(valueArr[8]),
        stof(valueArr[9]));
    col.normalize();

    float difCo = stof(valueArr[10]);
    float espCo = stof(valueArr[11]);
    float ambCo = stof(valueArr[12]);
    float refCo = stof(valueArr[13]);
    float tranCo = stof(valueArr[14]);
    float rugCo = stof(valueArr[15]);

    Plane* pln = new Plane(p, normal, col,
        difCo, espCo, ambCo, refCo, tranCo, rugCo);
    
    return pln;
}

Camera* extractCamera(vector<string> valueArr){
    int width = stoi(valueArr[1]);
    int height = stoi(valueArr[2]);
    
    float distScreen = stof(valueArr[3]);
    Vector up(stof(valueArr[4]),
        stof(valueArr[5]),
        stof(valueArr[6])); 
    
    Point center(stof(valueArr[7]),
        stof(valueArr[8]),
        stof(valueArr[9]));
    
    Point focus(stof(valueArr[10]),
        stof(valueArr[11]),
        stof(valueArr[12]));

    return new Camera(height, width, distScreen, up, center, focus);
}

Light extractLight(vector<string> valueArr){
    Point center(stof(valueArr[1]),
        stof(valueArr[2]),
        stof(valueArr[3]));
    
    Color intensity(stof(valueArr[4]),
        stof(valueArr[5]),
        stof(valueArr[6]));
    intensity.normalize();

    Light light(center, intensity);

    return light;
}

int main() {
    vector<Light> lightList;
    vector<Object*> objectList;

    Camera *globalCam;
    Scene *globalScene;

    string value;
    vector<string> valueArr;
    string line;

    ifstream inputFile("input.txt");

    //itera o arquivo, linha pro linha
    while (getline(inputFile, line)) {
        valueArr.clear();
        stringstream streamLine(line);

        //coloca cada valor da linha em uma array
        while (getline(streamLine, value, ' ')) {
            valueArr.push_back(value);
        }

        //esfera
        if(valueArr[0] == "s"){
            objectList.push_back(extractSphere(valueArr));
        }
        //plano
        else if(valueArr[0] == "p"){
            objectList.push_back(extractPlane(valueArr));
        }
        //malha
        else if(valueArr[0] == "t"){
            int triCount = stoi(valueArr[1]);
            int vertCount = stoi(valueArr[2]);

            vector<Point> vertices;
            vector<tuple<int, int, int>> triangles;

            //extrai os vetores de cada linha
            for (int i=0;i<vertCount;i++){
                valueArr.clear();

                getline(inputFile, line);
                stringstream streamLine(line);

                //coloca cada valor da linha em uma array
                while (getline(streamLine, value, ' ')) {
                    valueArr.push_back(value);
                }

                Point vertPoint(
                    stof(valueArr[0]),
                    stof(valueArr[1]),
                    stof(valueArr[2]));

                vertices.push_back(vertPoint);
            }

            //extrai os triangulos de cada linha
            for (int i=0;i<triCount;i++){
                valueArr.clear();

                getline(inputFile, line);
                stringstream streamLine(line);

                //coloca cada valor da linha em uma array
                while (getline(streamLine, value, ' ')) {
                    valueArr.push_back(value);
                }

                tuple<int, int, int> vertIndex{
                    stof(valueArr[0])-1,
                    stof(valueArr[1])-1,
                    stof(valueArr[2])-1};

                triangles.push_back(vertIndex);
            }

            valueArr.clear();

            getline(inputFile, line);
            stringstream streamLine(line);

            //coloca cada valor da linha em uma array
            while (getline(streamLine, value, ' ')) {
                valueArr.push_back(value);
            }

            Color col(
                stof(valueArr[0]),
                stof(valueArr[1]),
                stof(valueArr[2]));
            col.normalize();

            float difCo = stof(valueArr[3]);
            float espCo = stof(valueArr[4]);
            float ambCo = stof(valueArr[5]);
            float refCo = stof(valueArr[6]);
            float tranCo = stof(valueArr[7]);
            float rugCo = stof(valueArr[8]);

            Mesh* mesh = new Mesh(triCount, vertCount, vertices, triangles, col,
                difCo, espCo, ambCo, refCo, tranCo, rugCo);

            objectList.push_back(mesh);
        }
        //camera
        else if(valueArr[0] == "c"){
            globalCam = extractCamera(valueArr);
        }
        //luzes
        else if(valueArr[0] == "l"){
            lightList.push_back(extractLight(valueArr));
        }
        //ambiente
        else if(valueArr[0] == "a"){
            Color col(stof(valueArr[1]),
                stof(valueArr[2]),
                stof(valueArr[3]));
            col.normalize();

            globalScene = new Scene(col, lightList);
        }
    }

    inputFile.close();

    Vector translate(-20, 0, 0);
    TranslationTransform t(translate);
    RotationTransform rt = RotationTransform(33 , 'z');

    //rotate cam around 10 0 0
    Point rotCenter(80, 4, 0);
    TranslationTransform translCenter(rotCenter);
    Transform rotation = t*translCenter.inverse()*rt*translCenter;

    trace(*globalCam, *globalScene, objectList, "image0");
    objectList[0]->apply(t);
    trace(*globalCam, *globalScene, objectList, "image1");
    globalCam->apply(rotation);
    trace(*globalCam, *globalScene, objectList, "image2");
    
    return 0;
};