#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <sstream>
#include <tuple>
#include <omp.h>
#include "Objects.h"
#include <bits/stdc++.h>
#define PI 3.14159265f
#define EPSILON 0.001f
#define MAXBOUNCE 7
#define ENV_REFINDEX 1
#define FOV 90
#define SHADOWS true
#define REFLECTION true
#define REFRACTION true

using namespace std;

Color intersectRay(Scene& scn, vector<Object*>& objects, Vector& dir, Point& origin, int bounceCount);

Color phong(vector<Object*>& objects, Scene& scn, Object& obj, Point& interPoint, Point& specPoint, Vector& normal, int bounceCount){
    Color finalColor = scn.ambient * obj.ambCo;
    Vector interSpectVec = (specPoint - interPoint).normalized();
    Color reflColor;
    Color tranColor;

    bounceCount -= 1;
    
    //recursão de reflexão
    if(REFLECTION){
        if(obj.reflCo != 0 && bounceCount >= 0){
            Vector reflection = interSpectVec.reflect(normal);
            Color Ir = intersectRay(scn, objects, reflection, interPoint, bounceCount);
            reflColor = Ir * obj.reflCo;
        }
    }

    //recursão de refração
    if(REFRACTION){
        if(obj.tranCo != 0 && bounceCount >= 0){
            Vector refraction = interSpectVec.refract(normal, obj.refrInd);

            Color It = intersectRay(scn, objects, refraction, interPoint, bounceCount);
            tranColor = It * obj.tranCo;
        }
    }

    for(Light light : scn.lights){
        Color partialColor;
        Color diffColor;
        Color specColor;
        Vector Li = (light.center - interPoint);
        float distLight = Li.length();
        float lightPassed = 1;

        Li.normalize();

        Vector Ri = (normal * 2 * Li.dot(normal)) - Li;
        
        //checa se a luz esta bloqueada por algum objeto
        if(SHADOWS){
            vector<Object*>::iterator iter;
            for(iter = objects.begin(); iter != objects.end(); iter++) 
            {
                //ponto de interseção
                tuple<Point, Vector, float> inter = (*iter)->intersect(interPoint, Li);
                
                float dist = get<2>(inter);

                //se dist > distLight ponto esta atras da luz
                if(dist >= EPSILON && dist + EPSILON < distLight){
                    lightPassed *= (*iter)->tranCo;

                    //se a luz está totalmente bloqueada
                    if(lightPassed <= EPSILON)
                        break;
                }
            }
        }

        if(lightPassed > EPSILON){
            float RiDotV = Ri.dot(interSpectVec);
            //impede que RiDotV seja negativo
            if(RiDotV < 0)
                RiDotV = 0;

            //phong
            diffColor = light.color * obj.color * obj.difCo * normal.dot(Li);
            diffColor.clamp();

            specColor = light.color * obj.espCo * pow(RiDotV, obj.rugCo);
            specColor.clamp();
            
            partialColor = diffColor + specColor;
            partialColor.clamp();

            partialColor = partialColor * lightPassed;
        }

        finalColor = finalColor + partialColor;
    }

    finalColor = finalColor + reflColor + tranColor;
    finalColor.clamp();
    
    return finalColor;
}

Color intersectRay(Scene& scn, vector<Object*>& objects, Vector& dir, Point& origin, int bounceCount){
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
        if(dist >= EPSILON && dist < closestDist){
            closestDist = dist;
            closestObj = *iter;
            closestInter = inter;
        }
        
    }
    
    Color finalColor;

    if(closestObj != NULL)
        finalColor = phong(objects, scn, *closestObj, get<0>(closestInter), origin, get<1>(closestInter), bounceCount);

    return finalColor;
}


void trace(Camera& cam, Scene& scn, vector<Object*>& objects, string fileName){
    vector<vector<vector<int>>> pixels(cam.height, vector<vector<int>>(cam.width, vector<int>(3, 0)));

    float gx = cam.distScreen * tan(FOV * PI / 180.0 / 2.0);
    float gy = gx * cam.height / cam.width;

    Vector pixWidth = (cam.orthoV * 2 * gx)/(cam.width-1);
    Vector pixHeight = (cam.orthoW * 2 * gy)/(cam.height-1);

    //pixel no canto superior esquerdo
    Vector firstPix = cam.orthoU * cam.distScreen - cam.orthoV * gx + cam.orthoW * gy;

    #pragma omp parallel for
    for(int h=0; h<cam.height;h++){
        for(int w=0; w<cam.width;w++){
            //vector que vai do foco pro pixel
            Vector pixVector = firstPix + pixWidth * (w-1) - pixHeight * (h-1);

            Color finalColor = intersectRay(scn, objects, pixVector, cam.center, MAXBOUNCE);

            pixels[h][w][0] = (int)(finalColor.R*255);
            pixels[h][w][1] = (int)(finalColor.G*255);
            pixels[h][w][2] = (int)(finalColor.B*255);
        }
    }

    fileName = fileName + ".ppm";
    ofstream imagePpm(fileName);

    imagePpm << "P3\n"<< cam.width << " " << cam.height << "\n255\n";

    for(int h=0; h<cam.height;h++){
        for(int w=0; w<cam.width;w++){
            imagePpm << pixels[h][w][0] << " ";
            imagePpm << pixels[h][w][1] << " ";
            imagePpm << pixels[h][w][2] << " ";
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
    float reflCo = stof(valueArr[11]);
    float tranCo = stof(valueArr[12]);
    float rugCo = stof(valueArr[13]);
    float refrInd = stof(valueArr[14]);

    Sphere* sph= new Sphere(center, radius, col,
        difCo, espCo, ambCo, reflCo, tranCo, rugCo, refrInd);

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
    float reflCo = stof(valueArr[13]);
    float tranCo = stof(valueArr[14]);
    float rugCo = stof(valueArr[15]);
    float refrInd = stof(valueArr[16]);

    Plane* pln = new Plane(p, normal, col,
        difCo, espCo, ambCo, reflCo, tranCo, rugCo, refrInd);
    
    return pln;
}

Paraboloid* extractParaboloid(vector<string> valueArr){
    Point focus(
        stof(valueArr[1]),
        stof(valueArr[2]),
        stof(valueArr[3]));

    Point planePoint(
        stof(valueArr[4]),
        stof(valueArr[5]),
        stof(valueArr[6]));
    
    Vector normal(
        stof(valueArr[7]),
        stof(valueArr[8]),
        stof(valueArr[9]));

    Color col(
        stof(valueArr[10]),
        stof(valueArr[11]),
        stof(valueArr[12]));
    col.normalize();

    float difCo = stof(valueArr[13]);
    float espCo = stof(valueArr[14]);
    float ambCo = stof(valueArr[15]);
    float reflCo = stof(valueArr[16]);
    float tranCo = stof(valueArr[17]);
    float rugCo = stof(valueArr[18]);
    float refrInd = stof(valueArr[19]);

    Paraboloid* prbl = new Paraboloid(focus, planePoint, normal, col,
        difCo, espCo, ambCo, reflCo, tranCo, rugCo, refrInd);
    
    return prbl;
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

    //itera o arquivo, linha por linha
    while (getline(inputFile, line)) {
        valueArr.clear();
        stringstream streamLine(line);

        //coloca cada valor da linha em uma array
        while (getline(streamLine, value, ' ')) {
            valueArr.push_back(value);
        }
        
        //caso seja uma linha vazia
        if(valueArr.size() == 0)
            continue;

        //esfera
        if(valueArr[0] == "s"){
            objectList.push_back(extractSphere(valueArr));
        }
        //plano
        else if(valueArr[0] == "p"){
            objectList.push_back(extractPlane(valueArr));
        }
        //paraboloide
        else if(valueArr[0] == "b"){
            objectList.push_back(extractParaboloid(valueArr));
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
            float reflCo = stof(valueArr[6]);
            float tranCo = stof(valueArr[7]);
            float rugCo = stof(valueArr[8]);
            float refrInd = stof(valueArr[9]);

            Mesh* mesh = new Mesh(triCount, vertCount, vertices, triangles, col,
                difCo, espCo, ambCo, reflCo, tranCo, rugCo, refrInd);

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

    Vector translate(-10, 0, -10);
    TranslationTransform t(translate);
    RotationTransform rt = RotationTransform(-50 , 'y');

    clock_t start, end;
    start = clock();
 
    trace(*globalCam, *globalScene, objectList, "image0");
    globalCam->apply(rt);
    trace(*globalCam, *globalScene, objectList, "image1");
    objectList[0]->apply(t);
    trace(*globalCam, *globalScene, objectList, "image2");
    
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << fixed
         << time_taken << setprecision(5);
    cout << " sec " << endl;
    return 0;
};