#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <sstream>
#include <tuple>
#include "Objects.h"

#define PI 3.14159265
using namespace std;

Color phong(Scene scn, Object obj, Point interPoint, Point specPoint, Vector normal){
    Color finalColor = scn.ambient * obj.ambCo;

    Vector V = (specPoint - interPoint).normalized();

    for(Light light : scn.lights){
        Color calcColor;
        Vector Li = (light.center - interPoint).normalized();
        Vector Ri = (normal * 2 * Li.dot(normal)) - Li;

        //phong
        calcColor = light.color * obj.color * obj.difCo * normal.dot(Li) + light.color * obj.espCo * pow(Ri.dot(V), obj.rugCo) ;

        calcColor.clamp();
        
        finalColor = finalColor + calcColor;
    }

    finalColor.clamp();
    
    return finalColor;
}

void trace(Camera cam, Scene scn, vector<Object*> objects){
    ofstream imagePpm("image.ppm");

    imagePpm << "P3\n"<< cam.width << " " << cam.height << "\n255\n";

    Vector t = (cam.target - cam.center).normalized();
    Vector b = cam.up.cross(t).normalized();
    Vector v = t.cross(b).normalized();
    float fov = 90;

    float gx = cam.distScreen * tan(fov * PI / 180.0 / 2.0);
    float gy = gx * cam.height / cam.width;

    Vector pixWidth = (b * 2 * gx)/(cam.width-1);
    Vector pixHeight = (v * 2 * gy)/(cam.height-1);

    //pixel no canto superior esquerdo
    Vector firstPix = t * cam.distScreen - b * gx + v * gy;

    for(int i=0; i<cam.height;i++){
        for(int j=0; j<cam.width;j++){
            float closestDist = numeric_limits<float>::infinity();
            tuple<Point, Vector, float> closestInter;
            Object* closestObj = NULL;

            //itera sobre todos os objetos
            vector<Object*>::iterator iter;
            for(iter = objects.begin(); iter != objects.end(); iter++) 
            {
                //vector que vai do foco pro pixel
                Vector pixVector = firstPix + pixWidth * (j-1) - pixHeight * (i-1);

                //ponto de interseção
                tuple<Point, Vector, float> inter = (*iter)->intersect(cam.center, pixVector);
                
                float dist = get<2>(inter);

                //se dist < que tamanho do pixVector ponto está entre tela e foco
                if(dist >= pixVector.length() && dist < closestDist){
                    closestDist = dist;
                    closestObj = *iter;
                    closestInter = inter;
                }
            }
            
            Color finalColor;

            if(closestObj != NULL)
                finalColor = phong(scn, *closestObj, get<0>(closestInter), cam.center, get<1>(closestInter));

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
        else if(valueArr[0] == "a"){
            Color col(stof(valueArr[1]),
                stof(valueArr[2]),
                stof(valueArr[3]));
            col.normalize();

            globalScene = new Scene(col, lightList);
        }
    }

    inputFile.close();

    trace(*globalCam, *globalScene, objectList);

    return 0;
};