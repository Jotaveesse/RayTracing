#include <fstream>
#include <cmath>

#define PI 3.14159265f

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