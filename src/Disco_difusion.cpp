#include <iostream>
#include <cmath>
#include "vector.h"
#include <random> //considero que random es mas robusta y estable que random 64 de JD
using namespace std;

//* ------------------------------ Constantes globales ------------------------------
// Constantes Temporales
const double dt = 1e-3;
const double T = 25;

// Constantes Globales
const double Gamma = 16.0;
const double alpha = 1 - exp(-Gamma * dt);
const int N = 10;

// Constantes para fuerzas
const double Cd = 0.47; // Coeficiente de arrastre (esfera)
const double rho = 1.225; // Densidad del aire (kg/m^3)
const double S = 1.0; // Constante de Magnus

//* ------------------------------ Declaraion de clases, metodos de clases y funciones globales ------------------------------
// Declaracion de clases

class Cuerpo;

// Implementacion de las clases

class Cuerpo
{
private:
    vector3D r, v, F;
    double m, R;

public:
    void Inicie(double x0, double y0, double z0, double vx0, double vy0, double vz0, double m0, double R0); // Inicializa las variables de una particula
    void CalculateFuerzas(void); // Calcula las fuerzas que actuan sobre solo una particula (no calcula las fuerzas entre particulas)
    void Muevase(double dt, double KT, mt19937 &gen, normal_distribution<> &dis); // Mueve una particula en un intervalo de tiempo dt
    void Arranque(double dt); // funcion para arrancar el algoritmo de leapfrog
    void Dibujese(void); // Dibuja la particula para gnuplot
    double get_x() { return r.x(); }
    double get_y() { return r.y(); }
    double get_z() { return r.z(); }
    double get_vx() { return v.x(); }
    double get_vy() { return v.y(); }
    double get_vz() { return v.z(); }
    double get_R() { return R; }
};

// Declaracion de funciones globales

void InicieAnimacion(const std::string& str); // Inicia la animacion en gnuplot
void InicieCuadro(void); // Inicia un cuadro en gnuplot
void TermineCuadro(void); // Termina un cuadro en gnuplot
bool stringToBool(const std::string& str); // Convierte un string a un booleano
void Imprimase(double t, Cuerpo *Polen, bool gnuplot); // Imprime los datos de las particulas en un tiempo t tanto para gnuplot como datos crudos

//* ------------------------------ main ------------------------------

int main(int argc, char *argv[]) 
{
    std::string input = argv[1];
    bool gnu = stringToBool(input);

    double r0 = 0.0;
    double KT = 4.0;
    double v0 = r0;
    double m0 = 1.0;
    double R0 = 1.0;

    double Ncuadros = 5000; //TODO: Revisar impresion, a 5000 se etsan generadndo menos datos (4167 en realidad)
    double tcuadro = T / Ncuadros;

    random_device rd;
    mt19937 rng(rd()); // Generador con semilla aleatoria
    // mt19937 rng(1); // Generador con semilla fija
    normal_distribution<> dis(0.0, 1.0); // Distribucion normal con media 0 y desviacion estandar 1

    Cuerpo Polen[N]; // Inicializa un arreglo de particulas

    //Arranque del sistema

    for (int i = 0; i < N; i++) Polen[i].Inicie(r0, 0.0, 0.0, v0, 0.0, 0.0, m0, R0);

    for (int i = 0; i < N; i++) Polen[i].CalculateFuerzas();

    for (int i = 0; i < N; i++) Polen[i].Arranque(dt);

    // Simulacion

    for (double t = 0.0, tdibujo = 0.0; t < T + dt; t += dt, tdibujo += dt)
    {
        if (tdibujo > tcuadro || t == 0)
        {
            Imprimase(t, Polen, gnu);
            tdibujo = 0;
        }

        for (int i = 0; i < N; i++) Polen[i].CalculateFuerzas();
        
        for (int i = 0; i < N; i++) Polen[i].Muevase(dt, KT, rng, dis);
        
    }
    return 0;
}

//* ------------------------------ Implementacion de metodos de clases y funciones globales ------------------------------
// implementacion de los metodos de clases

void Cuerpo::Inicie(double x0, double y0, double z0, double vx0, double vy0, double vz0, double m0, double R0)
{
    // inicializa vectores
    r.load(x0, y0, z0);
    v.load(vx0, vy0, vz0);
    // inicializa variables
    m = m0;
    R = R0;
}

void Cuerpo::CalculateFuerzas(void)
{
    double A = M_PI * R * R;
    // borrar fuerza 
    F.load(0, 0, 0);
    // No hay fuerzas externas a calcular xd
}

void Cuerpo::Muevase(double dt, double KT, mt19937 &gen, normal_distribution<> &dis)
{
    // algoritmo de leapfrog estocastico
    vector3D vp = v + F * dt;
    vector3D epsilon;
    epsilon.load(dis(gen), dis(gen), 0);
    vector3D deltav = -alpha * vp + sqrt(alpha * (2 - alpha) * (KT / m)) * epsilon;
    r += (vp + deltav * 0.5) * dt;
    v = vp + deltav;
}

void Cuerpo::Arranque(double dt)
{
    v -= F * dt / (2 * m);
}

void Cuerpo::Dibujese(void) 
{
    cout << " , " << r.x() << "+" << R << "*cos(t)," << r.y() << "+" << R << "*sin(t)";
}

// implementacion de funciones globales

void InicieAnimacion(const std::string& str) 
{
    cout << "set terminal gif animate" << endl;
    cout << "set output '"<< str <<"'" << endl;
    cout << "unset key" << endl;
    cout << "set xrange[-5:5]" << endl;
    cout << "set yrange[-5:5]" << endl;
    cout << "set size ratio -1" << endl;
    cout << "set parametric" << endl;
    cout << "set trange [0:7]" << endl;
    cout << "set isosamples 12" << endl;
}

void InicieCuadro(void) 
{
    cout << "plot 0,0 ";
}

void TermineCuadro(void) 
{
    cout << endl;
}

bool stringToBool(const std::string& str) 
{
    return (str == "true" || str == "1");
}

void Imprimase(double t, Cuerpo *Polen, bool gnuplot) 
{
    ////cout << gnuplot << endl;
    if (gnuplot) 
    {
        if (t == 0) InicieAnimacion("dif.gif"); //! Cambiar nombre del gif cuando se vaya a ejecutar
        InicieCuadro();
        for (int i = 0; i < N; i++) Polen[i].Dibujese();
        TermineCuadro();
    }
    else
    {
        for (int i = 0; i < N; i++) 
        {
            cout << t << " " << i << " " << Polen[i].get_x() << " " << Polen[i].get_y() << " " << Polen[i].get_vx() << " " << Polen[i].get_vy() << endl;
        }
    }
}
