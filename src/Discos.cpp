#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//* ------------------------------ Constantes globales ------------------------------
// Constantes Temporales
const double dt = 6e-4;
const double tf = 100;

// Constantes Globales
const double Gamma = 16.0;
const double alpha = 1 - exp(-Gamma * dt);

const int Nx=5, Ny=5; //Número de Particulas (cuadrícula)
const int N=Nx*Ny;
const double Lx=Nx*10,Ly=Ny*10, dx=Lx/Nx, dy=Ly/Ny;

// Constantes para fuerzas
const double Cd = 0.47; // Coeficiente de arrastre (esfera)
const double rho = 1.225; // Densidad del aire (kg/m^3)
const double S = 1.0; // Constante de Magnus
const double KH=1e4; //Constante de Hertz

//* ------------------------------ Declaraion de clases, metodos de clases y funciones globales ------------------------------

//Declarar clases----------------
class Cuerpo;
class Colisionador;

//----------------Declaraion de la clase Cuerpo--------------
class Cuerpo{
    private:
     vector3D r, V, w, F;   //Vectores posición, velocidad, velocidad ángular y fuerza
     double m, R, omega; // Masa, radio y velocidad ángular
     double A =2*R; // Longitud característica del disco

    public:
        void Inicie(double x0, double y0, double vx0, double vy0, double omega0, double m0, double R0);
        void BorreFuerza(void){F.load(0,0,0);};
        void SumeFuerza(vector3D dF){F+=dF;};
        void CalculeFuerza(void);
        void Muevase(double dt, double kT, Crandom & ran64);
        void Arranque(double dt);
        double Get_x(void){return r.x();};
        double Get_y(void){return r.y();};
        double Get_z(void){return r.z();};
        double Get_Vx(void){return V.x();};
        double Get_Vy(void){return V.y();};
        double Get_Vz(void){return V.z();};
        void Dibujese(void);
        friend class Colisionador;           
};

//------------  Declaraion de la clase Colisionador  ----------------
class Colisionador{
    private:
        
    public:
        void CalculeAllF(Cuerpo*Particulas);
        void CalculeF_colision(Cuerpo & Particula1,Cuerpo & Particula2);
        void CalculeF_pared(Cuerpo & Particula1,Cuerpo & pared);
        void CalculeF_magnus(Cuerpo & Particula);
};

// Declaracion de funciones globales

void InicieAnimacion(const std::string& str); // Inicia la animacion en gnuplot
void InicieCuadro(void); // Inicia un cuadro en gnuplot
void TermineCuadro(void); // Termina un cuadro en gnuplot
bool stringToBool(const std::string& str); // Convierte un string a un booleano
void Imprimase(double t, Cuerpo *Polen, bool gnuplot); // Imprime los datos de las particulas en un tiempo t tanto para gnuplot como datos crudos

//----------------------------------   Programa principal   -------------------------------

int main(int argc, char *argv[]){
    Crandom ran64(1);  //generador números aleatorios
    Colisionador Col;

    //lector 
    std::string input = argv[1];
    bool gnu = stringToBool(input);

    //Constantes iniciales
    double KT = 4.0;
    double m0 = 1.0;
    double R0 = 1.0;
    double omega0 = 100.0;

    //Variable auxiliares para la animación
    int i, Ncuadros=500; double t,tdibujo=0,tcuadro=tf/Ncuadros; 

    //Iniciar pared
    Cuerpo Particula[N+1];

    double Rpared=50*R0, Mpared=10*R0;
    Particula[N].Inicie(0,0,0,0,0,Mpared,Rpared); //pared circular

    //Iniciar los cuerpos
    double theta, x0, y0;
    for(int ix=0;ix<Nx;ix++){
        for(int iy=0;iy<Ny;iy++){
            theta=2*M_PI*ran64.r();
            x0=(ix-2)*dx; y0=(iy-2)*dy; 
            //------------------------(x0, y0, Vx0, Vy0,  w0,  m0, R0)
            Particula[iy*Nx+ix].Inicie(x0, y0, 0, 0.0, omega0, m0, R0);

        } 
    }

    //for(i=0;i<N;i++) Particula[i].CalculeFuerza();
    Col.CalculeAllF(Particula); 
    for(i=0;i<N;i++) Particula[i].Arranque(dt);

    for(t=0;t<tf;t+=dt,tdibujo+=dt){
        
        if(tdibujo > tcuadro || t == 0)
        {
            Imprimase(t, Particula, gnu);
            tdibujo = 0;
        }
    
        //for(i=0;i<N;i++) Particula[i].CalculeFuerza();
        Col.CalculeAllF(Particula); 
        for(i=0;i<N;i++) Particula[i].Muevase(dt,KT,ran64);
    }
    return 0;
}

//----------  Implementar funciones de la clase Cuerpo  --------------

void Cuerpo::Inicie(double x0, double y0, double vx0, double vy0, double omega0, double m0, double R0){
    // inicializa vectores posición y velocidad
    r.load(x0, y0, 0.0); V.load(vx0, vy0, 0.0);
    // inicializa variables
    m=m0; R=R0; omega=omega0;
}

void Cuerpo::CalculeFuerza(void){
    //Declarar velocidad ángular
    vector3D w;
    w.load(0, 0, omega);

    F.load(0, 0, 0); // borrar fuerza 

    // Agrergar fuerzas

    vector3D Fm = -0.5 * Cd * rho * A * R * (w ^ V); // Fuerza de Magnus

    // Sumando las fuerzas
    F += Fm;
}

void Cuerpo::Muevase(double dt, double kT, Crandom & ran64){
    //Algotimo Browniano leap-frog estocástico 
    vector3D Vprime=V+F*dt;
    vector3D epsilon; epsilon.load(ran64.gauss(0,1), ran64.gauss(0,1), 0); //Variable aleatoria
    vector3D deltaV= -alpha*Vprime+sqrt(alpha*(2-alpha)*kT/m)*epsilon;
    r+=(Vprime+deltaV*0.5)*dt; 
    V=Vprime+deltaV;
}

void Cuerpo::Arranque(double dt){
    //Algotimo leap-frog
    V-=F*(dt/(2*m));
}

void Cuerpo::Dibujese(void){
    cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//----------  Implementar funciones de la clase Colisionador  --------------

void Colisionador::CalculeAllF(Cuerpo*Particulas){
    int i,j;
    //Borrar la fuerza de todas los Particulas
    for(i=0;i<N;i++){
        Particulas[i].BorreFuerza();
    };
    //Calcular fuerza de magus para cada partícula
    for(i=0;i<N;i++){
        CalculeF_magnus(Particulas[i]);
    };
    //Calcular Fuerza entre colisiones de las partículas
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            CalculeF_colision(Particulas[i],Particulas[j]);
        } 
    }
    //Calcular Fuerza entre cada partícula y la pared
    for(i=0;i<N;i++){
        CalculeF_pared(Particulas[i],Particulas[N]);
    }
}

void Colisionador::CalculeF_colision(Cuerpo & Particula1,Cuerpo & Particula2){
    //Determinar si hay colisión
    vector3D r21=Particula2.r-Particula1.r;
    double rn=r21.norm();
    double R1=Particula1.R, R2=Particula2.R; 
    double s=(R1+R2)-rn;
    if(s>0){//Si hay colisión
        //Calcular el vector normal
        vector3D n=r21*(1/rn);
        //Calcular fuerza
        vector3D F=n*(KH*pow(s,1.5));
        Particula2.SumeFuerza(F); Particula1.SumeFuerza(-1*F);
    };
}

void Colisionador::CalculeF_pared(Cuerpo & Particula,Cuerpo & pared){
    //Fuerza elástica de la pared circular
    double rn=Particula.r.norm();
    double Rp=pared.R, R=Particula.R; 
    double s=(R+rn)-Rp;
    if(s>0){//Si hay colisión
        //Calcular el vector normal
        vector3D n=Particula.r*(1/rn);
        //Calcular fuerza
        vector3D F=n*(KH*pow(s,1.5));
        Particula.SumeFuerza(-1*F); 
    };
}
void Colisionador::CalculeF_magnus(Cuerpo & Particula){
    vector3D Fm = -0.5 * Cd * rho * Particula.A * Particula.R * (Particula.w ^ Particula.V); // Fuerza de Magnus
    Particula.SumeFuerza(Fm); 
}
//--------------------------
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
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<1e-6<<"+"<<50<<"*cos(t),"<<1e-6<<"+"<<50<<"*sin(t)"; //Pared Circular
}
void TermineCuadro(void){
    cout<<endl;
}

bool stringToBool(const std::string& str) 
{
    return (str == "true" || str == "1");
}

void Imprimase(double t, Cuerpo *Polen, bool gnuplot) 
{
    if (gnuplot) 
    {
        if (t == 0) InicieAnimacion("fcen.gif"); //! Cambiar nombre del gif cuando se vaya a ejecutar
        InicieCuadro();
        for (int i = 0; i < N; i++) Polen[i].Dibujese();
        TermineCuadro();
    }
    else
    {
        for (int i = 0; i < N; i++) 
        {
            cout << t << " " << i << " " << Polen[i].Get_x() << " " << Polen[i].Get_y() << " " << Polen[i].Get_Vx() << " " << Polen[i].Get_Vy() << endl;
        }
    }
}