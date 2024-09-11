#include <iostream>
#include <cmath>
#include  <fstream>
#include "vector.h"
#include "Random64.h"
using namespace std;

//* ------------------------------ Constantes globales ------------------------------
// Constantes Temporales
const double dt = 6e-4;
const double tf = 100;

// Constantes Globales
const double Gamma = 10.0;
const double alpha = 1 - exp(-Gamma * dt);

const int Nx=5, Ny=5; //Número de Particulas (cuadrícula)
const int N=Nx*Ny;
const double Lx=Nx*3,Ly=Ny*3, dx=Lx/Nx, dy=Ly/Ny;

// Constantes del sistema
const double KT = 4.0; // "temperatura"
const double Cd = 0.47; // Coeficiente de arrastre (esfera)
const double rho = 1.225; // Densidad del aire (kg/m^3)
const double S = 1.0; // Constante de Magnus
const double KH=1e4; //Constante de Hertz
const double Kcundall = 500; //Constante de Cundall
const double mu = 0.4 ; // Coeficiente de fricción de disco-disco //TODO: Revisar valor, lo dejo como 0.4 por ahora
const double K = 1.0; // Constante de resorte
const double GammaHertz = 50; // Constante de amortiguamiento de Hertz
//* ------------------------------ Declaraion de clases, metodos de clases y funciones globales ------------------------------

//Declarar clases----------------
class Cuerpo;
class Colisionador;

//----------------Declaraion de la clase Cuerpo--------------
class Cuerpo{
    private:
     vector3D r, V, w, F;   //Vectores posición, velocidad, velocidad ángular y fuerza
     double m, R, theta, omega; // Masa, radio y velocidad ángular
     double A, I, tau; // Longitud característica del disco

    public:
        void Inicie(double x0, double y0, double vx0, double vy0, double omega0, double m0, double R0);
        void BorreFuerza(void)
        {
            F.load(0,0,0);
            tau = 0;
        };
        void SumeFuerza(vector3D dF, double dtau)
        {
            F+=dF;
            tau += dtau;
        };
        void Muevase(double dt, double kT, Crandom & ran64);
        void Arranque(double dt);
        double Get_x(void){return r.x();};
        double Get_y(void){return r.y();};
        double Get_z(void){return r.z();};
        double Get_Vx(void){return V.x();};
        double Get_Vy(void){return V.y();};
        double Get_Vz(void){return V.z();};
        double Get_omega(void){return omega;};
        void Dibujese(void);
        friend class Colisionador;           
};

//------------  Declaraion de la clase Colisionador  ----------------
class Colisionador{
    private:
        double xCundall[N+1][N+1], sold[N+1][N+1];
    public:
        void Inicie(void);
        void CalculeAllF(Cuerpo*Particulas);
        void CalculeF_colision(Cuerpo & Particula1,Cuerpo & Particula2, double &xCundall, double &sold, double dt);
        void CalculeF_pared(Cuerpo & Particula1,Cuerpo & pared, double &xCundall, double &sold, double dt);
        void CalculeF_magnus(Cuerpo & Particula);
        void CalculeF_ceentral(Cuerpo & Particula);
        void CalculeT_ind(Cuerpo & Particula);
};

// Declaracion de funciones globales

void InicieAnimacion(const std::string& str); // Inicia la animacion en gnuplot
void InicieCuadro(void); // Inicia un cuadro en gnuplot
void TermineCuadro(void); // Termina un cuadro en gnuplot
int stringToBool(std::string str); // Convierte un string a un booleano
void Imprimase(double t, Cuerpo *Polen, int gnuplot); // Imprime los datos de las particulas en un tiempo t tanto para gnuplot como datos crudos

//----------------------------------   Programa principal   -------------------------------

int main(int argc, char *argv[]){
    Crandom ran64(1);  //generador números aleatorios
    Colisionador Col;

    //lector 
    std::string input = argv[1];
    //std::string input = "0"; //! Solo para pruebas en debuggin
    int gnu = stringToBool(input);

    //Constantes iniciales
    double m0 = 1.0;
    double R0 = 1.0;
    double omega0 = 0.0;

    //Variable auxiliares para la animación
    int i, Ncuadros=500; double t,tdibujo=0,tcuadro=tf/Ncuadros; 

    //Iniciar pared
    Cuerpo Particula[N+1];

    double Rpared=10*R0, Mpared=50*R0;
    Particula[N].Inicie(0,0,0,0,0,Mpared,Rpared); //pared circular

    //Iniciar los cuerpos
    double theta, x0, y0;
    for(int iy=0;iy<Ny;iy++){
        for(int ix=0;ix<Nx;ix++){
            theta=2*M_PI*ran64.r();
            x0=(ix-2)*dx; y0=(iy-2)*dy; 
            //------------------------(x0, y0, Vx0, Vy0,  w0,  m0, R0)
            Particula[iy*Nx+ix].Inicie(x0, y0, 0.0, 0.0, omega0, m0, R0);

        } 
    }

    Col.CalculeAllF(Particula); 
    for(i=0;i<N;i++) Particula[i].Arranque(dt);

    for(t=0;t<tf;t+=dt,tdibujo+=dt){
        
        if(tdibujo > tcuadro || t == 0)
        {
            Imprimase(t, Particula, gnu);
            tdibujo = 0;
        }
    
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
    m=m0; R=R0; 
    theta=0; omega = omega0;
    A = 2*R;
    I = 0.5 * m * R * R; 
}

void Cuerpo::Muevase(double dt, double kT, Crandom & ran64){
    //Algotimo Browniano leap-frog estocástico 
    vector3D Vprime=V+F*dt/m; //? antes no se devidia por la masa, en clase aparecia asi, aun asi busque y en varias partes se divide por la masa como en el LF tradicional
    vector3D epsilon; epsilon.load(ran64.gauss(0,1), ran64.gauss(0,1), 0); //Variable aleatoria
    vector3D deltaV= -alpha*Vprime+sqrt(alpha*(2-alpha)*kT/m)*epsilon;
    r+=(Vprime+deltaV*0.5)*dt; 
    V=Vprime+deltaV;

    //Algortimo leap-frog para rotación "determinista"
    omega += tau*dt/I;
    theta += omega*dt;
}

void Cuerpo::Arranque(double dt){
    //Algotimo leap-frog
    V-=F*(dt/(2*m));
    omega -= tau*(dt/(2*I));
}

void Cuerpo::Dibujese(void){
    cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
    <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

//----------  Implementar funciones de la clase Colisionador  --------------

void Colisionador::Inicie(void)
{
    for (int i = 0; i < N+1; i++)
    {
        for (int j = 0; j < N+1; j++)
        {
        xCundall[i][j] = sold[i][j] = 0;
        }
    }
}

void Colisionador::CalculeAllF(Cuerpo*Particulas){
    int i,j;
    //Borrar la fuerza de todas los Particulas
    for(i=0;i<N;i++){
        Particulas[i].BorreFuerza();
    };
    //Calcular fuerza de magus, fuerza central y fuerza contra la pared para cada partícula
    for(i=0;i<N;i++){
        CalculeF_magnus(Particulas[i]);
        CalculeF_ceentral(Particulas[i]);
        CalculeF_pared(Particulas[i],Particulas[N], xCundall[i][N], sold[i][N], dt);
        CalculeT_ind(Particulas[i]);
    };
    //Calcular Fuerza entre colisiones de las partículas
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            CalculeF_colision(Particulas[i],Particulas[j], xCundall[i][j], sold[i][j], dt);
        } 
    }
}

void Colisionador::CalculeF_colision(Cuerpo & Particula1,Cuerpo & Particula2, double &xCundall, double &sold, double dt){
    //Determinar si hay colisión
    vector3D r21=Particula2.r-Particula1.r;
    double rn=r21.norm();
    double R1=Particula1.R, R2=Particula2.R; 
    double s=(R1+R2)-rn;
    if(s>0){//Si hay colisión

        //Calcular el vector normal
        vector3D n = r21 / rn, t, k;
        t.load(n.y(), -n.x(), 0);
        k.load(0, 0, 1);

        // calculo de velocidades de contacto
        vector3D Rw;
        Rw.load(0, 0, R2 * Particula2.omega + R1 * Particula1.omega);
        vector3D Vc = (Particula2.V - Particula1.V) - (Rw ^ n);
        double Vcn = Vc * n, Vct = Vc * t;

        //Calcular fuerza Herzt
        double m1 = Particula1.m, m2 = Particula2.m;
        double m12 = m1 * m2 / (m1 + m2);
        double Fn=KH*pow(s,1.5) - GammaHertz * sqrt(s) * m12 * Vcn;

        // calculo fuerza tangencial (Cundall)
        xCundall += Vct * dt;
        double Ft = -Kcundall * xCundall;
        double Ftmax = mu * fabs(Fn);
        if (fabs(Ft) > Ftmax)
        {
        Ft = Ftmax * Ft / fabs(Ft);
        }

        // calcular y cargar fuerzas
        vector3D F1, F2, tau1, tau2;
        F2 = n * Fn + t * Ft; 
        tau2 = ((n * (-R2)) ^ F2); 
        F1 = F2 * -1;
        tau1 = ((n * R1) ^ F1);

        Particula2.SumeFuerza(F2, tau2 * k); 
        Particula1.SumeFuerza(F1, tau1 * k);
    }
    if(sold>=0 && s<0){xCundall=0;}  
    sold=s;
}

void Colisionador::CalculeF_pared(Cuerpo & Particula,Cuerpo & pared, double &xCundall, double &sold, double dt){
    //Fuerza elástica de la pared circular
    double rn=Particula.r.norm();
    double Rp=pared.R, R=Particula.R; 
    double s=(R+rn)-Rp;
    if(s>0){//Si hay colisión
        //Calcular el vector normal
        vector3D n=Particula.r*(1/rn), t, k;
        t.load(n.y(), -n.x(), 0);
        k.load(0, 0, 1);

        // calculo de velocidades de contacto
        vector3D Rw;
        Rw.load(0, 0, R * Particula.omega);
        vector3D Vc = (pared.V - Particula.V) - (Rw ^ n);
        double Vcn = Vc * n, Vct = Vc * t;

        //Calcular fuerza Herzt
        double m1 = pared.m, m2 = Particula.m;
        double m12 = m1 * m2 / (m1 + m2);
        double Fn=KH*pow(s,1.5) - GammaHertz * sqrt(s) * m12 * Vcn;

        // calculo fuerza tangencial (Cundall)
        xCundall += Vct * dt;
        double Ft = -Kcundall * xCundall;
        double Ftmax = mu * fabs(Fn);
        if (fabs(Ft) > Ftmax)
        {
        Ft = Ftmax * Ft / fabs(Ft);
        }

        // calcular y cargar fuerzas
        vector3D F1, F2, tau;
        F2 = n * Fn + t * Ft; 
        F1 = F2 * -1;
        tau = ((n * R) ^ F1);

        Particula.SumeFuerza(F1, tau * k);
    }
}
void Colisionador::CalculeF_magnus(Cuerpo & Particula)
{
    vector3D w;
    w.load(0,0,Particula.omega); //* para el resto de calculos es más conveniente trabajar con omega
    vector3D Fm = -0.5 * Cd * rho * Particula.A * Particula.R * (w ^ Particula.V); // Fuerza de Magnus
    Particula.SumeFuerza(Fm, 0); 
}

void Colisionador::CalculeF_ceentral(Cuerpo & Particula)
{
    vector3D Fc = -K * Particula.r; // Fuerza central
    Particula.SumeFuerza(Fc, 0); 
}

void Colisionador::CalculeT_ind(Cuerpo & Particula)
{
    vector3D T, Tv, Tfric, k, F_zero;
    F_zero.load(0, 0, 0);
    double tau0 = -KT * Particula.R * 10;

    //Torque por vemtilador
    Tv.load(0, 0, tau0);
    //Torque de fricción
    k.load(0, 0, 1);
    Tfric.load(0, 0, -Gamma * Particula.omega);

    T = Tv + Tfric;
    Particula.SumeFuerza(F_zero, T * k);
}
//--------------------------
void InicieAnimacion(const std::string& str) 
{
    cout << "set terminal gif animate" << endl;
    cout << "set output '"<< str <<"'" << endl;
    cout << "unset key" << endl;
    cout<<"set xrange["<<-Lx-2<<":"<<Lx+2<<"]"<<endl;
    cout<<"set yrange["<<-Ly-2<<":"<<Ly+2<<"]"<<endl;
    cout<<"set size ratio -1"<<endl;
    cout<<"set parametric"<<endl;
    cout<<"set trange [0:7]"<<endl;
    cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<1e-6<<"+"<<10<<"*cos(t),"<<1e-6<<"+"<<10<<"*sin(t)"; //Pared Circular
}
void TermineCuadro(void){
    cout<<endl;
}

int stringToBool(std::string str) 
{
    if (str == "False" || str == "0")
    {
        return 0;
    }

    else if (str == "true" || str == "1")
    {
        return 1;
    }

    else if (str == "Both" || str == "2")
    {
        return 2;
    }
    return 0;
}

void Imprimase(double t, Cuerpo *Polen, int gnuplot) 
{
    if ( gnuplot > 2 || gnuplot < 0)
    {
        gnuplot = 0;
    }
    switch (gnuplot)
    {
        case 0:
            for (int i = 0; i < N; i++) 
            {
                cout << t << " " << i << " " << Polen[i].Get_x() << " " << Polen[i].Get_y() << " " 
                << Polen[i].Get_Vx() << " " << Polen[i].Get_Vy() << Polen[i].Get_omega() << endl;
            }
            break;
        case 1:
            if (t == 0) InicieAnimacion("mdisc.gif");
            InicieCuadro();
            for (int i = 0; i < N; i++) Polen[i].Dibujese();
            TermineCuadro();
            break;
        case 2:
            if (t == 0) InicieAnimacion("mdisc.gif");
            InicieCuadro();
            for (int i = 0; i < N; i++) Polen[i].Dibujese();
            TermineCuadro();

            std::ofstream MyFile("mdisc.dat", std::ios::app);
            for (int i = 0; i < N; i++) 
            {
                MyFile << t << " " << i << " " << Polen[i].Get_x() << " " << Polen[i].Get_y() << " " 
                << Polen[i].Get_Vx() << " " << Polen[i].Get_Vy() << Polen[i].Get_omega() << endl;
            }
            MyFile.close();
            break;
    }
}