//Semi analytical solution of the two body problem
//Kepler's equation is numerically solved in order to find E(t)
//Analytical equations of the two body problem are then used to find r(t),f(t),x(t),y(t)

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define pi 3.14159265358979

//Initial conditions
const double mu = 0.9; //mass parameter mu = G(m1 + m2)
const double a = 1.0; //semi-major axis
const double e = 0.3; //eccentricity
const double tau = 0.0; //time of pericentre passage (phase)
const double n() //mean motion. Is being calculated through the 3rd law of Kepler
{
    return sqrt(mu/pow(a,3));
}

//time interval [t0,tmax] and time step dt for the numerical solution of Kepler's equation
const double t0 = 0.0;
const double dt = 0.01;
//Total number of rotations
const long N_rotations = 1000;

double atan2pi(double number)
{
    if (number < 0) return number + 2*pi;
    return number;
}

//Energy integral C
double E0()
{
    return -mu/(2*a);
}

//f(E) = E - e*sin(E) - M(t)
double f(double E, double t)
{
    double M = fmod(n()*(t - tau),2*pi);
    return E - e*sin(E) - M;
}

double df(double E)
{
    return 1 - e*cos(E);
}

//Numerical solution of Kepler's equations via Newton-Raphson method
double Kepler(double t)
{
    if (t < 1E-14) return t;

    short counter = 0;
    const short MAX_COUNTER = 50;
    double M = fmod(n()*(t - tau),2*pi);
    double E0,E1,error;
    E0 = M + (sin(M)/fabs(sin(M)))*0.85*e; //Danby's proposal
    do
    {
        E1 = E0 - f(E0,t)/df(E0);
        error = fabs(E0 - E1);
        E0 = E1;
        counter++;
    }
    while (error >= 1E-14 && counter < MAX_COUNTER);
    if (counter >= MAX_COUNTER) //Newton-Raphson did not converge
    {
        printf("Error in solving Kepler's equation. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    return E1;
}

int main()
{
    FILE *data_txy = fopen("analytic_txy.txt","w");
    FILE *data_trf = fopen("analytic_trf.txt","w");
    double Energy = E0();
    double t,E;
    double tmax = 2.0*pi*N_rotations;
    double r,f;
    double x,y,xdot,ydot;
    fprintf(data_txy,"t\t\tx\t\t\ty\t\t\tdE\n\n");
    fprintf(data_trf,"t\t\tr\t\t\tf\t\t\tdE\n\n");
    printf("Calculating orbit. Please wait...\n");
    clock_t begin, end;
    begin = clock();
    for (t = 0.0; t <= tmax; t += dt)
    {
        E = Kepler(t);
        r = a*(1 - e*cos(E));
        f = atan2pi(2*atan(sqrt((1+e)/(1-e))*tan(E/2.0)));
        x = a*(cos(E) - e);
        y = a*sqrt(1 - e*e)*sin(E);
        fprintf(data_txy,"%f\t%.16f\t%.16f\t%.16f\n",t,x,y,fabs(E0()-Energy));
        fprintf(data_trf,"%f\t%.16f\t%.16f\t%.16f\n",t,r,f,fabs(E0()-Energy));
        xdot = -a*n()*sin(E)/(1 - e*cos(E));
        ydot = a*sqrt(1-e*e)*cos(E)*n()/(1 - e*cos(E));
        Energy = 0.5*(xdot*xdot + ydot*ydot) - mu/r;
    }
    end = clock();
    double delta_t = (double)(end-begin)/CLOCKS_PER_SEC;
    printf("tmax >> %f\n",tmax);
    printf("Estimated calculation time : delta_t = %f\n",delta_t);
    printf("Data are written in 'analytic_txy.txt' and 'analytic_trf.txt'\n");
    fclose(data_txy);
    fclose(data_trf);
    return 0;
}
