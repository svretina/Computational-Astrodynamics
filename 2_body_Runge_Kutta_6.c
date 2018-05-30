//Numerical solution of the two body problem (Runge-Kutta 4)

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#define pi 3.14159265358979

//initial conditions
const double mu = 0.9; //mass parameter mu = G(m1 + m2)
const double x0 = 0.7;
const double y0 = 0.0;
const double vx0 = 0.0;
const double vy0 = 1.292837411057;
//time parameters
const double t0 = 0.0;
const double dt = 0.01;
//Total number of iterations
const long N_rotations = 1000;

//Energy integral C
double E0()
{
    return 0.5*(vx0*vx0 + vy0*vy0) - mu/sqrt(x0*x0 + y0*y0);
}

//ODEs system
////////////////////////////////////////////////////////////////
double f1(double t, double x, double y, double vx, double vy) //
{                                                             //
    return vx;                                                //
}                                                             //
                                                              //
double f2(double t, double x, double y, double vx, double vy) //
{                                                             //
    return vy;                                                //
}                                                             //
                                                              //
double f3(double t, double x, double y, double vx, double vy) //
{                                                             //
    double temp = pow(sqrt(x*x + y*y),3);                     //
    return -(mu/temp)*x;                                      //
}                                                             //
                                                              //
double f4(double t, double x, double y, double vx, double vy) //
{                                                             //
    double temp = pow(sqrt(x*x + y*y),3);                     //
    return -(mu/temp)*y;                                      //
}                                                             //
////////////////////////////////////////////////////////////////

//Runge-Kutta 4 method for a system of 4 ODEs
void RK4(double t, double *x, double *y, double *vx, double *vy, double dt)
{
    double k1=f1(t,*x,*y,*vx,*vy);
    double l1=f2(t,*x,*y,*vx,*vy);
    double m1=f3(t,*x,*y,*vx,*vy);
    double n1=f4(t,*x,*y,*vx,*vy);

    double k2=f1(t+dt/2, *x+(dt/2)*k1, *y+(dt/2)*l1, *vx+(dt/2)*m1, *vy+(dt/2)*n1);
    double l2=f2(t+dt/2, *x+(dt/2)*k1, *y+(dt/2)*l1, *vx+(dt/2)*m1, *vy+(dt/2)*n1);
    double m2=f3(t+dt/2, *x+(dt/2)*k1, *y+(dt/2)*l1, *vx+(dt/2)*m1, *vy+(dt/2)*n1);
    double n2=f4(t+dt/2, *x+(dt/2)*k1, *y+(dt/2)*l1, *vx+(dt/2)*m1, *vy+(dt/2)*n1);

    double k3=f1(t+dt/2, *x+(dt/2)*k2, *y+(dt/2)*l2, *vx+(dt/2)*m2, *vy+(dt/2)*n2);
    double l3=f2(t+dt/2, *x+(dt/2)*k2, *y+(dt/2)*l2, *vx+(dt/2)*m2, *vy+(dt/2)*n2);
    double m3=f3(t+dt/2, *x+(dt/2)*k2, *y+(dt/2)*l2, *vx+(dt/2)*m2, *vy+(dt/2)*n2);
    double n3=f4(t+dt/2, *x+(dt/2)*k2, *y+(dt/2)*l2, *vx+(dt/2)*m2, *vy+(dt/2)*n2);

    double k4=f1(t+dt, *x+dt*k3, *y+dt*l3, *vx+dt*m3, *vy+dt*n3);
    double l4=f2(t+dt, *x+dt*k3, *y+dt*l3, *vx+dt*m3, *vy+dt*n3);
    double m4=f3(t+dt, *x+dt*k3, *y+dt*l3, *vx+dt*m3, *vy+dt*n3);
    double n4=f4(t+dt, *x+dt*k3, *y+dt*l3, *vx+dt*m3, *vy+dt*n3);

    *x = *x + (dt/6)*(k1+2*k2+2*k3+k4);
    *y = *y + (dt/6)*(l1+2*l2+2*l3+l4);
    *vx = *vx + (dt/6)*(m1+2*m2+2*m3+m4);
    *vy = *vy + (dt/6)*(n1+2*n2+2*n3+n4);
}

int main()
{

    FILE *data_txy = fopen("RK4.txt","w");
    fprintf(data_txy,"t\t\tx\t\t\ty\t\t\tdE\n\n");

    double Energy = E0();
    double x, y, vx, vy, t;
    double tmax = 2.0*pi*N_rotations;
    //initialize RK4 method
    t = t0;
    x = x0;
    y = y0;
    vx = vx0;
    vy = vy0;

    printf("Calculating orbit. Please wait...\n");
    clock_t begin, end;
    begin = clock();
    for (t = t0; t <= tmax; t += dt)
    {
        fprintf(data_txy,"%f\t%.16f\t%.16f\t%.16f\n",t,x,y,fabs(E0()-Energy));
        RK4(t, &x, &y, &vx, &vy, dt);
        Energy = 0.5*(vx*vx + vy*vy) - mu/sqrt(x*x + y*y);
    }
    end = clock();
    double delta_t = (double)(end-begin)/CLOCKS_PER_SEC;
    fprintf(data_txy,"%f\t%.16f\t%.16f\t%.16f\n",t,x,y,fabs(E0()-Energy));
    printf(">> tmax = %f\n",tmax);
    printf("Estimated calculation time : delta_t = %f\n",delta_t);
    printf("Results were written in 'RK4.txt'\n");
    fclose(data_txy);
    return 0;
}
