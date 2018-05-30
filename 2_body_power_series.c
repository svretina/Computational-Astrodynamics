//Integrating the 2-body problem using Taylor series method

#include<stdio.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>


#define pi 3.14159265358979

//initial conditions
const double mu = 0.9;
double x0;
const double y_0 = 0.0;
const double z0 = 0.0; //vx0
double k0; //vy0
int N; //number of Taylor terms
const double t0 = 0.0;
double dt;
//Total number of rotations
int N_rotations;
float e;
//energy integral
double E0()
{
    return 0.5*(z0*z0 + k0*k0) - mu/sqrt(x0*x0 + y_0*y_0);
}

int main(int argc, char **argv)
{

  if (argc!=7){
    printf("Wrong number of arguments.");
    return 0;
  }
  else{
    e = atof(argv[1]);
    x0 = atof(argv[2]);
    k0 = atof(argv[3]);
    dt = atof(argv[4]);
    N_rotations = atoi(argv[5]);
    N = atoi(argv[6]);
  }
    double a[N+1], b[N+1], c[N+1], d[N+1];
    double A[N+1], B[N+1], C[N+1], D[N+1], E[N+1], F[N+1];
    a[0] = x0;
    b[0] = y_0;
    c[0] = z0;
    d[0] = k0;
    double Energy = E0();
    double tmax = 2.0*pi*N_rotations;


    char method[] = "series";
    char extension[] = ".txt";
    char *eccentricity;
    eccentricity = argv[1];
    char *step;
    char *rotations;
    rotations = argv[5];
    step = argv[4];
    char *terms;
    terms = argv[6];
    char outputfile[strlen(method)+strlen(extension)+strlen(eccentricity)+strlen(step)+strlen(rotations)+strlen(terms)+5];
    snprintf( outputfile , sizeof(outputfile ), "%s_%s_%s_%s_%s%s", method,eccentricity,step,rotations,terms,extension );


    
    FILE *fp1 = fopen(outputfile,"w");
    fprintf(fp1,"t\t\tx\t\t\ty\t\t\tdE\n\n");
    fprintf(fp1,"%f\t%.16f\t%.16f\t%.16f\n",t0,x0,y_0,fabs(Energy-E0()));
    //printf("Calculating orbit. Please wait...\n");
    clock_t begin, end;
    begin = clock();
    for (double t = t0; t <= tmax; t += dt)
    {
        //initialize Taylor expansions
        double sumx = a[0];
        double sumy = b[0];
        double sumz = c[0];
        double sumk = d[0];

        a[1] = c[0];
        b[1] = d[0];
        A[0] = a[0]*a[0];
        B[0] = b[0]*b[0];
        C[0] = A[0]+B[0];
        D[0] = pow(C[0],-3.0/2.0);
        E[0] = a[0]*D[0];
        F[0] = b[0]*D[0];
        c[1] = -mu*E[0];
        d[1] = -mu*F[0];

        for (int n = 1; n < N; ++n)
        {
            a[n+1] = c[n]/(n+1.);
            b[n+1] = d[n]/(n+1.);

            //calculate An,Bn,Cn,Dn,En,Fn
            double sumA = 0.0;
            double sumB = 0.0;
            double sumC = 0.0;
            double sumD = 0.0;
            double sumE = 0.0;
            double sumF = 0.0;
            for (int i = 0; i <= n; ++i)
            {
                sumA += a[i]*a[n-i]; //build x^2
                sumB += b[i]*b[n-i]; //build y^2
            }
            sumC = sumA + sumB; //build x^2 + y^2
            A[n] = sumA;
            B[n] = sumB;
            C[n] = sumC;
            for (int i = 0; i <= n-1; ++i)
            {
                sumD += ((n-i)*(-3.0/2.0) - i)*C[n-i]*D[i]; //(build x^2 + y^2)^(-3/2)
            }
            D[n] = sumD/(n*C[0]);
            for (int i = 0; i <= n; ++i)
            {
                sumE += a[i]*D[n-i]; //build x*(x^2 + y^2)^(-3/2)
                sumF += b[i]*D[n-i]; //build y*(x^2 + y^2)^(-3/2)
            }
            E[n] = sumE;
            F[n] = sumF;

            c[n+1] = -mu*E[n]/(n+1.0);
            d[n+1] = -mu*F[n]/(n+1.0);

            //append to Taylor expansions
            sumx += a[n]*pow(dt,n);
            sumy += b[n]*pow(dt,n);
            sumz += c[n]*pow(dt,n);
            sumk += d[n]*pow(dt,n);
        }
        //keep track of the total energy
        Energy = 0.5*(sumz*sumz + sumk*sumk) - mu/sqrt(sumx*sumx + sumy*sumy);
        fprintf(fp1,"%f,%.16f,%.16f,%.16f\n",t+dt,sumx,sumy,fabs(E0()-Energy));
        //update initial conditions for the next t
        a[0]=sumx;
        b[0]=sumy;
        c[0]=sumz;
        d[0]=sumk;
    }
    end = clock();
    double delta_t = (double)(end-begin)/CLOCKS_PER_SEC;
    FILE *fp2 = fopen( "time.txt","a");
    fprintf(fp2,"%s,%s,%f,%s,%d,%f\n",method,eccentricity,dt,rotations,N,delta_t);

    
    //printf(">> tmax = %f\n",tmax);
    //printf("Estimated calculation time : delta_t = %f\n",delta_t);
    //printf("Data were written in 'series.txt'\n");
    fclose(fp2);
    fclose(fp1);
    return 0;
}
