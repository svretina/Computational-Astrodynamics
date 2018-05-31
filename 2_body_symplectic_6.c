//Simulating the 2-body problem via symplectic integration of 6th order (Yoshida)

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#define pi 3.14159265358979

//initial conditions
const double mu = 0.9;
double q10;  //x0
const double q20 = 0.0; //y0
const double p10 = 0.0; //vx0
double p20; //vy0
const double t0 = 0.0;
double dt;
//Total number of rotations
int N_rotations = 1000000;
float e;

double dTdp1(double p1)
{
    return p1;
}

double dTdp2(double p2)
{
    return p2;
}

double dVdq1(double q1, double q2)
{
    return (mu*q1)/((q1*q1 + q2*q2)*sqrt(q1*q1 + q2*q2));
}

double dVdq2(double q1, double q2)
{
    return (mu*q2)/((q1*q1 + q2*q2)*sqrt(q1*q1 + q2*q2));
}

double E0()
{
    return 0.5*(p10*p10 + p20*p20) - mu/sqrt(q10*q10 + q20*q20);
}

int main(int argc, char **argv)
{
  // external arguments passed:
  // e to write in filenames output data
  // q10
  // p20
  // dt
  // Number of rotations
  
  if (argc!=5){
    printf("Wrong number of arguments.");
      return 0;
    }
  else{
    e = atof(argv[1]);
    q10 = atof(argv[2]);
    p20 = atof(argv[3]);
    dt = atof(argv[4]);
    //N_rotations = atoi(argv[5]);
  }
  //printf("e=%.2f x0 = %f v_y0 = %f  dt= %f num_rot = %d\n", e,q10,p20,dt,N_rotations);
    //symplectic coefficients
    const double w1 = -1.1776799841788772038597;
    const double w2 = 0.23557321335935624273538;
    const double w3 = 0.78451361047756185129742;
    const double w0 = 1.0 - 2.0*(w1+w2+w3);

    const double c1 = w3/2.0;
    const double c8 = c1;

    const double c2 = (w3+w2)/2.0;
    const double c7 = c2;

    const double c3 = (w2+w1)/2.0;
    const double c6 = c3;

    const double c4 = (w1+w0)/2.0;
    const double c5 = c4;

    const double d1 = w3;
    const double d7 = d1;
    const double d2 = w2;
    const double d6 = d2;
    const double d3 = w1;
    const double d5 = d3;
    const double d4 = w0;

    double q1t,q2t,p1t,p2t;
    double q1_temp1, q2_temp1, p1_temp1, p2_temp1;
    double q1_temp2, q2_temp2, p1_temp2, p2_temp2;
    double q1_temp3, q2_temp3, p1_temp3, p2_temp3;
    double q1_temp4, q2_temp4, p1_temp4, p2_temp4;
    double q1_temp5, q2_temp5, p1_temp5, p2_temp5;
    double q1_temp6, q2_temp6, p1_temp6, p2_temp6;
    double q1_temp7, q2_temp7;
    double q1,p1,q2,p2;

    double tmax = 2*pi*N_rotations;

    q1t = q10;
    p1t = p10;
    q2t = q20;
    p2t = p20;
    
    char method[] = "symplectic";
    char extension[] = ".txt";
    char *eccentricity;
    eccentricity = argv[1];
    char *step;
    char *rotations;
    rotations = argv[5];
    step = argv[4];
    char outputfile[strlen(method)+strlen(extension)+strlen(eccentricity)+strlen(step)+strlen(rotations)+4];
    snprintf( outputfile , sizeof(outputfile ), "%s_%s_%s_%s%s", method,eccentricity,step,rotations,extension );

    int counter;
    FILE *fp = fopen( outputfile,"w");
    //fprintf(fp,"\n\n t\t\tx\t\t\ty\t\t\tdE\n\n");
    fprintf(fp,"%f,%.16f,%.16f,%.16f\n",t0,q10,p10,0.0);
    //printf("Calculating orbit. Please wait...\n");
    clock_t begin, end;
    begin = clock();
    for (double t = t0; t <= tmax; t += dt)
    {
        q1_temp1 = q1t + c1*dt*dTdp1(p1t);
        q2_temp1 = q2t + c1*dt*dTdp2(p2t);
        p1_temp1 = p1t - d1*dt*dVdq1(q1_temp1,q2_temp1);
        p2_temp1 = p2t - d1*dt*dVdq2(q1_temp1,q2_temp1);

        q1_temp2 = q1_temp1 + c2*dt*dTdp1(p1_temp1);
        q2_temp2 = q2_temp1 + c2*dt*dTdp2(p2_temp1);
        p1_temp2 = p1_temp1 - d2*dt*dVdq1(q1_temp2,q2_temp2);
        p2_temp2 = p2_temp1 - d2*dt*dVdq2(q1_temp2,q2_temp2);

        q1_temp3 = q1_temp2 + c3*dt*dTdp1(p1_temp2);
        q2_temp3 = q2_temp2 + c3*dt*dTdp2(p2_temp2);
        p1_temp3 = p1_temp2 - d3*dt*dVdq1(q1_temp3,q2_temp3);
        p2_temp3 = p2_temp2 - d3*dt*dVdq2(q1_temp3,q2_temp3);

        q1_temp4 = q1_temp3 + c4*dt*dTdp1(p1_temp3);
        q2_temp4 = q2_temp3 + c4*dt*dTdp2(p2_temp3);
        p1_temp4 = p1_temp3 - d4*dt*dVdq1(q1_temp4,q2_temp4);
        p2_temp4 = p2_temp3 - d4*dt*dVdq2(q1_temp4,q2_temp4);

        q1_temp5 = q1_temp4 + c5*dt*dTdp1(p1_temp4);
        q2_temp5 = q2_temp4 + c5*dt*dTdp2(p2_temp4);
        p1_temp5 = p1_temp4 - d5*dt*dVdq1(q1_temp5,q2_temp5);
        p2_temp5 = p2_temp4 - d5*dt*dVdq2(q1_temp5,q2_temp5);

        q1_temp6 = q1_temp5 + c6*dt*dTdp1(p1_temp5);
        q2_temp6 = q2_temp5 + c6*dt*dTdp2(p2_temp5);
        p1_temp6 = p1_temp5 - d6*dt*dVdq1(q1_temp6,q2_temp6);
        p2_temp6 = p2_temp5 - d6*dt*dVdq2(q1_temp6,q2_temp6);

        q1_temp7 = q1_temp6 + c7*dt*dTdp1(p1_temp6);
        q2_temp7 = q2_temp6 + c7*dt*dTdp2(p2_temp6);

        p1 = p1_temp6 - d7*dt*dVdq1(q1_temp7,q2_temp7);
        p2 = p2_temp6 - d7*dt*dVdq2(q1_temp7,q2_temp7);
        q1 = q1_temp7 + c8*dt*dTdp1(p1);
        q2 = q2_temp7 + c8*dt*dTdp2(p2);

        double Energy = 0.5*(p1*p1 + p2*p2) - mu/sqrt(q1*q1 + q2*q2);
	counter = (int) (t/(2.*pi));
	  if ((0<counter && counter<1) || (1000<counter && counter<1001) || (10000<counter && counter<10001) ||  (100000<counter && counter<100001) || (999999<counter && counter<1000000)) {
	    fprintf(fp,"%f,%.16f,%.16f,%.16f\n",t+dt,q1,q2,fabs(E0()-Energy));
	  }
        q1t = q1;
        p1t = p1;
        q2t = q2;
        p2t = p2;
    }
    end = clock();
    double delta_t = (double)(end-begin)/CLOCKS_PER_SEC;
    
    FILE *fp2 = fopen( "time.txt","a");
    fprintf(fp2,"%s,%s,%f,%s,%f\n",method,eccentricity,dt,rotations,delta_t);

    fclose(fp2);
    fclose(fp);
    return 0;
}
