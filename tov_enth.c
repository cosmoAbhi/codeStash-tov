#include<stdlib.h>
#include<stdio.h>
#include<conio.h>
#include<math.h>

double EOS(double, int);//returns P as a function of rho
double EOSinv(double, int);//returns rho as a function of P
double EOSder(double, double, int); //returns the value of dp/d(rho)
double rDer(double, double, double);
double mDer(double, double, double, double);
double PDer(double, double);
double betaDer(double, double, double, double, double, double);
double enthC=0.1451;
int numstep=1451;
double enthStep=-0.0001;
double enth=0.1451;
double rsmall=0.01;
double K=7.29935;
double tau=5.0/3.0;
double a0=1.0;
double expLambda;
int i;
//all units are in km
double r[1451];
double m[1451];
double rho[1451];
double P[1451];
double y[1451];
double H[1451];
double beta[1451];
double rRK[4];
double mRK[4];
double PRK[4];
double betaRK[4];
double rad, mass, press, dens, Beta, h;

int main()
{
    FILE *fptr;
    fptr=fopen("TOVoutputRK4new.txt","w");
    if(fptr == NULL)
    {
       printf("Error!");
       exit(1);
    }
    rho[0]=0.000741111;//km^(-2) central density
    r[0]=rsmall;
    P[0]=EOS(rho[0], 1);
    m[0]=4/3*M_PI*pow(rsmall, 3.0)*rho[0];
    H[0]=a0*pow(rsmall, 2.0);//a0 is arbitrarily chosen
    beta[0]=2*a0*rsmall;
    y[0]=rsmall*beta[0]/H[0];
    fprintf(fptr, "    h       r       m        P       rho       H        y");
    fprintf(fptr, "\n%f %f %f %f %f %f %f", enth, r[0], m[0], P[0], rho[0], H[0], y[0]);
    for(i=1;i<numstep;i++)
    {
        if(enth==0.0)
            break;
        //RK4 step 1: calculating 'k1's and m, P, rho, beta, H values to be used in step 2 (refer to http://mathworld.wolfram.com/Runge-KuttaMethod.html)
        rRK[0]=enthStep*rDer(r[i-1], m[i-1], P[i-1]);
        mRK[0]=enthStep*mDer(r[i-1], m[i-1], P[i-1], rho[i-1]);
        PRK[0]=enthStep*PDer(P[i-1], rho[i-1]);
        betaRK[0]=enthStep*betaDer(r[i-1], m[i-1], P[i-1], rho[i-1], beta[i-1], H[i-1]);
        rad=r[i-1]+rRK[0]/2.0;
        mass=m[i-1]+mRK[0]/2.0;
        press=P[i-1]+PRK[0]/2.0;
        Beta=beta[i-1]+betaRK[0]/2.0;
        dens=EOSinv(press, 1);
        h=H[i-1]+Beta*rDer(rad, mass, press)*enthStep;
        //RK4 step 2: calculating 'k2's and values for step 3
        rRK[1]=enthStep*rDer(rad, mass, press);
        mRK[1]=enthStep*mDer(rad, mass, press, dens);
        PRK[1]=enthStep*PDer(press, dens);
        betaRK[1]=enthStep*betaDer(rad, mass, press, dens, Beta, h);
        rad=r[i-1]+rRK[1]/2.0;
        mass=m[i-1]+mRK[1]/2.0;
        press=P[i-1]+PRK[1]/2.0;
        Beta=beta[i-1]+betaRK[1]/2.0;
        dens=EOSinv(press, 1);
        h=H[i-1]+Beta*rDer(rad, mass, press)*enthStep;
        //RK4 step 3: calculating 'k3's and values for step 4
        rRK[2]=enthStep*rDer(rad, mass, press);
        mRK[2]=enthStep*mDer(rad, mass, press, dens);
        PRK[2]=enthStep*PDer(press, dens);
        betaRK[2]=enthStep*betaDer(rad, mass, press, dens, Beta, h);
        rad=r[i-1]+rRK[1];
        mass=m[i-1]+mRK[2];
        press=P[i-1]+PRK[2];
        Beta=beta[i-1]+betaRK[2];
        dens=EOSinv(press, 1);
        h=H[i-1]+Beta*rDer(rad, mass, press)*enthStep;
        //RK4 step 4: calculating 'k4's
        rRK[3]=enthStep*rDer(rad, mass, press);
        mRK[3]=enthStep*mDer(rad,mass, press, dens);
        PRK[3]=enthStep*PDer(press, dens);
        betaRK[3]=enthStep*betaDer(rad, mass, press, dens, Beta, h);
        //calculating RK4 final results for this step
        r[i]=r[i-1]+rRK[0]/6.0+rRK[1]/3.0+rRK[2]/3.0+rRK[3]/6.0;
        m[i]=m[i-1]+mRK[0]/6.0+mRK[1]/3.0+mRK[2]/3.0+mRK[3]/6.0;
        P[i]=P[i-1]+PRK[0]/6.0+PRK[1]/3.0+PRK[2]/3.0+PRK[3]/6.0;
        beta[i]=beta[i-1]+betaRK[0]/6.0+betaRK[1]/3.0+betaRK[2]/3.0+betaRK[3]/6.0;
        rho[i]=EOSinv(P[i], 1);
        H[i]=H[i-1]+beta[i]*rDer(r[i], m[i], P[i])*enthStep;
        y[i]=r[i]*beta[i]/H[i];
        enth=enth+enthStep;
        fprintf(fptr, "\n%f %f %f %f %f %f %f", enth, r[i], m[i], P[i], rho[i], H[i], y[i]);
    }
    fclose(fptr);
}
double EOS(double dens, int choice)
{
    switch(choice)
    {
        case 1: return K*pow(dens, tau);
    }
    return 0.0;
}
double EOSinv(double pres, int choice)
{
    switch(choice)
    {
        case 1: return pow(pres/K, 1/tau);
    }
    return 0.0;
}
double EOSder(double P, double rho, int choice)
{
    switch(choice)
    {
        case 1: return tau*P/rho;
    }
    return 0.0;
}
double rDer(double r, double m, double P)
{
    return -r*(r-2.0*m)/(m+4.0*M_PI*pow(r, 3.0)*P);
}
double mDer(double r, double m, double P, double rho)
{
    return 4.0*M_PI*pow(r, 2.0)*rho*rDer(r, m, P);
}
double PDer(double P, double rho)
{
    return rho+P;
}
double betaDer(double r, double m, double P, double rho, double beta, double H)
{
    expLambda=r/(r-2*m);
    return (2*expLambda*H*(-2*M_PI*(5*rho+9*P+(rho+P)/EOSder(P, rho, 1))+3/pow(r, 2.0)+2*expLambda*pow((m/pow(r, 2.0)+4*M_PI*r*P),2.0))+2*beta/r*expLambda*(-1+m/r+2*M_PI*pow(r, 2.0)*(rho-P)))*rDer(r, m, P);
}
