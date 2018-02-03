#include<stdlib.h>
#include<stdio.h>
#include<conio.h>
#include<math.h>

double EOS(double, int);//returns P as a function of rho
double EOSinv(double, int);//returns rho as a function of P
double EOSder(double, double, int); //returns the value of dp/d(rho)
double mDer(double, double);
double PDer(double, double, double, double);
double betaDer(double, double, double, double, double, double);
double rstart=0.01;
double rstop=12.96;
double rstep=0.01;
double r=0.01;
double K=7.29935;
double tau=5.0/3.0;
double a0=1.0;
double expLambda;
int i;
//all units are in km
double m[1500];
double rho[1500];
double P[1500];
double y[1500];
double H[1500];
double beta[1500];
double mRK[4];
double PRK[4];
double betaRK[4];
double mass, press, dens, Beta, h;

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
    P[0]=EOS(rho[0], 1);
    m[0]=4/3*M_PI*pow(rstart, 3.0)*rho[0];
    H[0]=a0*pow(rstart, 2.0);//a0 is arbitrarily chosen
    beta[0]=2*a0*rstart;
    y[0]=r*beta[0]/H[0];
    fprintf(fptr, "    r       m        P       rho       H        y");
    fprintf(fptr, "\n%f %f %f %f %f %f", r , m[0], P[0], rho[0], H[0], y[0]);
    for(i=1;i<1296;i++)
    {
        if(r==rstop)
            break;
        //RK4 step 1: calculating 'k1's and m, P, rho, beta, H values to be used in step 2 (refer to http://mathworld.wolfram.com/Runge-KuttaMethod.html)
        mRK[0]=rstep*mDer(r, rho[i-1]);
        PRK[0]=rstep*PDer(r, m[i-1], P[i-1], rho[i-1]);
        betaRK[0]=rstep*betaDer(r, m[i-1], P[i-1], rho[i-1], beta[i-1], H[i-1]);
        mass=m[i-1]+mRK[0]/2.0;
        press=P[i-1]+PRK[0]/2.0;
        Beta=beta[i-1]+betaRK[0]/2.0;
        dens=EOSinv(press, 1);
        h=H[i-1]+Beta*rstep;
        //RK4 step 2: calculating 'k2's and values for step 3
        mRK[1]=rstep*mDer(r+rstep/2.0, dens);
        PRK[1]=rstep*PDer(r+rstep/2.0, mass, press, dens);
        betaRK[1]=rstep*betaDer(r+rstep/2.0, mass, press, dens, Beta, h);
        mass=m[i-1]+mRK[1]/2.0;
        press=P[i-1]+PRK[1]/2.0;
        Beta=beta[i-1]+betaRK[1]/2.0;
        dens=EOSinv(press, 1);
        h=H[i-1]+Beta*rstep;
        //RK4 step 3: calculating 'k3's and values for step 4
        mRK[2]=rstep*mDer(r+rstep/2.0, dens);
        PRK[2]=rstep*PDer(r+rstep/2.0, mass, press, dens);
        betaRK[2]=rstep*betaDer(r+rstep/2.0, mass, press, dens, Beta, h);
        mass=m[i-1]+mRK[2];
        press=P[i-1]+PRK[2];
        Beta=beta[i-1]+betaRK[2];
        dens=EOSinv(press, 1);
        h=H[i-1]+Beta*rstep;
        //RK4 step 4: calculating 'k4's
        mRK[3]=rstep*mDer(r+rstep, dens);
        PRK[3]=rstep*PDer(r+rstep, mass, press, dens);
        betaRK[3]=rstep*betaDer(r+rstep, mass, press, dens, Beta, h);
        //calculating RK4 final results for this step
        m[i]=m[i-1]+mRK[0]/6.0+mRK[1]/3.0+mRK[2]/3.0+mRK[3]/6.0;
        P[i]=P[i-1]+PRK[0]/6.0+PRK[1]/3.0+PRK[2]/3.0+PRK[3]/6.0;
        beta[i]=beta[i-1]+betaRK[0]/6.0+betaRK[1]/3.0+betaRK[2]/3.0+betaRK[3]/6.0;
        rho[i]=EOSinv(P[i], 1);
        H[i]=H[i-1]+beta[i]*rstep;
        r=r+rstep;
        y[i]=r*beta[i]/H[i];
        fprintf(fptr, "\n%f %f %f %f %f %f", r , m[i], P[i], rho[i], H[i], y[i]);
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
double mDer(double r, double rho)
{
    return 4.0*M_PI*pow(r, 2.0)*rho;
}
double PDer(double r, double m, double P, double rho)
{
    return -1.0*(rho+P)*(m+4.0*M_PI*pow(r, 3.0)*P)/(pow(r, 2.0)-2.0*m*r);
}
double betaDer(double r, double m, double P, double rho, double beta, double H)
{
    expLambda=r/(r-2*m);
    return 2*expLambda*H*(-2*M_PI*(5*rho+9*P+(rho+P)/EOSder(P, rho, 1))+3/pow(r, 2.0)+2*expLambda*pow((m/pow(r, 2.0)+4*M_PI*r*P),2.0))+2*beta/r*expLambda*(-1+m/r+2*M_PI*pow(r, 2.0)*(rho-P));
}
