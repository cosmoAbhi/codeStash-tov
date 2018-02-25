#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <math.h>
#include <sstream>
using namespace std;

double EOS(double, int);//returns P as a function of rho
double EOSinv(double, int);//returns rho as a function of P
double EOSder(double, double, int); //returns the value of dp/d(rho)

double rstart=0.1;
double rstop=15.0;
double rstep=0.1;
double r=0.1;
double K=7.29935;
double tau=5/3;
double a0=1.0;
double expLambda, C0, C1, nuPrime, yDer;
int i, readCount=1, mod, store=0, choice, EOScount=0;
double read, frac;
//all units are in km
double m[150];
double rho[150];
double P[150];
double y[150];
double H[150];
double beta[150];

int ser[1000];
double nB[1000];
double pEOS[1000];
double rhoEOS[1000];

int main()
{
    cout<<"Specify EOS choice: 1 for polytrope, 2 for EOS.txt: ";
    cin>>choice;
    if(choice==2)
    {
        ifstream infile;
        infile.open("EOS.txt");

        while (!infile.eof())
        {
            infile>> ser[EOScount];
            infile>> nB[EOScount];
            infile>> rhoEOS[EOScount];
            infile>> pEOS[EOScount];
            rhoEOS[EOScount]*=7.4237e-19;
            pEOS[EOScount++]*=8.2601e-40;
        }

        EOScount--;
        infile.close();
    }
    ofstream write;
    write.open("TOVtry_output.txt");

    rho[0]=0.000741111;//km^(-2) central density
    P[0]=EOS(rho[0], choice);
    m[0]=4/3*M_PI*pow(rstart, 3.0)*rho[0];
    //y[0]=2.0; //central value of y= l
    H[0]=a0*pow(rstart, 2.0);//a0 is arbitrarily chosen
    beta[0]=2*a0*rstart;
    y[0]=r*beta[0]/H[0];
    write<<" r      m           P               rho          H        y"<<endl;
    write<<r<<" "<< m[0]<<" "<<P[0]<<"  "<<rho[0]<<"    "<<H[0]<<"  "<<y[0]<<endl;
    for(i=1;i<150;i++)
    {
        if(r==rstop)
            break;
        m[i]=m[i-1]+4*M_PI*pow(r, 2.0)*rho[i-1]*rstep;
        P[i]=P[i-1]-(rho[i-1]+P[i-1])*(m[i-1]+4*M_PI*pow(r, 3.0)*P[i-1])/(pow(r, 2.0)-2*m[i-1]*r)*rstep;
        rho[i]=EOSinv(P[i], choice);
        expLambda=r/(r-2*m[i-1]);
        //C1=2/r+expLambda*(2*m[i-1]/pow(r, 2.0)+4*M_PI*r*(P[i-1]-rho[i-1]));
        //nuPrime=2*(m[i-1]+4*M_PI*pow(r, 3.0)*P[i-1])/(pow(r, 2.0)-2*m[i-1]*r);
        //C0=expLambda*(-2.0*(2.0+1.0)/pow(r, 2.0)+4*M_PI*(rho[i-1]+P[i-1])/der(i-1, 1)+4*M_PI*(5*rho[i-1]+9*P[i-1]))-pow(nuPrime, 2.0);
        //yDer=-y[i-1]*(y[i-1]-1)/r-C1*y[i-1]-r*C0*y[i-1];
        //y[i]=y[i-1]+yDer*rstep;
        H[i]=H[i-1]+beta[i-1]*rstep;
        beta[i]=beta[i-1]+2*expLambda*H[i-1]*(-2*M_PI*(5*rho[i-1]+9*P[i-1]+(rho[i-1]+P[i-1])/EOSder(P[i-1], rho[i-1], choice))+3/pow(r, 2.0)+2*expLambda*pow((m[i-1]/pow(r, 2.0)+4*M_PI*r*P[i-1]),2.0))+2*beta[i-1]/r*expLambda*(-1+m[i-1]/r+2*M_PI*pow(r, 2.0)*(rho[i-1]-P[i-1]));
        y[i]=r*beta[i]/H[i];
        r=r+rstep;
        write<<r<<" "<< m[i]<<" "<<P[i]<<"  "<<rho[i]<<"    "<<H[i]<<"  "<<y[i]<<endl;
    }
    write.close();
    return 0;
}
double EOS(double dens, int choice)
{
    switch(choice)
    {
        case 1: return K*pow(dens, tau);
        case 2: for(i=1;i<EOScount;i++)
                {
                    if(dens>rhoEOS[i-1]&&dens<rhoEOS[i+1])
                    {
                        frac=(dens-rhoEOS[i-1])/(rhoEOS[i]-rhoEOS[i-1]);
                        return (frac*(pEOS[i]-pEOS[i-1])+pEOS[i-1]);
                    }

                }
    }
    return 0.0;
}
double EOSinv(double pres, int choice)
{
    switch(choice)
    {
        case 1: return pow(pres/K, 1/tau);
        case 2: for(i=1;i<EOScount;i++)
                {
                    if(pres>pEOS[i-1]&&pres<pEOS[i+1])
                    {
                        frac=(pres-pEOS[i-1])/(pEOS[i]-pEOS[i-1]);
                        return (frac*(rhoEOS[i]-rhoEOS[i-1])+rhoEOS[i-1]);
                    }

                }
    }
    return 0.0;
}
double EOSder(double P, double rho, int choice)
{
    switch(choice)
    {
        case 1: return tau*P/rho;
        case 2: for(i=1;i<EOScount;i++)
                {
                    if(P>pEOS[i-1]&&P<pEOS[i+1])
                    {
                        frac=(P-pEOS[i-1])/(pEOS[i]-pEOS[i-1]);
                        if(frac<0.5)
                            return (P-pEOS[i-1])/(rho-rhoEOS[i-1]);
                        else
                            return (pEOS[i]-P)/(rhoEOS[i]-rho);
                    }

                }
    }
    return 0.0;
}
