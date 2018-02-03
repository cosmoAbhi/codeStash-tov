#include<stdlib.h>
#include<stdio.h>
#include<conio.h>
#include<math.h>
//trial comment
//This code calculates the critical enthalpy of a neutron star
double EOS(double, int);
double EOSinv(double, int);
double K=7.29935;
double tau=5.0/3.0;
double rhoC=0.000741111;
double hC=0.0;
double pC, psmall, pstep, p, denom;
int i;

int main()
{
    pC=EOS(rhoC, 1);
    psmall=0.0;//pC/10000000.0;
    pstep=(pC-psmall)/1000000000.0;
    p=psmall;
    for(i=0;i<1000000000n;i++)
    {
        if(p==pC)
            break;
        p=p+pstep/2.0;
        denom=EOSinv(p, 1)+p;
        hC=hC+pstep/denom;
        p=p+pstep/2.0;
    }
    printf("pC= %f", pC);
    printf("\nhC= %f", hC);
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
