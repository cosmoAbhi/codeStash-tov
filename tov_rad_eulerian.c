#include<stdio.h>
#include<conio.h>
#include<math.h>

double EOS(double, int);//returns P as a function of rho
double EOSinv(double, int);//returns rho as a function of P
double der(int, int); //returns the value of dp/d(rho)
double rstart=0.1;
double rstop=15.0;
double rstep=0.1;
double r=0.1;
double K=7.29935;
double tau=5/3;
double a0=1.0;
double expLambda, C0, C1, nuPrime, yDer;
int i, readCount=1, mod, store=0;
double read;
//all units are in km
double m[150];
double rho[150];
double P[150];
double y[150];
double H[150];
double beta[150];

double nB[1000];
double Peos[1000];
double rhoEos[1000];

int main()
{
    FILE *fptr;
    FILE *eosptr;
    fptr=fopen("TOVoutput.txt","w");
    eosptr=fopen("EOS.txt","r");

    //reading columns from EOS.txt and storing in arrays:
    while(fscanf(eosptr, "%e", &read)==1)
    {
            printf("\n %e",read);
            mod=readCount%4;
            switch(mod)
            {
                case 0: readCount++;
                        break;
                case 1: nB[store]=read;
                        readCount++;
                        //printf("\n %e",read);
                        break;
                case 2: Peos[store]=read;
                        readCount++;
                        //printf("\n %e",read);
                        break;
                case 3: rhoEos[store++]=read;
                        //printf("\n %e",read);
                        readCount++;
            }
    }
    //printf("%e", 2.720E-14);
    //solving the TOV equations utwards from the center:
    if(fptr == NULL)
    {
       printf("Error!");
       exit(1);
    }
    rho[0]=0.000741111;//km^(-2) central density
    P[0]=EOS(rho[0], 1);
    m[0]=4/3*M_PI*pow(rstart, 3.0)*rho[0];
    //y[0]=2.0; //central value of y= l
    H[0]=a0*pow(rstart, 2.0);//a0 is arbitrarily chosen
    beta[0]=2*a0*rstart;
    fprintf(fptr, "%f %f %f %f %f", r , m[0], P[0], rho[0], H[0]);
    for(i=1;i<150;i++)
    {
        if(r==rstop)
            break;
        m[i]=m[i-1]+4*M_PI*pow(r, 2.0)*rho[i-1]*rstep;
        P[i]=P[i-1]-(rho[i-1]+P[i-1])*(m[i-1]+4*M_PI*pow(r, 3.0)*P[i-1])/(pow(r, 2.0)-2*m[i-1]*r)*rstep;
        rho[i]=EOSinv(P[i], 1);
        expLambda=r/(r-2*m[i-1]);
        //C1=2/r+expLambda*(2*m[i-1]/pow(r, 2.0)+4*M_PI*r*(P[i-1]-rho[i-1]));
        //nuPrime=2*(m[i-1]+4*M_PI*pow(r, 3.0)*P[i-1])/(pow(r, 2.0)-2*m[i-1]*r);
        //C0=expLambda*(-2.0*(2.0+1.0)/pow(r, 2.0)+4*M_PI*(rho[i-1]+P[i-1])/der(i-1, 1)+4*M_PI*(5*rho[i-1]+9*P[i-1]))-pow(nuPrime, 2.0);
        //yDer=-y[i-1]*(y[i-1]-1)/r-C1*y[i-1]-r*C0*y[i-1];
        //y[i]=y[i-1]+yDer*rstep;
        H[i]=H[i-1]+beta[i-1]*rstep;
        beta[i]=beta[i-1]+2*expLambda*H[i-1]*(-2*M_PI*(5*rho[i-1]+9*P[i-1]+(rho[i-1]+P[i-1])/der(i-1, 1))+3/pow(r, 2.0)+2*expLambda*pow((m[i-1]/pow(r, 2.0)+4*M_PI*r*P[i-1]),2.0))+2*beta[i-1]/r*expLambda*(-1+m[i-1]/r+2*M_PI*pow(r, 2.0)*(rho[i-1]-P[i-1]));
        r=r+rstep;
        fprintf(fptr, "\n%f %f %f %f %f", r , m[i], P[i], rho[i], H[i]);
    }
    fclose(fptr);
}
double EOS(double dens, int choice)
{
    switch(choice)
    {
        case 1: return K*pow(dens, tau);
                break;
    }
}
double EOSinv(double pres, int choice)
{
    switch(choice)
    {
        case 1: return pow(pres/K, 1/tau);
                break;
    }
}
double der(int i, int choice)
{
    switch(choice)
    {
        case 1: return tau*P[i]/rho[i];
                break;
    }
}
