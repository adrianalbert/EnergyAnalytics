/* Sample C routines to accompany "Ratios of normal variables"

/* Below are C routines for studying the density and
distribution of the ratio of jointly normal variates z/w.
A sample main() is included. It can be used to relate
the distribution  of z/w to the standard form (a+x)/(b+y).
Compiled with gcc under DOS, the main() produces the
following output (given the default input):

    z/w is distributed as s+[(a+x)/(b+y)]/r, with
    a=0.500000,b=8.000000,r=-1.333333,s=1.000000

    Since b>4, the practical moments of z/w are:
    approx mean=0.951977, approx sig=0.095955

    To see if z/w is close to normal with those moments,
    compare g(t) with normal density having that mu and sig,
    or generate, say, 10,000 z/w's and see if
    Phi((zoverw()-mu)/sig) appears to be uniform in [0,1),
    using chisquare-99 test on frequencies in the intervals
         .00-.01, .01-.02,...,,99-1.00.
    Result: Chisquare_99 = 105.8200    (70<chisq<130 OK)

    To see if G is appropriate, generate 1,000,000 z/w's and see
    if G(zoverw()) seems to be uniform in [0,1) by chisquare test.
    Result: Chisquare_99 = 103.3132    (70<chisq<130 OK)



You may want to see if your platform produces that output,
and, more generally, write a main() to  use the  f,F g,G
functions and/or the zoverw() procedure to study particular
ratios of normal variates z/w of interest to you.
*/
/*
Because so many of the routines require the means, standard
deviations and correlation for z and w, they are made static
doubles under the names muz,muw,sigz,sigw and rho.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static unsigned long jc=123456789, jx=362436069;
static double muz,muw,sigz,sigw,rho,cov;

/* VNI is a macro for inline generation of a uniform v in -1<v<1.
  Combines congruential and xorshift RNGs
*/
#define VNI ( jc=69069*jc+123,jx^=jx<<13,jx^=jx>>17,jx^=jx<<5,\
             4.656612873e-10*(jc+jx)-1. )

/*-----------Procedure for generating random z/w------------- */
/* z and w are jointly normal, static muz,muw,sigz,sigw,rho */
double zoverw()
{double s=2,x,y,z,w;
 while(s>1) {x=VNI; y=VNI; s=x*x+y*y;}
 s=sqrt(-2*log(s)/s); x=s*x; y=s*y;
 z=muz+sigz*(x*sqrt(1-rho*rho)+rho*y);
 w=muw+sigw*y;
 return(z/w) ;
}

/*--- f(t): density of T=(a+x)/(b+y)-----------------------  */
/* checks for changes of aa and bb; repeated calls often vary only t.*/
double f(double aa, double bb, double t)
{long double s,st=0,x,w,v,i=1;
static double a=-1,b,c;
if(aa!=a | bb!=b)
{a=aa; b=bb; c=.31830988618379L*exp(-.5*(a*a+b*b));}
    x=(b+a*t)*(b+a*t)/(1+t*t);
    s=1+x; w=v=x; st=0;
    while(s!=st) s=(st=s)+(w*=v/(i+=2));
 return c*s/(1+t*t);
}

/*  Phi(x)=int(exp(-t^2/2)/sqrt(2*Pi),t=-infinity..x)       */
    double Phi(double x)
    {long double s=x,t=0,b=x,q=x*x,i=1;
    if(x<-8) return 0.; if (x>8)  return 1.;
    while(s!=t) s=(t=s)+(b*=q/(i+=2));
     return .5+s*exp(-.5*q-.91893853320467274178L);}

/*--F(z)=Prob((a+x)/(b+y)<z)=int(f(t),t=-infinity .. z),---------*/
/*  choose n, for Simpson's rule; default 100  */
 double F(double a, double b, double z)
{int i,n=100;
 double h=z/n,s=f(a,b,0)+f(a,b,z)+4*f(a,b,z-h);
 for(i=1;i<=n-3;i+=2) s+=4*f(a,b,i*h)+2*f(a,b,i*h+h);
 return Phi(a)+Phi(b)-2*Phi(a)*Phi(b)+s*h/3.;
}

/*----G(t)=Pr(z/w < t) based on assuming w always positive.-------
 Assumes static muz,muw,sigz,sigw,rho have been assigned.
 To test suitability of G, see if thousands of
 calls G(zoverw())  seem to yield uniform [0,1) values.
*/
double G(double t)
{double v=sqrt(sigz*sigz-2*t*sigz*sigw*rho+t*t*sigw*sigw);
return Phi((t*muw-muz)/v);
}

/*Density g(t)=G'(t) for z/w, assuming w always positive.
   Same input as G.   Resulting g often looks normal;
   if so, find mean m and variance sig^2 as done below,
   then compare that normal density with g(t), or generate
   many z/w's by means of the zoverw() routine and see if
   Phi((z/w-m)/sig) seems uniform in [0,1),
*/
double g(double t)
{double v=sqrt(sigz*sigz-2*t*sigz*sigw*rho+t*t*sigw*sigw);
 double w=(t*muw-muz)/v;
v=(sigz*(muw*sigz-muz*rho*sigw)+t*sigw*(muz*sigw-sigz*rho*muw))/(v*v*v);
 return .39894228*exp(-.5*w*w)*v;
}

int main()
{int i,k[100]={0};
double a,b,r,s,appmu,appsig;
/* assume, without loss of generality, that w has mean muw>=0. */
/* assign values to parameters for z and w; they are static doubles
   named    muz,muw,sigz,sigw,rho.   For example: */

 muz=30.5;muw=32;sigz=5;sigw=4;rho=.8;

/* get r,s,a,b so that  r*z/w-s=(a+x)/(b+y) */
 s=rho*sigz/sigw;
 r=sigw/(sigz*sqrt(1-rho*rho));
 b=muw/sigw;
 a=(muz/sigz-rho*muw/sigw)/sqrt(1-rho*rho);
 if(a<0) a=-a; r=-r;

 printf(" z/w is distributed as s+[(a+x)/(b+y)]/r, with\n");
 printf(" a=%f,b=%f,r=%f,s=%f\n\n",a,b,r,s);

 if(b>4){ appmu=a/(1.01*b-.2713);
     appsig=sqrt((a*a+1)/(b*b+.108*b-3.795)-appmu*appmu);
  printf(" Since b>4, the practical moments of z/w are:\n");
  printf(" approx mean=%f, approx sig=%f\n\n",
             (appmu=s+appmu/r),(appsig=appsig/fabs(r)));
  printf(" To see if z/w is close to normal with those moments,\n");
  printf(" compare g(t) with normal density having that mu and sig,\n");
  printf(" or generate, say, 10,000 z/w's and see if\n");
  printf(" Phi((zoverw()-mu)/sig) appears to be uniform in [0,1),\n");
  printf(" using chisquare-99 test on frequencies in the intervals\n");
  printf("         .00-.01, .01-.02,...,.99-1.00.\n");
  for(i=0;i<10000;i++) k[(int)(100*Phi((zoverw()-appmu)/appsig))]++;
  	  s=0; for(i=0;i<100;i++) {s+=(k[i]-100)*(k[i]-100)/100.; k[i]=0;}
	  printf("  Result: Chisquare_99 = %.4f    (70<chisq<130 OK)\n\n",s);
  printf(" To see if G is appropriate, generate 1,000,000 z/w's and see\n");
  printf(" if G(zoverw()) seems to be uniform in [0,1) by chisquare test.\n");
	  for(i=0;i<1000000;i++) k[(int)(100*G(zoverw()))]++;
	  s=0; for(i=0;i<100;i++) s+=(k[i]-10000)*(k[i]-10000)/10000.;
	  printf("  Result: Chisquare_99 = %.4f    (70<chisq<130 OK)\n\n",s);

	}

}





