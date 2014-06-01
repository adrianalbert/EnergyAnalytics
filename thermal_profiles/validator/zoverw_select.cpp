
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
 double F_ratio(double a, double b, double z)
{int i,n=100;
 double h=z/n,s=f(a,b,0)+f(a,b,z)+4*f(a,b,z-h);
 for(i=1;i<=n-3;i+=2) s+=4*f(a,b,i*h)+2*f(a,b,i*h+h);
 return Phi(a)+Phi(b)-2*Phi(a)*Phi(b)+s*h/3.;
}
