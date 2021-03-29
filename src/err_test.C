
#include  <iostream>
#include  <cmath>                          // to include sqrt(), etc.
#include  <cstdlib>                        // for atoi() and atof()
#include  <unistd.h>                       // for getopt()
#include  <fstream>                        // for file output
#include  <string>
using namespace std;

typedef double  real;                      // "real" as a general name for the
                                           // standard floating-point data type

real erff(real x);
real gammp(real a, real x);
void gcf(real *gammcf, real a, real x, real *gln);
void gser(real *gamser, real a, real x, real *gln);
real gammln(real xx);
real jfit(real x);

int main( )
{

  real xx = 0.5;
  printf("\n errorFn(%f) = %f", xx, erff(xx));

  return true;
}



//***********************************************************************
real erff(real x)
{

      return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);

}

//***********************************************************************
real gammp(real a, real x)
{

      real gammcf,gamser,gln;

      if (x < 0. || a <= 0.0) cerr << "bad arguments in gammp";

      if (x < a+1.) {
        gser(&gamser,a,x,&gln);
        return gamser;
      } else {
        gcf(&gammcf,a,x,&gln);
        return 1.0-gammcf;
      }
}


//***********************************************************************
void gcf(real *gammcf, real a, real x, real *gln)
{
      const int itmax = 100;
      const real eps = 3.e-7;
      const real fpmin = 1.0e-30;

      int i;

      real an,b,c,d,del,h;

      *gln = gammln(a);
      b = x+1.0-a;
      c = 1.0/fpmin;
      d = 1.0/b;
      h = d;

      for (i=1; i<=itmax; i++) {
        an = -i*(i-a);
        b += 2.0;
        d = an*d+b;
        if (fabs(d) < fpmin) d=fpmin;
        c = b+an/c;
        if (fabs(c) < fpmin) c=fpmin;
        d = 1.0/d;
        del = d*c;
        h *= del;
        if (fabs(del-1.0) < eps) break;
      }

      if (i> itmax) cerr << "a too large, ITMAX too small in gcf";

      *gammcf = exp(-x+a*log(x)-(*gln))*h;
}

//***********************************************************************
void gser(real *gamser, real a, real x, real *gln)
{
      const int  itmax = 100;
      const real eps = 3.0e-7;

      real ap,del,sum;

      *gln = gammln(a);

      if (x <= 0.0) {
        if (x < 0.0) cerr << "x < 0 in gser";
        *gamser = 0.0;
        return;
      } else {

        ap = a;
        sum = 1.0/a;
        del = sum;

        for (int n=1; n <= itmax; n++) {
          ++ap;
          del *= x/ap;
          sum += del;
          if (fabs(del) < fabs(sum)*eps) {
            *gamser = sum*exp(-x+a*log(x)-(*gln));
            return;
          }
        }

        cerr << "a too large, ITMAX too small in gser";

      }

}

//***********************************************************************
real gammln(real xx)
{

        real x,y,tmp,ser;
 
        const real cof[6]={76.18009172947146,-86.50532032941677, 
                      24.01409824083091,-1.231739572450155, 
                      0.1208650973866179e-2,-0.5395239384953e-5}; 

        //cof[0] = 76.18009173;
        //cof[1] = -86.50532033;
        //cof[2] = 24.01409822;
        //cof[3] = -1.231739516;
        //cof[4] = 0.120858003E-2;
        //cof[5] = -0.536382E-5;

        y = x = xx;
        tmp = x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;

        for (int j=0; j<=5; j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x); 

}

//***********************************************************************
real jfit(real x)
{

real jfitv, num, den, fx;
real a,b,c,d,e,f,g,h,i,j;
//      REAL*8 jfit, x, num, den, fx
//      REAL*8 a,b,c,d,e,f,g,h,i,j

// Fitting parameters.

     // PARAMETER (a=0.2021d0,b=0.04474d0,c=3.1621d0,
     //& d=0.5488d0,e=2.9965d0,f=0.7375d0,g=0.9217d0,
     //& h=0.2749d0,i=1.2583d0,j=2.3069d0)

     a=0.2021;
     b=0.04474;
     c=3.1621;
     d=0.5488;
     e=2.9965;
     f=0.7375;
     g=0.9217;
     h=0.2749;
     i=1.2583;
     j=2.3069;

      //num = DLOG(1.d0+a*x)*(b*(x**c)+d*(x**e))
      num = log(1.0+a*x)*(b*( pow(x,c) )+d*( pow(x,e) ));
      den = pow( 1.0 + f*( pow(x,g) ) + h*( pow(x,i) ) ,j);
      fx = log(1.0+x)-x/(1.0+x);

      jfitv = num/den;
      jfitv = jfitv/fx/fx;

      return jfitv;
//      END
}
