// CodeCogs Adapted GNU General Public License Agreement
// Copyright (C) 2004-2005 CodeCogs, Zyba Ltd, Broadwood, Holford, TA5 1DU, England.
// This program is free software; you can redistribute it and/or modify it under
// the terms of the Adapted GNU General Public License as published by CodeCogs
// (www.codecogs.com/cart-3.htm). You must retain a copy of this licence in all copies.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the Adapted GNU General Public License for more details.
// 
// *** THIS SOFTWARE CAN NOT BE USED FOR COMMERCIAL GAIN.***
// ---------------------------------------------------------------------------------

#ifndef _CC_CALCULUS_SPECIAL_ERROR_FN_H_ 
   #define _CC_CALCULUS_SPECIAL_ERROR_FN_H_ 
   #ifndef cc_error 
   #include <assert.h> 
 #ifndef NDEBUG 
   #include <iostream> 
 #define cc_error(a) std::cerr << (a) << std::endl; 
   #else 
   #define cc_error(a) 
   #endif 
   #endif 

 #include <math.h> 
 #include <codecogs/algebra/polynomial/polynomial.h> 

 #define MINLOG 7.08396418532264106224E2  /* log 2^1022 */ 

 namespace Maths
{
 namespace Special
{

 // forward declaration of co-dependent function 
 double  errorFnC( double  a);

 //

double  errorFn(  double  x )
{
   static   double  T[] = {
   9.60497373987051638749E0,
   9.00260197203842689217E1,
   2.23200534594684319226E3,
   7.00332514112805075473E3,
   5.55923013010394962768E4
  };
   static   double  U[] = {
   /* 1.00000000000000000000E0,*/ 
   3.35617141647503099647E1,
   5.21357949780152679795E2,
   4.59432382970980127987E3,
   2.26290000613890934246E4,
   4.92673942608635921086E4
  };
   double  y, z;
   if ( fabs(x) > 1.0 )
     return ( 1.0 - errorFnC(x) );
  z = x * x;
  y = x * Algebra::Polynomial::polyEval( z, T, 4 )
        / Algebra::Polynomial::polyEval1( z, U, 5 );
   return ( y );
}


inline   double  errorFnC_exp( double  x)
{
   static   double  P[] = {
   2.46196981473530512524E-10,
   5.64189564831068821977E-1,
   7.46321056442269912687E0,
   4.86371970985681366614E1,
   1.96520832956077098242E2,
   5.26445194995477358631E2,
   9.34528527171957607540E2,
   1.02755188689515710272E3,
   5.57535335369399327526E2
  };
   static   double  Q[] = {
   /* 1.00000000000000000000E0,*/ 
   1.32281951154744992508E1,
   8.67072140885989742329E1,
   3.54937778887819891062E2,
   9.75708501743205489753E2,
   1.82390916687909736289E3,
   2.24633760818710981792E3,
   1.65666309194161350182E3,
   5.57535340817727675546E2
  };
   static   double  R[] = {
   5.64189583547755073984E-1,
   1.27536670759978104416E0,
   5.01905042251180477414E0,
   6.16021097993053585195E0,
   7.40974269950448939160E0,
   2.97886665372100240670E0
  };
   static   double  S[] = {
   /* 1.00000000000000000000E0,*/ 
   2.26052863220117276590E0,
   9.39603524938001434673E0,
   1.20489539808096656605E1,
   1.70814450747565897222E1,
   9.60896809063285878198E0,
   3.36907645100081516050E0
  };
   double  p,q;
   if ( x < 8.0 )
  {
    p = Algebra::Polynomial::polyEval( x, P, 8 );
    q = Algebra::Polynomial::polyEval1( x, Q, 8 );
  }
   else 
  {
    p = Algebra::Polynomial::polyEval( x, R, 5 );
    q = Algebra::Polynomial::polyEval1( x, S, 6 );
  }
   return  (p/q);
}


double  errorFnC( double  a)
{
   double  p,q,x,y,z;
   if ( a < 0.0 )
    x = -a;
   else 
    x = a;
   if ( x < 1.0 )
     return ( 1.0 - errorFn(a) );

  z = -a * a;
   if ( z < MINLOG )
  {
under:
    cc_error(  "errorFnC: UNDERFLOW"  );
     if ( a < 0 )
       return ( 2.0 );
     else 
       return ( 0.0 );
  }

  z = Maths::Arithmetic::expx2( a, -1 );

  y = errorFnC_exp(x) * z;
   if ( a < 0 )
    y = 2.0 - y;

   if ( y == 0.0 )
     goto  under;

   return (y);
}



};   // namespace Special 
};   // namespace Maths 
 #endif 
   //
