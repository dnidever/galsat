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

#ifndef _CC_MATHS_ARITHMETIC_EXPX2_H_ 
   #define _CC_MATHS_ARITHMETIC_EXPX2_H_ 
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
 #include <float.h> 

 #define MAXLOG   7.08396418532264106224E2 

 #define M 128.0 
   #define MINV .0078125 

 //

double  expx2(  double  x,  int  sign=1)   // Cephes name: expx2 
{
   double  u, u1, m, f;
  x = fabs(x);
   if  (sign < 0)
    x = -x;

 /* Represent x as an exact multiple of M plus a residual. 
      M is a power of 2 chosen so that exp(m * m) does not overflow 
      or underflow and so that |x - m| is small.  */ 
  m = MINV * floor(M * x + 0.5);
  f = x - m;

   /* x^2 = m^2 + 2mf + f^2 */ 
  u = m * m;
  u1 = 2*m*f  +  f*f;
   if  (sign < 0)
   {
      u = -u;
      u1 = -u1;
    }

   if  ((u+u1) > MAXLOG)
     return  DBL_MAX;

   /* u is exact, u1 is small.  */ 
  u = exp(u) * exp(u1);
   return (u);
}



};   // namespace Arithmetic 
};   // namespace Maths 
 #endif 
   //
