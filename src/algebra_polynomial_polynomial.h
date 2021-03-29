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

#ifndef _CC_ALGEBRA_POLYNOMIAL_POLYNOMIAL_H_ 
   #define _CC_ALGEBRA_POLYNOMIAL_POLYNOMIAL_H_ 

 #ifndef cc_error 
   #include <assert.h> 
 #ifndef NDEBUG 
   #include <iostream> 
 #define cc_error(a) std::cerr << (a) << std::endl; 
   #else 
   #define cc_error(a) 
   #endif 
   #endif 

 namespace Algebra
{

 namespace Polynomial
{

 //

double  polyEval( double  x,  const   double  coef[],  int  N)
{
   const   double  *p = coef;
   double  ans = *p++;
   int  i = N;

   do 
    ans = ans * x  +  *p++;
   while ( --i );

   return ( ans );
}


double  polyEval1( double  x,  const   double  coef[],  int  N)
{
   const   double  *p = coef;
   double  ans = x + *p++;
   int  i = N-1;

   do 
    ans = ans * x  + *p++;
   while ( --i );

   return ( ans );
}



};   // namespace Poly 
};   // namespace Maths 

 #endif 
   //
