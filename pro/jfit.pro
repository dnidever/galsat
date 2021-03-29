function jfit, x

;real jfitv, num, den, fx;
;real a,b,c,d,e,f,g,h,i,j;

;// Fitting parameters.

;     // PARAMETER (a=0.2021d0,b=0.04474d0,c=3.1621d0,
;     //& d=0.5488d0,e=2.9965d0,f=0.7375d0,g=0.9217d0,
;     //& h=0.2749d0,i=1.2583d0,j=2.3069d0)

      a=0.2021
      b=0.04474
      c=3.1621
      d=0.5488
      e=2.9965
      f=0.7375
      g=0.9217
      h=0.2749
      i=1.2583
      j=2.3069

      num = alog(1.0+a*x)*(b*( x^c )+d*( x^e ))
      den = ( 1.0 + f*( x^g ) + h*( x^i ))^j
      fx = alog(1.0+x)-x/(1.0+x)

      jfitv = num/den
      jfitv = jfitv/fx/fx

      return, jfitv
end