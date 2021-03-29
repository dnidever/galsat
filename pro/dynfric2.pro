pro dynfric2

    G = 0.000004300915
    ;G = 4.4967e-6
    mfric = 2e10
    vh2 = 121.*121.
    vel = [100.,100.,100.]
    v2 = norm(vel)^2.
    acc = vel*0.
    msat = 2e10
    ds = 13.
    r = 50.
    r3 = r^3.
    r2 = r^2.
    epsilon = 3.0
    PI = 3.1415926535

            ; Calculate mass enclosed within r (=(1/r^2 G)(d phi /dr))
            mr = 2.0*vh2*r3/(r2+ds*ds)/G 
            ; Calculate speed - note: strictly should be updated so speed coincides with position
            v = sqrt(v2);
      
            ; Coulomb logarithm via James' fitting function
            ot = 1.0/3.0;
            rtide = r*(ot*mfric/mr)^ot
            xlim = rtide
            coulog = alog(r/rtide)+jfit(xlim)+0.135
      
            ; coulog is the Coulomb Logarithm
            ; determined from Hashimoto 2003 
            ; note: coulog becomes 0 if satellite radius is less than 1.4 epsilon
            ; this is so we don't get dynamical acceleration at small radii
            ; epsilon is the plummer scale length of the satellite
            ;epsilon = soft
            if (r ge (1.4*epsilon)) then begin
              coulog_hashi = alog(r / (1.4 * epsilon))
            endif else begin 
              coulog_hashi = 0.
            endelse

            if (coulog gt 0.0) then begin
      
              ; Local density (=(1/4 pi r^2)(d mr/dr))
              ; rho = (vh2/2.0d0/pi/G)*(r*r+3.0d0*c*c)/(r*r+c*c)/(r*r+c*c)
               rho = (vh2/2.0/pi/G)*(r2+3.0*ds*ds)/(r2+ds*ds)/(r2+ds*ds)
      
               sigma2 = vh2
               xx = v/sqrt(2.0*sigma2)  ; capital X in B+T
      
              ; B+T, eqn 7-18, df = dv/dt = acc without the vector part
              df = 4.0*pi*coulog_hashi*(G*G)*mfric*rho*(erf(xx)-2.0*xx*exp(-(xx*xx))/sqrt(pi))/(v*v*v)
       
            endif else begin
              df = 0.0
            endelse

       
            ;dfx = -df*x(4)
            ;dfy = -df*x(5)
            ;dfz = -df*x(6)
      
            ; Calculating the Acceleration
            acc(0) = -df*vel(0)
            acc(1) = -df*vel(1)
            acc(2) = -df*vel(2)

print,acc

stop

end