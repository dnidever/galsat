pro dynfric1

    PI = 3.1415926535		; Pi
    G = 0.000004300915  	; Grav const in sim units
    G2 = G*G			; Grav const squared

    v_halo2 = 121.*121.
    vsat = [100.,100.,100.]
    acc = vsat*0.
    msat = 2e10
    rsat = 50.
    epsilon = 3.0

    ; the square of the mag of the satellite velocity
    vsat2 = vsat(0)*vsat(0)+vsat(1)*vsat(1)+vsat(2)*vsat(2)

    ; this is the definition of sigma for a logarithmic potential
    sig = sqrt(v_halo2)

    ; Chi is ...
    ; Defined in both Zentner and KVJ code.
    chi = sqrt(vsat2)/(sqrt(2.0)*sig)	
    
    ; logarithmic potential-specific mass_int_r
    mass_int_r = 2 * v_halo2 * rsat / G

    ; logarithmic potential-specific rho
    rho = 2 * v_halo2 / (4. * PI * G * rsat*rsat)

    ; rtide is tidal radius of satellite? tidal radius of halo?
    ; ported from kvj code.
    rtide = rsat*(msat/(3.0*mass_int_r))^(1.0/3.0)

    ; coulog is the Coulomb Logarithm
    ; determined from Hashimoto 2003 
    ; note: coulog becomes 0 if satellite radius is less than 1.4 epsilon
    ; this is so we don't get dynamical acceleration at small radii
    if (rsat ge (1.4*epsilon)) then $
    	coulog_hashi = alog(rsat / (1.4 * epsilon)) $
    else $
	coulog_hashi = 0.

    ; force_dynfric is the force on the satellite due to dynamical friction.
    ; Same equation found in KVJ, Zentner, and Binney & Merrifield.
    ; Originally derived from Chandrasekhar.
    force_dynfric = 4*PI*coulog_hashi*G2*msat*rho*( erf(chi)-2*chi*exp(-(chi*chi))/sqrt(PI))/(vsat2*sqrt(vsat2))

    ; the true dynamical friction force opposes the motion of the body
    ; by directly opposing the velocity, not the acceleration.
    ; Acceleration = force divided by mass of satellite.
    for k = 0,2 do $
	acc(k) = -force_dynfric*vsat(k)

    print,acc

    stop

end
