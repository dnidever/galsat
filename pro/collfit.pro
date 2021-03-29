pro collfit

;  This program tries to fit the model sgr orbit
;  after the collision to the "observed" orbit

; From plotcoll.pro
; sgr is the one I want, I think

;restore,'coll_orbits.dat'
restore,'orbSGR_stitch.dat'

; want phi increasing
t = reverse(reform(sgr(0,*)))
x = reverse(reform(sgr(1,*)))
y = reverse(reform(sgr(2,*)))
z = reverse(reform(sgr(3,*)))

; FROM GFIT.PRO

    ;weights = dblarr(npts)+1.
    ;weights = dblarr(npts)+0.01
    ftol = 1d-10
    ;fa = {x:x,y:y,err:weights}
    fa = {t:t,x:x,y:y,z:z}

    par = [1.457, -0.295, 282.]
    ;par = [ 1.4569997  ,   -0.33624576 ]
    ;par = [ 1.50  ,   -0.20 ]
    npar = n_elements(par)

    parinfo = replicate({limited:[0,0],limits:[0.,0.],step:0.,fixed:0},npar)
    parinfo(0).step = 0.01
    parinfo(1).step = 0.01
    parinfo(2).fixed = 1

    p = mpfit('collfunc', par, functargs=fa, perror=perror,niter=iter,status=status,$
               bestnorm=chisq, /quiet, parinfo=parinfo, dof=dof, autoderivative=1,$
               ftol=ftol)   ;,/nocatch)

    if status eq 0 then print,'BAD INPUT PARMETERS'
    if status eq 0 then stop

    if not keyword_set(perror) then stop

    ; scaled sigpars (from mpfit.pro)
    sigpar = perror * sqrt(chisq/dof)

    fpar = p   ;final parameters

    ; total resid
   ; result = gfunc(x,fpar)
   ; resid = y-result
   ; ;rms = sqrt(total(resid^2)/(npts-n_elements(fpar)-1))
   ; rms = sqrt(total(weights*resid^2.)/(npts-1))  ; using Haud's method, pg. 92, eg.6


stop

end
