function collfunc,par,t=t,X=x,Y=y,Z=z,noderiv=noderiv,rms=rms,arr=arr,$
                  df=df,soft=soft,prog=prog,mlmc=mlmc,msgr=msgr,chi=chi

; This program is run by collgrid.pro
;
; INPUT
;   par      The input parameters, [mualph,mudelt,vcirc,ds]
;   t        Array of times for the Sgr orbit
;   x        Array of x-values for the Sgr orbit
;   y        Array of y-values for the Sgr orbit
;   z        Array of z-values for the Sgr orbit
;   /noderiv Do not compute the derivative
;   rms      The RMS fit to the pre-collision Sgr orbit
;   chi      The reduced chi squared fit to the pre-collision Sgr orbit
;   soft     The LMC Plummer potential scale length.
;   df       Use dynamical friction. df=1 -> yes, df=0 -> no. (default df=1)
;   prog     Galsat program name to use (default prog='galsat6')
;   mlmc     LMC mass
;   msgr     SGR mass
;
; OUTPUT
;   dev  Deviates
;
; x,y,z are the galactocentric coordinates of the "observed" orbit
;
; Created by David Nidever 2005

; Getting the parameters
mualph = par(0)
mudelt = par(1)
vcirc = par(2)
ds = par(3)
tmin = 5.0
;rvel = par(2)

; Running RUN_COLLPM.PRO
run_collpm,mualph,mudelt,arr,vcirc=vcirc,ds=ds,tmin=tmin,rvel=rvel,/noplot,soft=soft,$
           df=df,prog=prog,mlmc=mlmc,msgr=msgr

; If there there are softening parameters then get rid of them.
if (n_elements(arr(0,0,*)) eq 9) then begin
 dum = arr
 arr = arr(*,*,0:7)*0.
 arr(*,*,0:1) = dum(*,*,0:1)
 arr(*,*,2:7) = dum(*,*,3:8)
 dum = 0.
endif

; Creating the arrays
tm = reform(arr(*,1,0))
xm = reform(arr(*,1,2))
ym = reform(arr(*,1,3))
zm = reform(arr(*,1,4))
nm = n_elements(arr(*,1,0))

; Getting the longitudes
xyz2sph,xm,ym,zm,rm,thm,phim    ; model
xyz2sph,x,y,z,r,th,phi   ; observed
phi = phi-360.

; Getting the good points
;ind = where(phi lt -450. and phi gt -1000,nind)
ind = where(t lt -2.5 and t gt -4.5,nind)

; If phi acts normally then continue, otherwise bomb
if max(phim) lt 4000. then begin
 
  ; splining the model phi
  nxm = spline(phim,xm,phi(ind))
  nym = spline(phim,ym,phi(ind))
  nzm = spline(phim,zm,phi(ind))
 
  ; deviants
  rms = sqrt( total( (nxm-x(ind))^2. + (nym-y(ind))^2. + (nzm-z(ind))^2. )/(nind-1) )
  chi = total( ( (nxm-x(ind))^2. + (nym-y(ind))^2. + (nzm-z(ind))^2. )/sqrt(nxm^2.+nym^2.+nzm^2.) )/nm
  dev = sqrt( (nxm-x(ind))^2. + (nym-y(ind))^2. + (nzm-z(ind))^2. )

endif else begin

 dev = phi(ind)*0.+999999.
 rms = 999999.
 chi = 999999.

 ;stop

endelse

stop

return,dev

end
