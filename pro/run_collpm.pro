pro run_collpm,mualph,mudelt,output,vcirc=vcirc,dhalo=dhalo,dist=dist,save=save,$
               tmin=tmin,over=over,filename=filename,colrange=colrange,$
               nodots=nodots,nomin=nomin,title=title,rvel=rvel,$
               noplot=noplot,noprint=noprint,R0=R0,tfile=tfile,soft=soft,$
               df=df,prog=prog,mlmc=mlmc,msgr=msgr,stp=stp

; This program runs galsat with the LMC and Sgr and then plots the
; results.
;
; INPUTS:
;  mualph   proper motion in ra * cos(delta) in mas/yr
;  mudelt   proper motion in dec. in mas/yr
;  vcirc    circular velocity of the disk (default=220 km/s)
;  dhalo    halo softening parameter, originally 13 kpc
;  dist     distance to LMC
;  df       Use dynamical friction. df=1 -> yes, df=0 -> no. (default df=1)
;  prog     Galsat program name to use (default prog='galsat6')
;  R0       the distance from the sun to the galactic center (8.5 kpc)
;  tmin     how far back to run the simulation (default=4 Gyrs)
;           This number should be positive (but can handle negative).
;  /save    save the plot to postscript
;  filename name of postscript file (without the .ps)
;  colrange time range for minmum separation
;  /over    overplot stiched Sgr, solitary Sgr, and old LMC orbits.
;  /nodots  don't overplot the equally spaced points
;  /nomin   don't overplot the minimum separation points
;  /noplot  don't plot the results
;  /noprint  don't print anything to the screen
;  title    overall title of the plot
;  mlmc     LMC mass
;  msgr     SGR mass
;
; OUTPUT:
;  output   Output array of the galsat results
;
;  Created by D.Nidever Jan. 2005

deg2rad = !dpi/180.d
rad2deg = (180.d)/!dpi

if not keyword_set(tmin) then tmin=4
strtmin = stringize(tmin,ndec=2)
if tmin lt 0. then strtmin = stringize(-tmin,ndec=2)

; Setting MW grav. potential parameters
if n_elements(vcirc) eq 0 then vcirc=220.
if n_elements(dhalo) eq 0 then dhalo=13.

; Getting RA and DEC
alpha = ten(05,24,00.0)*15.   ; 1950 ra
delta = ten(-69,48,00)        ; 1950 dec
if not keyword_set(rvel) then rvel = 262.2              ; v_helio, km/s
if not keyword_set(dist) then dist = 50.   ;49.             ; in kpc
if not keyword_set(soft) then soft = 1.0                ; LMC softening param in kpc
if not keyword_set(mlmc) then mlmc = 2e10
if not keyword_set(msgr) then msgr = 7.5e8

; Getting galactocentric velocities
pm2vel,alpha,delta,rvel,mualph,mudelt,dist,vx,vy,vz,vcirc=vcirc

; Getting galactocentric x,y,z coordinates
gl=280.5              ;current l,b
gb=-32.5
;glactc,alpha,delta,1950.,gl,gb,j,/fk4
lbd2xyz,gl,gb,dist,x,y,z,R0=R0,/noprint

; From van der Marel
;mualph = 1.68  ; +/- 0.16 mas/yr
;mudelt = 0.34  ; +/- 0.16 mas/yr
;r = [-0.78, -41.55, -26.95]
;v = [-56, -219, 186]   ; +/- 36, 23, 35
;rvel = 262.2 ; +/- 3.4 km/s,  v_helio

strvxvyvz = strtrim(vx,2)+'  '+strtrim(vy,2)+'  '+strtrim(vz,2)
strxyz = strtrim(x,2)+'  '+strtrim(y,2)+'  '+strtrim(z,2)
strvcirc = strtrim(vcirc,2)
strdhalo = strtrim(dhalo,2)
strsoft = strtrim(soft,2)
strmlmc = strtrim(mlmc,2)
strmsgr = strtrim(msgr,2)

; Printing info
if not keyword_set(noprint) then begin
  print,''
  print,'LMC PARAMETERS'
  print,'----------------------------------------'
  print,'L     = ',strtrim(gl,2),'   deg'
  print,'B     = ',strtrim(gb,2),'   deg'
  print,'D     = ',strtrim(dist,2),'   kpc'
  print,'Vhel  = ',strtrim(rvel,2),'   km/s'
  print,'mualph = ',strtrim(mualph,2),'  mas/yr'
  print,'mudelt = ',strtrim(mudelt,2),'  mas/yr'
  print,'----------------------------------------'
  print,''

  ; Printing the param inputs
  print,'PARAM INPUTS'
  print,'-----------------------------------------------------'
  print,'Tmin   = ',strtmin,'   Gyrs'
  print,'Vcirc  = ',strvcirc,'   km/s'
  print,'Dhalo  = ',strdhalo,'   kpc'
  print,'XYZ    = ',strxyz,'   kpc'
  print,'VxVyVz = ',strvxvyvz,'   km/s'
  print,'Soft   = ',strsoft,'   kpc'
  print,'Mlmc   = ',strmlmc,' Msun'
  print,'Msgr   = ',strmsgr,' Msun'
  print,'-----------------------------------------------------'
endif ; not noprint

;stop

; Running galsat.pro
input = fltarr(2,8)
input(0,*) = [mlmc,soft,x,y,z,vx,vy,vz]
input(1,*) = [msgr, 0.0, 16.1556, 2.27043,-5.88738, 237.832, -42.3568, 221.998]
galsat,input,output,step=0.01,ostep=0.01,dstep=tmin,tmin=tmin,df=df,prog=prog,dhalo=dhalo

; Running plotting program
;print,' '
;print,'Plotting ... 10 sec.'
if not keyword_set(noplot) then $
plotcoll,save=save,tmin=tmin,mualpha=mualph,mudelta=mudelt,vcirc=vcirc,$
          over=over,filename=filename,colrange=colrange,nodots=nodots,nomin=nomin,$
          title=title,lmcdist=dist,dhalo=dhalo,arr=output,stp=stp
;stop

end
