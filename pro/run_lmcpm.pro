pro run_lmcpm,mualph,mudelt,vcirc,ds,dist=dist,save=save,tmin=tmin,$
               over=over,filename=filename,colrange=colrange,$
               nodots=nodots,nomin=nomin,title=title,rvel=rvel,$
               noplot=noplot,noprint=noprint,R0=R0,tfile=tfile,soft=soft

;+
; This program runs galsat with the LMC and Sgr and then plots the
; results.
;
; INPUTS:
;
;  mualph   proper motion in ra * cos(delta) in mas/yr
;  mudelt   proper motion in dec. in mas/yr
;  vcirc    circular velocity of the disk (default=220 km/s)
;  ds       halo softening parameter, originally 13 kpc
;  dist     distance to LMC
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
;  tfile    the temporary filename to use (WITHOUT any ending,
;               e.g. "tcoll" not "tcoll.out")
;  Created by D.Nidever Jan. 2005
;-

deg2rad = !dpi/180.d
rad2deg = (180.d)/!dpi

if not keyword_set(tmin) then tmin=4
strtmin = stringize(tmin,ndec=2)
if tmin lt 0. then strtmin = stringize(-tmin,ndec=2)

if n_elements(vcirc) eq 0 then vcirc=220.
if n_elements(ds) eq 0 then ds=13.

alpha = ten(05,24,00.0)*15.   ; 1950 ra
delta = ten(-69,48,00)        ; 1950 dec
if not keyword_set(rvel) then rvel = 262.2              ; v_helio, km/s
if not keyword_set(dist) then dist = 50.   ;49.             ; in kpc
if not keyword_set(soft) then soft = 1.0                ; LMC softening param in kpc

pm2vel,alpha,delta,rvel,mualph,mudelt,dist,vx,vy,vz,vcirc=vcirc

gl=280.5              ;current l,b
gb=-32.5
;glactc,alpha,delta,1950.,gl,gb,j,/fk4   ;  getting galactic coordinates
lbd2xyz,gl,gb,dist,x,y,z,R0=R0,/noprint

;if not keyword_set(vcirc) then vcirc = 210.

; From van der Marel
;mualph = 1.68  ; +/- 0.16 mas/yr
;mudelt = 0.34  ; +/- 0.16 mas/yr
;r = [-0.78, -41.55, -26.95]
;v = [-56, -219, 186]   ; +/- 36, 23, 35
;rvel = 262.2 ; +/- 3.4 km/s,  v_helio

strvxvyvz = strtrim(vx,2)+'  '+strtrim(vy,2)+'  '+strtrim(vz,2)
strxyz = strtrim(x,2)+'  '+strtrim(y,2)+'  '+strtrim(z,2)
strvcirc = strtrim(vcirc,2)
strds = strtrim(ds,2)
strsoft = strtrim(soft,2)

; Printing info
if not keyword_set(noprint) then begin
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
  print,'Dsoft  = ',strds,'   kpc'
  print,'XYZ    = ',strxyz,'   kpc'
  print,'VxVyVz = ',strvxvyvz,'   km/s'
  print,'-----------------------------------------------------'
endif ; not noprint

;stop

; Using param file

if n_elements(tfile) eq 0 then tfile = 'tcoll'

; put output into param
openw,unit,/get_lun,tfile+'.in'
printf,unit,'1'
printf,unit,'0'
printf,unit,strvcirc
printf,unit,strds
printf,unit,'2.0e10 ',strsoft,' ',strxyz,' ',strvxvyvz
;printf,unit,'7.5e8 0.0 16.1556 2.27043 -5.88738  237.832 -42.3568 221.998'
;printf,unit,'2.0e10 ',strxyz,' ',strvxvyvz
;printf,unit,'7.5e8 16.1556 2.27043 -5.88738  237.832 -42.3568 221.998'
close,unit
free_lun,unit

; Running nbody
print,' '
print,'Running GALSAT program ... 3 sec.'
spawn,'rm '+tfile+'.out'
;spawn,'galsat2 -d 0.01 -o 0.01 -t 4 -i -e 4 < tcoll.in > tcoll.out',dum
;spawn,'galsat2 -d 0.01 -o 0.01 -t '+strtmin+' -i -e '+strtmin+' < tcoll.in > tcoll.out',dum
;spawn,'galsat3 -d 0.01 -o 0.01 -t '+strtmin+' -i -e '+strtmin+' < '+tfile+'.in > '+tfile+'.out',dum
;spawn,'galsat4 -d 0.01 -o 0.01 -t '+strtmin+' -i -e '+strtmin+' < '+tfile+'.in > '+tfile+'.out',dum
spawn,'galsat5 -d 0.01 -o 0.01 -t '+strtmin+' -i -e '+strtmin+' < '+tfile+'.in > '+tfile+'.out',dum

; calculating the minimum separation
;arr = importnbody('tcoll.out')
;ns = n_elements(arr(*,0,0))
;minsep = 999999.
;minl = 999999.
;mint = 999999.
;t = reform(arr(*,0,0))
;sep = fltarr(ns)
;
;for l=0,ns-1 do sep(l) = norm(reform(arr(l,0,2:4)-arr(l,1,2:4)))
;tbrack = where(t gt -2.3 and t lt -1.0)
;colind = [minloc(sep(tbrack),/first)]+min(tbrack)
;minsep = sep(colind)
;mint = t(colind)

;stop

; Running plotting program
;print,' '
;print,'Plotting ... 10 sec.'
if not keyword_set(noplot) then $
plotcoll,tfile+'.out',save=save,tmin=tmin,mualpha=mualph,mudelta=mudelt,vcirc=vcirc,$
          over=over,filename=filename,colrange=colrange,nodots=nodots,nomin=nomin,$
          title=title,lmcdist=dist,ds=ds

;; Run maketrail
;spawn,'maketrail',dum
;
;; Run Ics2SM
;spawn,'IcsToSM_new < icinput85',dum

;stop

end
