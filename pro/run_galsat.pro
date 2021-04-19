pro run_galsat,output,ra=ra,dec=dec,year=year,lon=lon,lat=lat,dist=dist,xgc=xgc,ygc=ygc,zgc=zgc,$
               mualph=mualph,mudelt=mudelt,rvel=rvel,vx=vx,vy=vy,vz=vz,mass=mass,$
               vcirc=vcirc,ds=ds,save=save,tmin=tmin,filename=filename,colrange=colrange,$
               noprint=noprint,R0=R0,tfile=tfile,soft=soft,df=df,prog=prog,$
               step=step,dstep=dstep,ostep=ostep,forward=forward,mw=mw,minstep=minstep,$
               maxstep=maxstep,fstep=fstep,xtra=xtra,stp=stp,noplot=noplot

; This program runs galsat for a MW satellite
; There are various options for giving the input parameters.
; LOCATION:  (1) ra, dec, dist
;            (2) lon, lat, dist
;            (3) xgc, ygc, zgc
; VELOCITY:  (1) mualph, mudelt, rvel
;            (2) vx, vy, vz
;
; INPUTS:
;  =mualph   proper motion in ra * cos(delta) in mas/yr
;  =mudelt   proper motion in dec. in mas/yr
;  =rvel     Radial velocity of the satellite (heliocentric) in km/s
;  =vx       Velocity in X_galactocentric (towards gal.center) in km/s
;  =vy       Velocity in Y_galactocentric (towards lon=90) in km/s     
;  =vz       Velocity in Z_galactocentric (towards NGP) in km/s
;  =ra       Right Ascension in hours.
;  =dec      Declination in degrees.
;  =year     Year of ra,dec coordinates (default year=2000.)
;  =lon      Galactic Longitude in degrees.
;  =lat      Galactic Latitude in degrees.
;  =dist     Distance from sun in kpc.
;  =xgc      X_galactocentric in kpc.
;  =ygc      Y_galactocentric in kpc.
;  =zgc      Z_galactocentric in kpc.
;  =mass     Mass of satellite in M_sun.
;  =vcirc    circular velocity of the disk (default=220 km/s)
;  =ds       halo softening parameter, originally 13 kpc
;  /df       Use dynamical friction. df=1 -> yes, df=0 -> no. (default df=1)
;  =step     GALSAT: Step size control parameter (default: step = 0.01)
;  =dstep    GALSAT: Diagnostic interval (default: dstep = tmin)
;  =ostep    GALSAT: Output interval (default: ostep = 0.01)
;  =minstep  GALSAT: Minimum timestep (positive), when using test particles
;  =maxstep  GALSAT: Maximum timestep (positive), when using test particles
;  =fstep    GALSAT: Fixed timestep (positive), when using test particles
;  /xtra     GALSAT: Output extra debugging information (default: xtra=0)
;  /forward  Integrate forward instead of backwards (default: forward=0)
;  /mw       Do NOT use the Milky Way potential (default: mw=0)
;  /df       USE dynamical friction (default df=1)
;  prog      Galsat program name to use (default prog='galsat6')
;  R0        The distance from the sun to the galactic center (default 8.5 kpc)
;  =tmin     How far back to run the simulation (default=4 Gyrs)
;            This number should be positive (but can handle negative).
;  /noprint  Don't print anything to the screen
;  /noplot   Don't plot the orbit
;  /stp      Stop at the end of the program
;
; OUTPUT:
;  output   Output array of the galsat results
;
; PROGRAMS USED:
;  pm2vel   Convert proper motions to xyz galactocentric velocities
;  glactc   Convert b/w ra/dec and lon/lat
;  lbd2xyz  Convert lon/lat/dist to xyz galactocentric distances.
;
;  Created by D.Nidever Jan. 2005

deg2rad = !dpi/180.d
rad2deg = (180.d)/!dpi

if n_elements(tmin) eq 0 then tmin=4
strtmin = stringize(tmin,ndec=2)
if tmin lt 0. then strtmin = stringize(-tmin,ndec=2)

; Setting MW grav. potential parameters
if n_elements(vcirc) eq 0 then vcirc=220.
if n_elements(ds) eq 0 then ds=13.

; OTHER PARMETERS
if n_elements(mass) eq 0 then mass = 1e8
if n_elements(soft) eq 0 then soft = 1.0
if n_elements(df) eq 0 then df = 1

; LOCATION
nra = n_elements(ra)
ndec = n_elements(dec)
nlon = n_elements(lon)
nlat = n_elements(lat)
nxgc = n_elements(xgc)
nygc = n_elements(ygc)
nzgc = n_elements(zgc)
ndist = n_elements(dist)

lcase = 0
if nra ne 0 and ndec ne 0 and ndist ne 0 then lcase=1   ;ra/dec
if nlon ne 0 and nlat ne 0 and ndist ne 0 then lcase=2  ;lon/lat
if nxgc ne 0 and nygc ne 0 and nzgc ne 0 then lcase=3  ;x/y/z

case lcase of

  ; RA/DEC
  1: begin

    if n_elements(year) eq 0 then year=2000.    

    ; Converting ra,dec to lon/lat
    glactc,ra,dec,year,lon,lat,1

    ; Converting lon,lat,dist to xyz
    lbd2xyz,lon,lat,dist,xgc,ygc,zgc,R0=R0,/noprint

  end 

  ; LON/LAT
  2: begin

    if n_elements(year) eq 0 then year=2000.    

    ; Converting lon,lat,dist to xyz
    lbd2xyz,lon,lat,dist,xgc,ygc,zgc,R0=R0,/noprint

    ; Converting lon,lat to ra/dec for velocity use
    glactc,ra,dec,year,lon,lat,2

  end

  ; X/Y/Z
  3: begin

    if n_elements(year) eq 0 then year=2000.

    ; Converting to ra/dec for velocity use
    xyz2lbd,xgc,ygc,zgc,lon,lat,dist,R0=R0,/noprint

    ; Converting lon,lat to ra/dec for velocity use
    glactc,ra,dec,year,lon,lat,2

  end

  ; NONE OF THE ABOVE
  else: begin
    print,'NOT ENOUGH LOCATION INFORMATION GIVEN'
    print,' LOCATION:  (1) ra, dec, dist'
    print,'            (2) lon, lat, dist'
    print,'            (3) xgc, ygc, zgc'
    return
  end

endcase


; VELOCITY
nmua = n_elements(mualph)
nmud = n_elements(mudelt)
nrv = n_elements(rvel)
nvx = n_elements(vx)
nvy = n_elements(vy)
nvz = n_elements(vz)

vcase = 0
if nmua ne 0 and nmud ne 0 and nrv ne 0 then vcase=1   ; mualph,mudelt
if nvx ne 0 and nvy ne 0 and nvz ne 0 then vcase=2   ; vx,vy,vz

case vcase of

  ; MUALPH, MUDELT, RVEL
  1: begin

    ; Precessing ra, dec.  They need to be in degrees B1950
    bprecess,ra*15.d,dec,ra1950,dec1950,epoch=year

    ; Getting galactocentric velocities, ra/dec need to be in dec degrees B1950.
    pm2vel,ra1950,dec1950,rvel,mualph,mudelt,dist,vx,vy,vz,vcirc=vcirc

  end

  ; VX, VY, VZ
  2: begin

    ; Don't need to do anything

  end

  ; NONE OF THE ABOVE
  else: begin
    print,'NOT ENOUGH VELOCITY INFORMATION GIVEN'
    print,' VELOCITY:  (1) mualph, mudelt, rvel;
    print,'            (2) vx, vy, vz;
    return
  end

endcase

; GETTING STRINGS
strvxvyvz = strtrim(vx,2)+'  '+strtrim(vy,2)+'  '+strtrim(vz,2)
strxyz = strtrim(xgc,2)+'  '+strtrim(ygc,2)+'  '+strtrim(zgc,2)
strvcirc = strtrim(vcirc,2)
strds = strtrim(ds,2)
strsoft = strtrim(soft,2)
strmmass = strtrim(mass,2)

; PRINTING INFO
if not keyword_set(noprint) then begin
  ; Printing the param inputs
  print,''
  print,'PARAM INPUTS'
  print,'-----------------------------------------------------'
  print,'Tmin   = ',strtmin,'   Gyrs'
  print,'Vcirc  = ',strvcirc,'   km/s'
  print,'Dsoft  = ',strds,'   kpc'
  print,'XYZ    = ',strxyz,'   kpc'
  print,'VxVyVz = ',strvxvyvz,'   km/s'
  print,'Soft   = ',strsoft,'   kpc'
  print,'Mass   = ',strmmass,' Msun'
  print,'-----------------------------------------------------'
endif ; not noprint

; RUNNING GALSAT.PRO
input = fltarr(1,8)
input(0,*) = [mass,soft,xgc,ygc,zgc,vx,vy,vz]
if n_elements(step) eq 0 then step=0.01
if n_elements(ostep) eq 0 then ostep=0.01
if n_elements(dstep) eq 0 then dstep=tmin

galsat,input,output,step=step,ostep=ostep,dstep=tmin,tmin=tmin,df=df,prog=prog,mw=mw,$
       forward=forward,minstep=minstep,maxstep=maxstep,fstep=fstep,xtra=xtra,vcirc=vcirc

; PLOTTING
if not keyword_set(noplot) then plotnbody,arr=output

if keyword_set(stp) then stop

end
