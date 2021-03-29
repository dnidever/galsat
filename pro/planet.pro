pro planet,input,output,step=step,dstep=dstep,ostep=ostep,tmax=tmax,$
           noprint=noprint,xtra=xtra

; This program runs the c++ planet program.  It's probably best to
; call it from a batch file. 
;
; NOTE: All velocities should be input in KM/S.  They will also be
;       output this way.  The velocity internal to the C++ planet
;       program is AU/yr, but planet.pro and importplanet.pro 
;       automatically makes the right conversions.
;
; NOTE: The first nbody is assumed to be the STAR and have mass in 
;       units of M_sun, and all other bodies to be planets have masses
;       in units of M_jup.  All appropriate conversions will be made
;       internally.  1 M_sun = 1047.34 M_jup, 1 M_jup = 9.5480E-4 M_sun
;
; NOTE: All distances should be input in units of AU.
;
; INPUT
;   input   Array of inputs [Nbodies, [mass,x,y,z,Vx,Vy,Vz]]
;   step    Step size control parameter (default: step = 0.01)
;   dstep   Diagnostic interval (default: dstep = tmin)
;   ostep   Output interval (default: ostep = 0.01)
;   tmax    Maximum time to integral to (default tmin = 4.0)
;   /xtra   Output extra debugging information
;
; OUTPUT
;   output  Array of galsat output. [Nsnap,Nbodies,[t,mass,x,y,z,Vx,Vy,Vz]]
;
; PROGRAMS CALLED:
;  importplanet   Imports the output of planet
;  maketemp       Makes temporary filename
;
; Created by David Nidever August 2005

; Bad Input Parameters
if n_params() eq 0 or n_elements(input) eq 0 then begin
  print,'syntax - galsat,input,output,step=step,dstep=dstep,ostep=ostep,tmax=tmax'
  print,'                noprint=noprint,xtra=xtra'
endif

; CONVERSION FACTORS
; 1 AU = 1.4960E8 kms, 1 yr = 3.155815E7 sec
; 1 AU/yr = 1.4960E8/3.155815E7 = 4.740455
kms2auyr = 2.10950D-1
auyr2kms = 4.740455D0

; 1 M_sun = 1.989E33 g, 1 M_jup = 1.8991E30 g
mjup2msun = 9.54799778D-4
msun2mjup = 1047.34D0


prog='planet'
if n_elements(tmax) eq 0 then tmax=4.0
if n_elements(step) eq 0 then step=0.01
if n_elements(dstep) eq 0 then dstep=tmax
if n_elements(ostep) eq 0 then ostep=0.01

; # of Bodies
sz = size(input)
case sz(0) of
  1: nbody = 1
  2: nbody = sz(1)
endcase

; Getting inputs ready
if tmax lt 0. then strtmax = stringize(-tmax,ndec=2) else strtmax=strtrim(tmax,2)
strstep = strtrim(step,2)
strdstep = strtrim(dstep,2)
strostep = strtrim(ostep,2)
strnbody = strtrim(nbody,2)

; Setting flags
flags = ''
if keyword_set(xtra) then flags = flags+' -x'

; CONVERTING VELOCITIES (km/s -> AU/yr)
input2 = input
input2(*,4) = input2(*,4)*kms2auyr
input2(*,5) = input2(*,5)*kms2auyr
input2(*,6) = input2(*,6)*kms2auyr

; CONVERTING FIRST BODY'S MASS (STAR; M_sun -> M_jup
input2(0,0) = input2(0,0)*msun2mjup

; CREATING INPUT FILE
tfile = maketemp('tplan')
openw,unit,/get_lun,tfile+'.in'
printf,unit,strnbody
printf,unit,'0'
if nbody ge 1 then $
    for i=0,nbody-1 do printf,unit,strtrim(reform(input2(i,*)),2)
if nbody eq 1 then printf,unit,strtrim(input2,2)
close,unit
free_lun,unit

; Printing info
if not keyword_set(noprint) then begin
  ; Printing the param inputs
  print,'PARAM INPUTS'
  print,'--------------------------'
  print,'Prog = ',prog
  print,'Nbodies = ',strnbody
  print,'Tmax   = ',strtmax,'   Gyrs'
  print,'Step   = ',strstep
  print,'Ostep  = ',strostep,' Gyrs'
  print,'Dstep  = ',strdstep,' Gyrs'
  print,'--------------------------'

  print,''
  print,'MASS (M_sun/jup)  X    Y    Z  (AU)   Vx    Vy    Vz (km/s)'
  print,'--------------------------------------------------------'

  for i=0,nbody-1 do begin
    out = strjoin(strmid(strtrim(reform(input(i,*)),2),0,7),' ')
    print,out
  end


  print,'--------------------------------------------------------'
endif ; not noprint

; RUNNING PLANET
print,' '
print,'Running PLANET program ... a few sec.'
cmd = '( ./'+prog+' -d '+strstep+' -o '+strostep+' -t '+strtmax+' -i -e '+strdstep
cmd = cmd + flags+' < '+tfile+'.in > '+tfile+'.out ) > & '+tfile+'.diag'
spawn,cmd,dum

; GETTING THE OUTPUT
output = importplanet(tfile+'.out')

; REMOVING TEMPORARY FILES
spawn,'\rm '+tfile+'.in'
spawn,'\rm '+tfile+'.out'
spawn,'\rm '+tfile+'.diag'

;stop

end
