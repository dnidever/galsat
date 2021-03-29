function importplanet,file

; This function import the output data from an nbody run of planet
; The output is a 3D array: [Nsnaps, Nbodies, [t,mass,x,y,z,Vx,Vy,Vz]]
; NOTE:  Velocities will be returned in km/s.  They are automatically converted

; CONVERSION FACTORS
; 1 AU = 1.4960E8 kms, 1 yr = 3.155815E7 sec
; 1 AU/yr = 1.4960E8/3.155815E7 = 4.740455
kms2auyr = 2.10950D-1
auyr2kms = 4.740455D0

; 1 M_sun = 1.989E33 g, 1 M_jup = 1.8991E30 g
mjup2msun = 9.54799778D-4
msun2mjup = 1047.34D0


d = findfile(file)
if d eq '' then return,-1

spawn,'wc '+file,out
nlines = long(first_el(getwrd(out,0)))

; figuring out how many parameters there are
openr,unit,file,/get_lun
dum = ''
readf,unit,dum   ; # bodies
readf,unit,dum   ; time
readf,unit,dum   ; pars for 1st body
idum = strsplit(dum)
npar = n_elements(idum)
close,unit
free_lun,unit

; opening file
openr,unit,file,/get_lun

; some variables
n = 0
t = 0.0d
arr = dblarr(npar)
;arr = dblarr(7)
;arr = dblarr(13)
str = ''

; read number of bodies
readf,unit,n

; number of snaps
nsnap = nlines/(n+2.)

; creating array
;dum = {mass:0.d,x:0.d,y:0.d,z:0.d,vx:0.d,vy:0.d,vz:0.d}
;dum = {mass:0.d,x:dblarr(3),vel:dblarr(3)}
;dumn = replicate(dum,n)
;dum = {t:0.d,data:dblarr(7)}
;array = replicate(dum,nsnap)
;array = dblarr(nsnap,n,14)
;array = dblarr(nsnap,n,8)
array = dblarr(nsnap,n,npar+1)

; looping through the snaps
for i=0,nsnap-1 do begin
  if i gt 0 then readf,unit,str     ; n

  readf,unit,t

  ; looping through the bodies
  for j=0,n-1 do begin

    readf,unit,arr
    array(i,j,*) = [t,arr]

  end      ; for j

end     ; for i

close,unit
free_lun,unit



; CONVERTING VELOCITIES (AU/yr -> km/s)
array(*,*,5) = array(*,*,5)*auyr2kms
array(*,*,6) = array(*,*,6)*auyr2kms
array(*,*,7) = array(*,*,7)*auyr2kms

; CONVERTING FIRST BODY'S MASS (STAR; M_jup -> M_sun)
array(*,0,1) = array(*,0,1)*mjup2msun

;stop

return,array

end
