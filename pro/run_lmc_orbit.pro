pro run_lmc_orbit,psi,vperp,d=d

; I would like to be able to run maketrail by using
; three parameters for the velocity:
;  V  total velocity
;  and two angles, one in the plane of the sky (from B=90
;      meridian) and the other out of the plane of the sky.
;  This should make it easier to fit the orbit, I think.
;  Maybe also add "d" as a parameter.
;
; I think it would actually be better to use the radial
; velocity as one of the paramters, because that's fixed
; and then use the angle in the plane and the velocity
; in the plane.  Once the angle in the plane is fixed
; (and we have the right orbital plan) then there's
; only on parameter let (or possibly the distance).

deg2rad = !dpi/180.d
rad2deg = (180.d)/!dpi

psirad = psi*deg2rad

;l=0
;b=0
;psi = -45.
;vperp = 100.

; Hardwired values for LMC
;l = 293.
;b = -39.
;if not keyword_set(d) then d = 55.
l = 280.
b = -33.
if not keyword_set(d) then d = 49.
vgsr = 88.   ;84.  ;90.
tint = 15.
vcirc = 210.

; Getting position vector
lbd2xyz,l,b,d,x,y,z,/noprint
pvec = [x,y,z]  ; position vector

; Getting Vx, Vy, and Vz from Vgsr, Vperp and psi.

phi = l*deg2rad                ;phi in normal sph. coords.
th = 0.5*!dpi - b*deg2rad      ;theta in normal sph. coords.

rhat = [sin(th)*cos(phi),sin(th)*sin(phi),cos(th)]   ; r hat
that = [cos(th)*cos(phi),cos(th)*sin(phi),-sin(th)]  ; theta hat
phat = [-sin(phi),cos(phi),0.]                          ; phi hat

vr_vec = Vgsr*rhat     ;mag. times unit vector

; psi = 0 is along phi_hat, psi opens up towards north pole
; psi = 90 points towards the north pole
; psi = -90 points towards the sout pole
; psi = 180 points along -phi_hat
vperp_hat = cos(psirad)*phat - sin(psirad)*that  
vperp_vec = vperp*vperp_hat        ;mag. times unit vector

vtot_vec = vr_vec + vperp_vec      ;add up radial and perp components

; make strings
strtint = strtrim(tint,2)
strvcirc = strtrim(vcirc,2)
strxyz = strtrim(x,2)+',  '+strtrim(y,2)+',  '+strtrim(z,2)
strvxvyvz = strtrim(vtot_vec(0),2)+',  '+strtrim(vtot_vec(1),2)+',  '+strtrim(vtot_vec(2),2)

; Printing info
print,'LMC PARAMETERS'
print,'----------------------------------------'
print,'L     = ',strtrim(l,2),'   deg'
print,'B     = ',strtrim(b,2),'   deg'
print,'D     = ',strtrim(d,2),'   kpc'
print,'Vgsr  = ',strtrim(vgsr,2),'   km/s'
print,'Vperp = ',strtrim(vperp,2),'   km/s'
print,'Psi   = ',strtrim(psi,2),'   deg'
print,'----------------------------------------'
print,''

; Printing the param inputs
print,'PARAM INPUTS'
print,'-----------------------------------------------------'
print,'Time   = ',strtint,'   Gyrs'
print,'Vcirc  = ',strvcirc,'   km/s'
print,'XYZ    = ',strxyz,'   kpc'
print,'VxVyVz = ',strvxvyvz,'   km/s'
print,'-----------------------------------------------------'

;stop

; Using param file

; put output into param
openw,unit,/get_lun,'param'
printf,unit,'1.0e10'        ; mass    (not in current use)
printf,unit,'0.1'           ; fractional mass loss
printf,unit,'0.0'           ; mass of satellite which acts for dynamical friction
printf,unit,strtint         ; t_int
printf,unit,strvcirc        ; vcirc
printf,unit,strxyz             ; x,y,z
printf,unit,strvxvyvz          ; vx,vy,vz
close,unit
free_lun,unit

; Run maketrail
spawn,'maketrail',dum

; Run Ics2SM
spawn,'IcsToSM_new < icinput85',dum

;stop

end
