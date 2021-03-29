pro plotplanet,arr,dfile,over=over,twoplot=twoplot,stp=stp

; This program plots the nbody orbits

ns = n_elements(arr(*,0,0))
nb = n_elements(arr(0,*,0))
; R for LMC
r = dblarr(ns,nb)
for i=0,nb-1 do begin
  for j=0,ns-1 do r(j,i) = norm(reform(arr(j,i,2:4)))
end

window,0,xsize=650,ysize=600

; colors for my screen
red = 250
lred = 210
green = 190000
orange = 310000
yellow = 450000
blue = -25000
lblue = -15000
purple = -20000
white = -1

thk = 2      ; line thickness

;linestyle: 0-solid line, 1-dotted, 2-dashed, 3-dot dash, 4-dot dot
;          dot dash, 5-long dashes

psym8,1.5

co1 = green
co2 = red
co3 = orange
co4 = blue

; setting the ranges
dx=max(arr(*,*,2))-min(arr(*,*,2))
dy=max(arr(*,*,3))-min(arr(*,*,3))
dz=max(arr(*,*,4))-min(arr(*,*,4))
if dx eq 0.0 then dx=max([dx,dy,dz])
if dy eq 0.0 then dy=max([dx,dy,dz])
if dz eq 0.0 then dz=max([dx,dy,dz])
dr=max(r)-min(r)
dt=max(arr(*,*,0))-min(arr(*,*,0))
xr=[min(arr(*,*,2))-0.2*dx,max(arr(*,*,2))+0.2*dx]
yr=[min(arr(*,*,3))-0.2*dy,max(arr(*,*,3))+0.2*dy]
zr=[min(arr(*,*,4))-0.2*dz,max(arr(*,*,4))+0.2*dz]
rr=[min(r)-0.2*dr,max(r)+0.2*dr]
tr=[min(arr(*,*,0))-0.2*dt,max(arr(*,*,0))]

co1 = white
thk = 0.5

; Z vs X  (upper-left)
!p.multi=[4,2,2]
plot,dist(5),/nodata,xtit='X (AU)',ytit='Z (AU)',xr=xr,yr=zr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,2),arr(*,i,4),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,2)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
end

; R vs T (upper-right)
!p.multi=[3,2,2]
plot,dist(5),/nodata,xtit='t (yr)',ytit='R (AU)',xr=tr,yr=rr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,0),r(*,i),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,0)],[r(0,i)],linestyle=0,co=co1,thick=thk,ps=8        ; current position
end

; Y vs X (lower-left)
!p.multi=[2,2,2]
plot,dist(5),/nodata,xtit='X (AU)',ytit='Y (AU)',xr=xr,yr=yr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,2),arr(*,i,3),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,2)],[arr(0,i,3)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
end

; Z vs Y  (lower-right)
!p.multi=[1,2,2]
plot,dist(5),/nodata,xtit='Y (AU)',ytit='Z (AU)',xr=yr,yr=zr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,3),arr(*,i,4),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,3)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
end

; diagnostics
if keyword_set(dfile) then begin
  diag = importdiag(dfile)
  etot = reform(diag(1,*))
  ekin = reform(diag(2,*))
  epot = reform(diag(3,*))
  steps = reform(diag(4,*))
endif

if keyword_set(stp) then stop

end
