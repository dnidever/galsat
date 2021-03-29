pro plotlmc,file,dfile,over=over,twoplot=twoplot,save=save

; This program plots the nbody orbits

file = 'lmc.out'

; importing nbody data
arr = importnbody(file)
ns = n_elements(arr(*,0,0))
nb = n_elements(arr(0,*,0))
; R for LMC
r = dblarr(ns,nb)
for i=0,nb-1 do begin
  for j=0,ns-1 do r(j,i) = norm(reform(arr(j,i,2:4)))
end

; importing maketrail data
;marr = importdat('orbSGR.q9.dat')      ;t,x,y,z,vx,vy,vz,a,b,c
marr = importdat('orbLMC_121304.dat')      ;t,x,y,z,vx,vy,vz,a,b,c
nt = n_elements(marr(*,0))
mr = dblarr(nt)
for i=0,nt-1 do mr(i) = norm(reform(marr(i,1:3)))
gd = where(marr(*,0) le 0 and marr(*,0) ge -4.)

if not keyword_set(save) then erase  ; erase plot window

if keyword_set(save) then ps_open,'plotlmc',/color

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



; setting the ranges
dx=max(arr(*,*,2))-min(arr(*,*,2))
dy=max(arr(*,*,3))-min(arr(*,*,3))
dz=max(arr(*,*,4))-min(arr(*,*,4))
dr=max(r)-min(r)
dt=max(arr(*,*,0))-min(arr(*,*,0))
xr=[min(arr(*,*,2))-0.2*dx,max(arr(*,*,2))+0.2*dx]
yr=[min(arr(*,*,3))-0.2*dy,max(arr(*,*,3))+0.2*dy]
zr=[min(arr(*,*,4))-0.2*dz,max(arr(*,*,4))+0.2*dz]
rr=[min(r)-0.2*dr,max(r)+0.2*dr]
tr=[min(arr(*,*,0))-0.2*dt,max(arr(*,*,0))]

co1 = red
co2 = blue
thk = 0.5

; setting for postscript
if keyword_set(save) then begin
  loadct,13
  black=0
  purple=30
  blue=60
  aqua=80
  green=155   ;135
  yellow=200
  orange=225
  white=0
  red=300
  backg=white

  co1 = red
  co2 = blue

  ;sym8a = sym8a*0.8
  ;sym8b = sym8b*0.7
  ;thk = 3.0
  ;thk2 = 5.0
end

psym8

; Z vs X  (upper-left)
!p.multi=[4,2,2]
plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=xr,yr=zr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,2),arr(*,i,4),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,2)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
end
oplot,marr(gd,1),marr(gd,3),linestyle=2,co=co2

; R vs T (upper-right)
!p.multi=[3,2,2]
plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=tr,yr=rr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,0),r(*,i),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,0)],[r(0,i)],linestyle=0,co=co1,thick=thk,ps=8
end
oplot,marr(gd,0),mr(gd),linestyle=2,co=co2

; Y vs X (lower-left)
!p.multi=[2,2,2]
plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=xr,yr=yr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,2),arr(*,i,3),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,2)],[arr(0,i,3)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
end
oplot,marr(gd,1),marr(gd,2),linestyle=2,co=co2

; Z vs Y  (lower-right)
!p.multi=[1,2,2]
plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=yr,yr=zr,xs=1,ys=1
for i=0,nb-1 do begin
  oplot,arr(*,i,3),arr(*,i,4),linestyle=0,co=co1,thick=thk
  oplot,[arr(0,i,3)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
end
oplot,marr(gd,2),marr(gd,3),linestyle=2,co=co2

if keyword_set(save) then ps_close

; diagnostics
if keyword_set(dfile) then begin
  diag = importdiag(dfile)
  etot = reform(diag(1,*))
  ekin = reform(diag(2,*))
  epot = reform(diag(3,*))
endif

stop

end
