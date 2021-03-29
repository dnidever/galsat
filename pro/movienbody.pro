pro movienbody,file,dfile,over=over,twoplot=twoplot,nomw=nomw

; This program plots the nbody orbits

; importing nbody data
arr = importnbody(file)

; if there there are softening parameters then get rid of them.
if (n_elements(arr(0,0,*)) eq 9) then begin
  dum = arr
  arr = arr(*,*,0:7)*0.
  arr(*,*,0:1) = dum(*,*,0:1)
  arr(*,*,2:7) = dum(*,*,3:8)
  dum = 0.
endif

ns = n_elements(arr(*,0,0))
nb = n_elements(arr(0,*,0))
; R for LMC
r = dblarr(ns,nb)
for i=0,nb-1 do begin
  for j=0,ns-1 do r(j,i) = norm(reform(arr(j,i,2:4)))
end

erase  ; erase plot window

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
;
;psym8,1.5

co1 = green
co2 = red
co3 = orange
co4 = blue

; setting the ranges
;dx=max(arr(*,*,2))-min(arr(*,*,2))
;dy=max(arr(*,*,3))-min(arr(*,*,3))
;dz=max(arr(*,*,4))-min(arr(*,*,4))
;dr=max(r)-min(r)
;dt=max(arr(*,*,0))-min(arr(*,*,0))
;xr=[min(arr(*,*,2))-0.2*dx,max(arr(*,*,2))+0.2*dx]
;yr=[min(arr(*,*,3))-0.2*dy,max(arr(*,*,3))+0.2*dy]
;zr=[min(arr(*,*,4))-0.2*dz,max(arr(*,*,4))+0.2*dz]
;rr=[min(r)-0.2*dr,max(r)+0.2*dr]
;tr=[min(arr(*,*,0))-0.2*dt,max(arr(*,*,0))]

co1 = white
thk = 0.5

size = 2.*stdev(arr(*,*,2:4))
;size = 20.
xr = [-1.,1.]*size
yr = [-1.,1.]*size
zr = [-1.,1.]*size

!p.multi=0

; creating the disk
rdisk = 16.   ;25
diskx = rdisk*sin(findgen(1000)/999.*2.*!dpi)
disky = rdisk*cos(findgen(1000)/999.*2.*!dpi)
diskz = disky*0.
bd = where(abs(diskx) le 5. or abs(disky) le 5.,nbd)
remove,bd,diskx,disky,diskz

; creating the bulge
sphere,2.,bx,by,bz,/upper

; sphere
sphere,0.1,sphx,sphy,sphz

psym8,1.

for i=0,ns-1 do begin

  snap = reform(arr(i,*,*))
  xarr = (snap(*,2))(*)
  yarr = (snap(*,3))(*)
  zarr = (snap(*,4))(*)

  daz = 360./ns
  az0 = 360.  ;330.  ; 275.  ;-45.
  ax =  20.   ;45
  az = daz*i+az0   ;45.
  charsize = 1d-10   ; 2.5
  surface,dist(10),/nodata,az=az,ax=ax,/save,xr=xr,yr=yr,zr=zr,$
     xtit='X',ytit='Y',ztit='Z',charsize=2.5,xs=4,ys=4,zs=4

  ; x-axes
  plots,[-1,1]*size,[-1,-1]*size,[-1,-1]*size,/t3d
  plots,[-1,1]*size,[1,1]*size,[-1,-1]*size,/t3d
  plots,[-1,1]*size,[-1,-1]*size,[1,1]*size,/t3d
  plots,[-1,1]*size,[1,1]*size,[1,1]*size,/t3d

  ; y-axes
  plots,[-1,-1]*size,[-1,1]*size,[-1,-1]*size,/t3d
  plots,[1,1]*size,[-1,1]*size,[-1,-1]*size,/t3d
  plots,[-1,-1]*size,[-1,1]*size,[1,1]*size,/t3d
  plots,[1,1]*size,[-1,1]*size,[1,1]*size,/t3d

  ; z-axes
  plots,[1,1]*size,[1,1]*size,[-1,1]*size,/t3d
  plots,[-1,-1]*size,[1,1]*size,[-1,1]*size,/t3d
  plots,[1,1]*size,[-1,-1]*size,[-1,1]*size,/t3d
  plots,[-1,-1]*size,[-1,-1]*size,[-1,1]*size,/t3d

  ;overplot the milky way disk and bulge
  if not keyword_set(nomw) then begin
    plots,diskx,disky,diskz,/t3d,color=blue
    polyfill,diskx,disky,diskz,/t3d,color=blue    ;,/line_fill
    plots,bx,by,bz,color=yellow,/t3d  ; bulge
  endif

  ; plotting the stars
  if nb lt 500 then begin
    for j=0,nb-1 do plots,xarr(j)+sphx,yarr(j)+sphy,zarr(j)+sphz,/t3d
  endif else begin
    plots,xarr,yarr,zarr,/t3d,ps=8
  endelse

  ;; overplot LMC orbit
  ;plots,lxarr,lyarr,lzarr,/t3d,color=green
  ;plots,lxarr+lx,lyarr+ly,lzarr+lz,/t3d,color=green
  ;plots,lmc(0:i,2),lmc(0:i,3),lmc(0:i,4),/t3d,color=green
  ;
  ; overplot SGR orbit
  ;plots,sxarr,syarr,szarr,/t3d,color=red
  ;plots,sxarr+sx,syarr+sy,szarr+sz,/t3d,color=red
  ;plots,sgr(0:i,2),sgr(0:i,3),sgr(0:i,4),/t3d,color=red

  ;stop
  wait,0.3

end

;; Z vs X  (upper-left)
;!p.multi=[4,2,2]
;plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=xr,yr=zr,xs=1,ys=1
;for i=0,nb-1 do begin
;  oplot,arr(*,i,2),arr(*,i,4),linestyle=0,co=co1,thick=thk
;  oplot,[arr(0,i,2)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
;end
;
;; R vs T (upper-right)
;!p.multi=[3,2,2]
;plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=tr,yr=rr,xs=1,ys=1
;for i=0,nb-1 do begin
;  oplot,arr(*,i,0),r(*,i),linestyle=0,co=co1,thick=thk
;end
;
;; Y vs X (lower-left)
;!p.multi=[2,2,2]
;plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=xr,yr=yr,xs=1,ys=1
;for i=0,nb-1 do begin
;  oplot,arr(*,i,2),arr(*,i,3),linestyle=0,co=co1,thick=thk
;  oplot,[arr(0,i,2)],[arr(0,i,3)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
;end
;
;; Z vs Y  (lower-right)
;!p.multi=[1,2,2]
;plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=yr,yr=zr,xs=1,ys=1
;for i=0,nb-1 do begin
;  oplot,arr(*,i,3),arr(*,i,4),linestyle=0,co=co1,thick=thk
;  oplot,[arr(0,i,3)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
;end
;
;; diagnostics
;if keyword_set(dfile) then begin
;  diag = importdiag(dfile)
;  etot = reform(diag(1,*))
;  ekin = reform(diag(2,*))
;  epot = reform(diag(3,*))
;endif

stop

end
