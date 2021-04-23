pro plotnbodymc,file,dfile,over=over,twoplot=twoplot,arr=arr,stp=stp,$
              xr=xr,yr=yr,zr=zr,rr=rr,tr=tr,iso=iso,movie=movie,$
              notrail=notrail,anim=anim,afile=afile,dir=dir,colarr=colarr,$
              first=first,ps=ps,dotfirst=dotfirst,last=last,R0=R0

;+
; This program plots the nbody orbits
;
; INPUTS:
;  file      File containing the galsat output data
;  dfile     File containing the galsat diagnostic data
;  /over
;  /twoplot  Only two panels
;  arr=arr   Input the galsat output array
;  /stp      Stop at end
;  xr=xr     X-axis range
;  yr=yr     Y-axis range
;  zr=zr     Z-axis range
;  rr=rr     Radius-axis range
;  tr=tr     Time-axis range
;  /iso      Make plots isotropic, aspect ratio=1
;  /movie    Flip through the timesteps
;  /notrail  If /movie set then don't show the trail for each body
;  /anim     Make an animation, save snaps to files
;  afile=afile Suffix for animation files
;  colarr=colarr  Array of colors for the bodies
;  /first   Only plot the first snap
;  /last      Only plot the last snap
;  /dotfirst  Make a dot for the first particle (the main body)
;
;  By David Nidever
;-

if keyword_set(anim) then movie=1

; importing nbody data
if not keyword_set(arr) then arr = importnbody(file)

if not keyword_set(R0) then R0=8.5

origarr = arr

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
for i=0.,nb-1 do begin
  for j=0.,ns-1 do r(j,i) = norm(reform(arr(j,i,2:4)))
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

; setting for postscript
if !d.name eq 'PS' or keyword_set(anim) then begin
  loadct,39
  black=0
  purple=30
  blue=60
  aqua=80
  green=155   ;135
  yellow=200
  orange=225
  white=0
  red=250
  lred=240
  ;backg=white
  backg=yellow
  coarr = [green,aqua,yellow,blue,purple,lred,orange]
endif

thk = 2      ; line thickness

;linestyle: 0-solid line, 1-dotted, 2-dashed, 3-dot dash, 4-dot dot
;          dot dash, 5-long dashes

psym8,0.7  ;1.5

co1 = green
co2 = red
co3 = orange
co4 = blue

; setting the ranges
dx=max(arr(*,*,2))-min(arr(*,*,2))
dy=max(arr(*,*,3))-min(arr(*,*,3))
dz=max(arr(*,*,4))-min(arr(*,*,4))
dr=max(r)-min(r)
dt=max(arr(*,*,0))-min(arr(*,*,0))
if not keyword_set(xr) then xr=[min(arr(*,*,2))-0.2*dx,max(arr(*,*,2))+0.2*dx]
if not keyword_set(yr) then yr=[min(arr(*,*,3))-0.2*dy,max(arr(*,*,3))+0.2*dy]
if not keyword_set(zr) then zr=[min(arr(*,*,4))-0.2*dz,max(arr(*,*,4))+0.2*dz]
if not keyword_set(rr) then rr=[min(r)-0.2*dr,max(r)+0.2*dr]
if not keyword_set(tr) then tr=[min(arr(*,*,0))-0.2*dt,max(arr(*,*,0))]

if keyword_set(iso) then begin
  lo = min([xr(0),yr(0),zr(0)])
  hi = max([xr(1),yr(1),zr(1)])
  xr = [lo,hi]
  yr = [lo,hi]
  zr = [lo,hi]
endif

co1 = white
thk = 0.5
if n_elements(ps) eq 0 then ps=8

dotsize = 1.5
dotco = 250

; MOVIE
if keyword_set(movie) then begin

  co2 = white
  if keyword_set(notrail) then co1=white else co1=red

  nsnap = n_elements(arr(*,0,0))
  psym8,0.7

  ; Looping through snaps
  for j=0.,nsnap-1 do begin

    ; Opening file or erasing
    if keyword_set(anim) then begin
      if not keyword_set(dir) then dir='plotnbody'
      if file_search(dir,/test_dir) eq '' then file_mkdir,dir
      if not keyword_set(afile) then afile = 'plotnbody'
      num = strtrim(long(j),2)
      if num lt 10 then num='0'+num
      if num lt 100 then num='0'+num
      if num lt 1000 then num='0'+num
      filename = dir+'/'+afile+'_'+num
      if j mod 25 eq 0 then ps_open,filename,/color
    endif; else erase

    ; Z vs X  (upper-left)
    !p.multi=[4,2,2]
    plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=xr,yr=zr,xs=1,ys=1
    for i=0.,nb-1 do begin
      ps1 = ps & co1 = 255 & symsize = 1.0
      if keyword_set(dotfirst) and i eq 0 then begin
        ps1=8
        co1=dotco
        symsize=dotsize
      endif

      ;Overplotting the last position (to erase the position's dot)
      oplot,[arr((j-1)>0,i,2)],[arr((j-1)>0,i,4)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize  


      ;Overplotting the Current position 
      oplot,[arr(j,i,2)],[arr(j,i,4)],linestyle=0,co=co1,thick=thk,ps=ps1,symsize=symsize

      ; Overplotting the trail
      if not keyword_set(notrail) then $
        oplot,[arr(0:j,i,2)],[arr(0:j,i,4)],linestyle=0,co=co2,thick=thk

      ; Overplotting the GC and Sun
      oplot,[0],[0],ps=1  ;GC
      oplot,[-R0],[0],ps=2  ;SUN
    end
    ; Overplot the body last so you can see it
    if keyword_set(dotfirst) then begin
      oplot,[arr(j,0,2)],[arr(j,0,4)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.7*dotsize
      oplot,[arr(j,1,2)],[arr(j,1,4)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.5*dotsize
      ;if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
    endif

    ; R vs T (upper-right)
    !p.multi=[3,2,2]
    plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=tr,yr=rr,xs=1,ys=1
    for i=0.,nb-1 do begin
      ps1 = ps & co1 = 255 & symsize = 1.0
      if keyword_set(dotfirst) and i eq 0 then begin
        ps1=8
        co1=dotco
        symsize=dotsize
      endif

      ; Overplotting the last position
      oplot,[arr((j-1)>0,i,0)],[r((j-1)>0,i)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize

      ; Overplotting the current position
      oplot,[arr(j,i,0)],[r(j,i)],linestyle=0,co=co1,thick=thk,co=co1,ps=ps1,symsize=symsize
      if not keyword_set(notrail) then oplot,[arr(0:j,i,0)],[r(0:j,i)],linestyle=0,co=co2,thick=thk
    end
    ; Overplot the body last so you can see it
    if keyword_set(dotfirst) then begin
      oplot,[arr(j,0,0)],[r(j,0)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.7*dotsize
      oplot,[arr(j,1,0)],[r(j,1)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.5*dotsize
      ;if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
    endif

    ; Y vs X (lower-left)
    !p.multi=[2,2,2]
    plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=xr,yr=yr,xs=1,ys=1
    for i=0.,nb-1 do begin
      ps1 = ps & co1 = 255 & symsize = 1.0
      if keyword_set(dotfirst) and i eq 0 then begin
        ps1=8
        co1=dotco
        symsize=dotsize
      endif

      ;Overplotting last position
      oplot,[arr((j-1)>0,i,2)],[arr((j-1)>0,i,3)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize

      ; Overplotting the current position
      oplot,[arr(j,i,2)],[arr(j,i,3)],linestyle=0,co=co1,thick=thk,ps=ps1,symsize=symsize

      if not keyword_set(notrail) then oplot,[arr(0:j,i,2)],[arr(0:j,i,3)],linestyle=0,co=co2,thick=thk
      oplot,[0],[0],ps=1  ; GC
      oplot,[-R0],[0],ps=2  ;SUN
    end
    ; Overplot the body last so you can see it
    if keyword_set(dotfirst) then begin
      oplot,[arr(j,0,2)],[arr(j,0,3)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.7*dotsize
      oplot,[arr(j,1,2)],[arr(j,1,3)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.5*dotsize
      ;if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
    endif

    ; Z vs Y  (lower-right)
    !p.multi=[1,2,2]
    plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=yr,yr=zr,xs=1,ys=1
    for i=0.,nb-1 do begin
      ps1 = ps & co1 = 255 & symsize = 1.0
      if keyword_set(dotfirst) and i eq 0 then begin
        ps1=8
        co1=dotco
        symsize=dotsize
      endif

      ; Overplotting the last position
      oplot,[arr((j-1)>0,i,3)],[arr((j-1)>0,i,4)],linestyle=0,co=0,thick=thk,ps=ps1,symsize=symsize

      ; Overplotting the current position
      oplot,[arr(j,i,3)],[arr(j,i,4)],linestyle=0,co=co1,thick=thk,ps=ps1,symsize=symsize

      if not keyword_set(notrail) then oplot,[arr(0:j,i,3)],[arr(0:j,i,4)],linestyle=0,co=co2,thick=thk
      oplot,[0],[0],ps=1  ; GC
      oplot,[0],[0],ps=2  ;SUN
    end
    ; Overplot the body last so you can see it
    if keyword_set(dotfirst) then begin
      oplot,[arr(j,0,3)],[arr(j,0,4)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.7*dotsize
      oplot,[arr(j,0,3)],[arr(j,0,4)],linestyle=0,co=co1,thick=thk,co=dotco,ps=8,symsize=0.5*dotsize
      ;if not keyword_set(notrail) then oplot,[arr(0:j,0,0)],[r(0:j,0)],linestyle=0,co=co2,thick=thk
    endif

    ; Closing file or waiting
    if keyword_set(anim) then begin
      if !d.name eq 'PS' then ps_close
    endif else begin
      ;if keyword_set(notrail) then wait,0.01 else wait,0.04
    endelse

    ;stop

  end ; for j

endif else begin

  ; NORMAL PLOTTING

  if keyword_set(last) then ind=ns-1 else ind = 0

  ; Z vs X  (upper-left)
  !p.multi=[4,2,2]
  plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=xr,yr=zr,xs=1,ys=1
  for i=0.,nb-1 do begin
    if keyword_set(colarr) then co1=colarr(i)

    ; All snaps
    if not keyword_set(first) and not keyword_set(last) then $
      oplot,arr(*,i,2),arr(*,i,4),linestyle=0,co=co1,thick=thk

    ; Only first or last snap
    if keyword_set(first) or keyword_set(last) then $
      oplot,[arr(ind,i,2)],[arr(ind,i,4)],co=co1,thick=thk,ps=ps

    if not keyword_set(last) then $
      oplot,[arr(0,i,2)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=ps    ;current position 
    oplot,[0],[0],ps=1  ;GC
    oplot,[-8,5],[0],ps=2  ;SUN

  end

  ; R vs T (upper-right)
  !p.multi=[3,2,2]
  plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=tr,yr=rr,xs=1,ys=1
  for i=0.,nb-1 do begin
    if keyword_set(colarr) then co1=colarr(i)

    ; All snaps
    if not keyword_set(first) and not keyword_set(last) then $
      oplot,arr(*,i,0),r(*,i),linestyle=0,co=co1,thick=thk

    ; Only first or last snap
    if keyword_set(first) or keyword_set(last) then $
      oplot,[arr(ind,i,0)],[r(ind,i)],co=co1,thick=thk,ps=ps

    if not keyword_set(last) then $
      oplot,[arr(0,i,0)],[r(0,i)],linestyle=0,co=co1,thick=thk,ps=ps        ; current position
  end

  ; Y vs X (lower-left)
  !p.multi=[2,2,2]
  plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=xr,yr=yr,xs=1,ys=1
  for i=0.,nb-1 do begin
    if keyword_set(colarr) then co1=colarr(i)

    ; All snaps
    if not keyword_set(first) and not keyword_set(last) then $
      oplot,arr(*,i,2),arr(*,i,3),linestyle=0,co=co1,thick=thk

    ; Only first or last snap
    if keyword_set(first) or keyword_set(last) then $
      oplot,[arr(ind,i,2)],[arr(ind,i,3)],co=co1,thick=thk,ps=ps

    if not keyword_set(last) then $
      oplot,[arr(0,i,2)],[arr(0,i,3)],linestyle=0,co=co1,thick=thk,ps=ps    ;current position
    oplot,[0],[0],ps=1  ; GC
    oplot,[-R0],[0],ps=2  ;SUN
  end

  ; Z vs Y  (lower-right)
  !p.multi=[1,2,2]
  plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=yr,yr=zr,xs=1,ys=1
  for i=0.,nb-1 do begin
    if keyword_set(colarr) then co1=colarr(i)

    ; All snaps
    if not keyword_set(first) and not keyword_set(last) then $
      oplot,arr(*,i,3),arr(*,i,4),linestyle=0,co=co1,thick=thk

    ; Only first or last snap
    if keyword_set(first) or keyword_set(last) then $
      oplot,[arr(ind,i,3)],[arr(ind,i,4)],co=co1,thick=thk,ps=ps

    if not keyword_set(last) then $
      oplot,[arr(0,i,3)],[arr(0,i,4)],linestyle=0,co=co1,thick=thk,ps=ps    ;current position

    oplot,[0],[0],ps=1  ; GC
    oplot,[0],[0],ps=2  ;SUN
  end

endelse

; diagnostics
if keyword_set(dfile) then begin
  diag = importdiag(dfile)
  etot = reform(diag(1,*))
  ekin = reform(diag(2,*))
  epot = reform(diag(3,*))
  steps = reform(diag(4,*))
endif

if keyword_set(stp) then stop

arr = origarr

end
