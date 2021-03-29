pro plotcoll,file,dfile,over=over,twoplot=twoplot,save=save,$
             tmin=tmini,mualpha=mualpha,mudelta=mudelta,vcirc=vcirc,$
             filename=filename,colrange=colrange,nodots=nodots,$
             nomin=nomin,title=title,dhalo=dhalo,lmcdist=lmcdist,arr=arr,$
             stp=stp

; This program plots the results of galsat with the LMC and Sgr.
;
; INPUTS:
;  arr      The galsat output array
;  file     the galsat output file
;  dfile    the galsat diagnostics file
;  /twoplot only plot two of the plots (under construction)
;  mualph   proper motion in ra * cos(delta) in mas/yr (for plotting)
;  mudelt   proper motion in dec. in mas/yr (for plotting)
;  vcirc    circular velocity of the disk (default=220 km/s)
;  dhalo    halo softening parameter, originally 13 kpc
;  lmcdist  distance to LMC
;  tmin     how far back to run the simulation (default=4 Gyrs)
;           This number should be positive (but can handle negative).
;  /save    save the plot to postscript
;  filename name of postscript file (without the .ps)
;  colrange time range for minmum separation
;  /over    overplot stiched Sgr, solitary Sgr, and old LMC orbits.
;  /nodots  don't overplot the equally spaced points
;  /nomin   don't overplot the minimum separation points
;  title    overall title of the plot
;
;  Created by D.Nidever Jan. 2005


;over = 1
if n_elements(twoplot) eq 0 then twoplot = 0

tmin = -4.
if keyword_set(tmini) then tmin = tmini
if tmin gt 0. then tmin=-tmin

; importing nbody data
if not keyword_set(arr) then arr = importnbody(file)
; if there there are softening parameters then get rid of them.
if (n_elements(arr(0,0,*)) eq 9) then begin
  dum = arr
  arr = arr(*,*,0:7)*0.
  arr(*,*,0:1) = dum(*,*,0:1)
  arr(*,*,2:7) = dum(*,*,3:8)
  dum = 0.
endif
ns = n_elements(arr(*,0,0))
; R for LMC
r = dblarr(ns)
for i=0,ns-1 do r(i) = norm(reform(arr(i,0,2:4)))
; R for SGR
r2 = dblarr(ns)
for i=0,ns-1 do r2(i) = norm(reform(arr(i,1,2:4)))

;; TIME and 0.1 Gyr steps
;t = arr(*,0,0)
;tind = lonarr(6)
;for i=0,5 do begin
;  ti = -i*0.1-1.6
;  tind(i) = first_el(minloc(abs(t-ti),/first))
;end
;gdind = where(t gt tmin)

; TIME in 1 Gyr steps
t = arr(*,0,0)
nt1 = long(abs(tmin))-1
t1ind = lonarr(nt1)
for i=0,nt1-1 do begin
  ti = -i-1
  t1ind(i) = first_el(minloc(abs(t-ti),/first))
end
; TIME and in 0.2 Gyr steps
nt2 = long(abs(tmin)/0.2)-1
t2ind = lonarr(nt2)
for i=0,nt2-1 do begin
  ti = -(i+1)*0.2
  t2ind(i) = first_el(minloc(abs(t-ti),/first))
end
gdind = where(t gt tmin)

; SEPARATION between LMC and SGR
sep = dblarr(ns)
for i=0,ns-1 do sep(i) = norm(reform(arr(i,0,2:4)-arr(i,1,2:4)))
;print,'Minimum Separation Distance = ',stringize(min(sep),ndec=2),' kpc'

; time of closest approach
if not keyword_set(colrange) then colrange=[-2.3,-1.0]
if min(colrange) gt 0. then colrange=-colrange
tbrack = where(t gt colrange(0) and t lt colrange(1))
;tbrack = where(t gt -2.3 and t lt -1.0)
colind = [minloc(sep(tbrack),/first)]+min(tbrack)
print,'Minimum Separation Distance = ',stringize(min(sep(colind)),ndec=2),' kpc'

; printing stats at time of closest approach
lmcpos = reform(arr(colind,0,2:4))
sgrpos = reform(arr(colind,1,2:4))
print,''
print,'Positions at time of closest approach'
print,'LMC: ',stringize(lmcpos(0),ndec=2),' ',stringize(lmcpos(1),ndec=2),' ',stringize(lmcpos(2),ndec=2)
print,'SGR: ',stringize(sgrpos(0),ndec=2),' ',stringize(sgrpos(1),ndec=2),' ',stringize(sgrpos(2),ndec=2)

; importing normal SGR orbit
sgrarr = importnbody('sgr.out')
; if there there are softening parameters then get rid of them.
if (n_elements(sgrarr(0,0,*)) eq 9) then begin
  dum = sgrarr
  sgrarr = sgrarr(*,*,0:7)*0.
  sgrarr(*,*,0:1) = dum(*,*,0:1)
  sgrarr(*,*,2:7) = dum(*,*,3:8)
  dum = 0.
endif
ns2 = n_elements(sgrarr(*,0,0))
rs = dblarr(ns2)
for i=0,ns2-1 do rs(i) = norm(reform(sgrarr(i,0,2:4)))

; importing maketrail data
;marr = importdat('orbSGR.q9.dat')      ;t,x,y,z,vx,vy,vz,a,b,c
;;marr = importdat('orbLMC_121304.dat')      ;t,x,y,z,vx,vy,vz,a,b,c
;nt = n_elements(marr(*,0))
;mr = dblarr(nt)
;for i=0,nt-1 do mr(i) = norm(reform(marr(i,1:3)))
;gd = where(marr(*,0) le 0 and marr(*,0) ge -2.)

; overplotting the LMC and SGR orbits
if keyword_set(over) then begin

  ; Loading the data
  sgrarr1 = importdat2('orb1.0.dat')
  sgrarr2 = importdat2('orb321.dat')
  lmcarr = importdat2('orbLMC_121304.dat')

  sgrarr1(0,*) = sgrarr1(0,*)-2.0  ;outer Sgr kludge

  ; selecting the right time ranges
  sgrind1 = where(sgrarr1(0,*) gt -6.1 and sgrarr1(0,*) lt -2.0)
  sgrind2 = where(sgrarr2(0,*) ge -2.0 and sgrarr2(0,*) le 0.1)
  lmcind = where(lmcarr(0,*) gt -6.1 and lmcarr(0,*) le 0.1)

  ; pasting the two Sgr orbits together
  sgr = [[sgrarr1(*,sgrind1)],[sgrarr2(*,sgrind2)]]
  lmc = lmcarr(*,lmcind)

  ; LMC stuff
  nlmc = n_elements(lmc(0,*))
  rlmc = dblarr(nlmc)
  for i=0,nlmc-1 do rlmc(i) = norm(reform(lmc(1:3,i)))
  gdlmc = where(lmc(0,*) le 0. and lmc(0,*) ge tmin,ngdlmc)

  ; SGR stuff
  nsgr = n_elements(sgr(0,*))
  rsgr = dblarr(nsgr)
  for i=0,nsgr-1 do rsgr(i) = norm(reform(sgr(1:3,i)))
  gdsgr = where(sgr(0,*) le 0. and sgr(0,*) ge tmin,ngdsgr)

endif

;save,lmc,nlmc,rlmc,gdlmc,sgr,nsgr,rsgr,gdsgr,sgrarr,ns2,rs,file='coll_orbits.dat'
;restore,'coll_orbits.dat'

if not keyword_set(save) then erase  ; erase plot window

sym8a = 1.2
sym8b = 1.2

; setting for postscript
if keyword_set(save) then begin
  loadct,13
  black=0
  purple=30
  blue=60
  aqua=80
  green=155   ;135
  yellow=205  ;200
  orange=225
  white=0
  red=250    ;300
  backg=white
  foreg=black

  co1 = green
  co2 = red
  co3 = orange
  co4 = blue

  sym8a = sym8a*0.8
  sym8b = sym8b*0.7
  thk = 3.0
  thk2 = 5.0

endif else begin

  ;; colors for my screen
  ;red = 250
  ;lred = 210
  ;green = 190000
  ;orange = 310000
  ;yellow = 450000
  ;blue = -25000
  ;lblue = -15000
  ;purple = -20000
  ;white = -1

  ; decomposed screen colors
  device,decomposed=0
  loadct,39
  black=1
  purple=30
  dblue = 70
  blue=90  ;60
  lblue = 85
  aqua=80
  green=155   ;135
  yellow=195  ;190
  orange=210
  white=!p.color
  red=250    ;300
  lred=217
  backg=black
  foreg=white

endelse

; Setting colors
co1 = green  ;green   ;lmc
co2 = red            ; sgr
co3 = blue   ;orange   ; previous lmc orbit
co4 = orange   ;blue    ; sgr to fit
co5 = yellow  ;purple  ;lred  ; solitary sgr orbit

; Setting line thicknesses
thk = 2      ; line thickness
thk2 = 2
charthk = 1
xthk = 1
ythk = 1
if keyword_set(save) then begin
  charthk = 3
  xthk = 3
  ythk = 3
  thk = 4
  thk2 = 4
endif

; character size
charsize = 1.3
if keyword_set(save) then charsize=charsize*0.8
;if keyword_set(save) then charsize=charsize*0.9

; symbol size for equally spaced points
sym1 = 1.3
sym2 = 0.9
;if keyword_set(save) then sym1 = sym1*0.8
;if keyword_set(save) then sym2 = sym2*0.8


;linestyle: 0-solid line, 1-dotted, 2-dashed, 3-dot dash, 4-dot dot
;          dot dash, 5-long dashes

if not keyword_set(filename) then filename = 'plotcoll'
if keyword_set(save) then ps_open,filename,/color

; Setting the scale limits
if keyword_set(over) then begin
  xx = [arr(gdind,0,2),arr(gdind,1,2),reform(lmc(1,gdlmc)),reform(sgr(1,gdsgr)),sgrarr(*,0,2)]
  yy = [arr(gdind,0,3),arr(gdind,1,3),reform(lmc(2,gdlmc)),reform(sgr(2,gdsgr)),sgrarr(*,0,3)]
  zz = [arr(gdind,0,4),arr(gdind,1,4),reform(lmc(3,gdlmc)),reform(sgr(3,gdsgr)),sgrarr(*,0,4)]
  rr = [r(gdind),r2(gdind),rlmc(gdlmc),rsgr(gdsgr),rs(*),sep]
  xrange = [-5.+min(xx),max(xx)+5.]
  yrange = [-5.+min(yy),max(yy)+5.]
  zrange = [-5.+min(zz),max(zz)+5.]
  rrange = [-5.+min(rr),max(rr)+5.]
endif else begin
  xx = [arr(gdind,0,2),arr(gdind,1,2)]
  yy = [arr(gdind,0,3),arr(gdind,1,3)]
  zz = [arr(gdind,0,4),arr(gdind,1,4)]
  rr = [r(gdind),r2(gdind),sep]
  xrange = [-5.+min(xx),max(xx)+5.]
  yrange = [-5.+min(yy),max(yy)+5.]
  zrange = [-5.+min(zz),max(zz)+5.]
  rrange = [-5.+min(rr),max(rr)+5.]
endelse


;stop

; Z vs X  (upper-left)
if not keyword_set(twoplot) then begin
!p.multi=[4,2,2]
; alter the !x.window and !y.window to set the position
;!y.window(0) = !y.window(0)-0.05
;!y.window(1) = !y.window(1)-0.025
psym8,sym8a,/fill,/square
; xr=[-70,70],yr=[-110,110]
plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=xrange,yr=zrange,xs=1,ys=1,$
     /normal,position=[0.08,0.59-0.025,0.475,0.985-0.05],charth=charthk,xth=xthk,yth=ythk
oplot,[0],[0],ps=1       ;origin
oplot,arr(gdind,0,2),arr(gdind,0,4),linestyle=0,co=co1,thick=thk
oplot,[arr(0,0,2)],[arr(0,0,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,2),arr(gdind,1,4),linestyle=0,co=co2,thick=thk2
oplot,[arr(0,1,2)],[arr(0,1,4)],linestyle=0,co=co2,thick=thk2,ps=8    ;current position
;oplot,marr(gd,1),marr(gd,3),linestyle=2


;; equally spaced points in time
;psym8,sym8b
;oplot,arr(tind,0,2),arr(tind,0,4),co=co1,ps=8
;oplot,arr(tind,1,2),arr(tind,1,4),co=co2,ps=8

; equally spaced points in time
if not keyword_set(nodots) then begin
  psym8,sym8b
  ;sym1 = 0.8
  ;sym2 = 0.5
  oplot,arr(t1ind,0,2),arr(t1ind,0,4),co=co1,ps=8,symsize=sym1
  oplot,arr(t2ind,0,2),arr(t2ind,0,4),co=co1,ps=8,symsize=sym2
  oplot,arr(t1ind,1,2),arr(t1ind,1,4),co=co2,ps=8,symsize=sym1
  oplot,arr(t2ind,1,2),arr(t2ind,1,4),co=co2,ps=8,symsize=sym2
endif

if keyword_set(over) then begin
  ;oplot,lmc(1,gdlmc),lmc(3,gdlmc),linestyle=5,co=co3,thick=thk
  oplot,sgr(1,gdsgr),sgr(3,gdsgr),linestyle=5,co=co4,thick=thk
  ;oplot,sgrarr(*,0,2),sgrarr(*,0,4),co=co5,thick=0.2
endif

; point of closest approach
if not keyword_set(nomin) then begin
  starco1 = foreg
  starco2 = yellow
  starsym1 = 1.7
  starsym2 = 0.8
  oplot,arr(colind,0,2),arr(colind,0,4),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,0,2),arr(colind,0,4),co=starco2,ps=2,symsize=starsym2
  oplot,arr(colind,1,2),arr(colind,1,4),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,1,2),arr(colind,1,4),co=starco2,ps=2,symsize=starsym2
endif

endif ; not twoplot


; R vs T (upper-right)
!p.multi=[3,2,2]
;!y.window(0) = !y.window(0)-0.05
;!y.window(1) = !y.window(1)-0.025
if keyword_set(twoplot) then !p.multi=[2,2,0]
psym8,sym8a,/fill,/square
; xr=[tmin,0],yr=[0,110]
plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=[tmin,0],yr=rrange,xs=1,ys=1,$
     /normal,position=[0.58,0.59-0.025,0.975,0.985-0.05],charth=charthk,xth=xthk,yth=ythk
oplot,arr(gdind,0,0),r(gdind),linestyle=0,co=co1,thick=thk
;oplot,[arr(0,0,0)],[r(0)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,0),r2(gdind),linestyle=0,co=co2,thick=thk2
;oplot,[arr(0,1,0)],[r2(0)],linestyle=0,co=co2,thick=thk,ps=8    ;current position
;oplot,marr(gd,0),mr(gd),linestyle=2


;psym8,sym8b
;oplot,arr(tind,0,0),r(tind),co=co1,ps=8
;oplot,arr(tind,1,0),r2(tind),co=co2,ps=8

; equally spaced points in time
if not keyword_set(nodots) then begin
  psym8,sym8b
  ;sym1 = 0.8
  ;sym2 = 0.5
  oplot,arr(t1ind,0,0),r(t1ind),co=co1,ps=8,symsize=sym1
  oplot,arr(t2ind,0,0),r(t2ind),co=co1,ps=8,symsize=sym2
  oplot,arr(t1ind,1,0),r2(t1ind),co=co2,ps=8,symsize=sym1
  oplot,arr(t2ind,1,0),r2(t2ind),co=co2,ps=8,symsize=sym2
endif

oplot,t,sep,co=white

if keyword_set(over) then begin
  ;oplot,lmc(0,gdlmc),rlmc(gdlmc),linestyle=5,co=co3,thick=thk
  oplot,sgr(0,gdsgr),rsgr(gdsgr),linestyle=5,co=co4,thick=thk
  ;oplot,sgrarr(*,0,0),rs(*),co=co5,thick=0.2
endif

; point of closest approach
if not keyword_set(nomin) then begin
  oplot,arr(colind,0,0),r(colind),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,0,0),r(colind),co=starco2,ps=2,symsize=starsym2
  oplot,arr(colind,1,0),r2(colind),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,1,0),r2(colind),co=starco2,ps=2,symsize=starsym2
endif

; Y vs X (lower-left)
if not keyword_set(twoplot) then begin
!p.multi=[2,2,2]
;!y.window(1) = !y.window(1)-0.025
psym8,sym8a,/fill,/square
;xr=[-70,70],yr=[-110,110]
plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=xrange,yr=yrange,xs=1,ys=1,$
     /normal,position=[0.08,0.07,0.475,0.49-0.025],charth=charthk,xth=xthk,yth=ythk

oplot,[0],[0],ps=1       ;origin
oplot,arr(gdind,0,2),arr(gdind,0,3),linestyle=0,co=co1,thick=thk
oplot,[arr(0,0,2)],[arr(0,0,3)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,2),arr(gdind,1,3),linestyle=0,co=co2,thick=thk2
oplot,[arr(0,1,2)],[arr(0,1,3)],linestyle=0,co=co2,thick=thk,ps=8    ;current position
;oplot,marr(gd,1),marr(gd,2),linestyle=2

;psym8,sym8b
;oplot,arr(tind,0,2),arr(tind,0,3),co=co1,ps=8
;oplot,arr(tind,1,2),arr(tind,1,3),co=co2,ps=8

; equally spaced points in time
if not keyword_set(nodots) then begin
  psym8,sym8b
  ;sym1 = 0.8
  ;sym2 = 0.5
  oplot,arr(t1ind,0,2),arr(t1ind,0,3),co=co1,ps=8,symsize=sym1
  oplot,arr(t2ind,0,2),arr(t2ind,0,3),co=co1,ps=8,symsize=sym2
  oplot,arr(t1ind,1,2),arr(t1ind,1,3),co=co2,ps=8,symsize=sym1
  oplot,arr(t2ind,1,2),arr(t2ind,1,3),co=co2,ps=8,symsize=sym2
endif

if keyword_set(over) then begin
  ;oplot,lmc(1,gdlmc),lmc(2,gdlmc),linestyle=5,co=co3,thick=thk
  oplot,sgr(1,gdsgr),sgr(2,gdsgr),linestyle=5,co=co4,thick=thk
  ;oplot,sgrarr(*,0,2),sgrarr(*,0,3),co=co5,thick=0.2
endif

; point of closest approach
if not keyword_set(nomin) then begin
  oplot,arr(colind,0,2),arr(colind,0,3),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,0,2),arr(colind,0,3),co=starco2,ps=2,symsize=starsym2
  oplot,arr(colind,1,2),arr(colind,1,3),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,1,2),arr(colind,1,3),co=starco2,ps=2,symsize=starsym2
endif

end

; Z vs Y  (lower-right)
!p.multi=[1,2,2]
if keyword_set(twoplot) then !p.multi=[1,2,0]
;!y.window(1) = !y.window(1)-0.025
psym8,sym8a,/fill,/square
;xr=[-110,110],yr=[-110,110]
;plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=[-70,70],yr=[-70,70],xs=1,ys=1
plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=yrange,yr=zrange,xs=1,ys=1,$
     /normal,position=[0.58,0.07,0.975,0.49-0.025],charth=charthk,xth=xthk,yth=ythk
oplot,[0],[0],ps=1       ;origin
oplot,arr(gdind,0,3),arr(gdind,0,4),linestyle=0,co=co1,thick=thk
oplot,[arr(0,0,3)],[arr(0,0,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,3),arr(gdind,1,4),linestyle=0,co=co2,thick=thk2
oplot,[arr(0,1,3)],[arr(0,1,4)],linestyle=0,co=co2,thick=thk,ps=8    ;current position
;oplot,marr(gd,2),marr(gd,3),linestyle=2

;psym8,sym8b
;oplot,arr(tind,0,3),arr(tind,0,4),co=co1,ps=8
;oplot,arr(tind,1,3),arr(tind,1,4),co=co2,ps=8

; equally spaced points in time
if not keyword_set(nodots) then begin
  psym8,sym8b
  ;sym1 = 0.8
  ;sym2 = 0.5
  oplot,arr(t1ind,0,3),arr(t1ind,0,4),co=co1,ps=8,symsize=sym1
  oplot,arr(t2ind,0,3),arr(t2ind,0,4),co=co1,ps=8,symsize=sym2
  oplot,arr(t1ind,1,3),arr(t1ind,1,4),co=co2,ps=8,symsize=sym1
  oplot,arr(t2ind,1,3),arr(t2ind,1,4),co=co2,ps=8,symsize=sym2
endif

if keyword_set(over) then begin
  ;oplot,lmc(2,gdlmc),lmc(3,gdlmc),linestyle=5,co=co3,thick=thk
  oplot,sgr(2,gdsgr),sgr(3,gdsgr),linestyle=5,co=co4,thick=thk
  ;oplot,sgrarr(*,0,3),sgrarr(*,0,4),co=co5,thick=0.2
endif

; point of closest approach
if not keyword_set(nomin) then begin
  oplot,arr(colind,0,3),arr(colind,0,4),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,0,3),arr(colind,0,4),co=starco2,ps=2,symsize=starsym2
  oplot,arr(colind,1,3),arr(colind,1,4),co=starco1,ps=2,symsize=starsym1
  oplot,arr(colind,1,3),arr(colind,1,4),co=starco2,ps=2,symsize=starsym2
endif

; overplotting labels
align=0.5
normal = 1

;xyouts,0.20,0.485,'Model LMC',charsize=charsize,normal=normal,alignment=align,co=co1,charth=charthk
;xyouts,0.40,0.485,'Model SGR',charsize=charsize,normal=normal,alignment=align,co=co2,charth=charthk
;;xyouts,0.60,0.485,'Stitched SGR',charsize=charsize,normal=normal,alignment=align,co=co4,charth=charthk
;xyouts,0.80,0.485,'SGR-LMC Sep.',charsize=charsize,normal=normal,alignment=align,co=white,charth=charthk
;title=' '

if keyword_set(over) then begin
  xyouts,0.30,0.485,'Model LMC',charsize=charsize,normal=normal,alignment=align,co=co1,charth=charthk
  xyouts,0.45,0.485,'Model SGR',charsize=charsize,normal=normal,alignment=align,co=co2,charth=charthk
  xyouts,0.60,0.485,'Stitched SGR',charsize=charsize,normal=normal,alignment=align,co=co5,charth=charthk
  xyouts,0.78,0.485,'SGR-LMC Sep.',charsize=charsize,normal=normal,alignment=align,co=white,charth=charthk

  ;xyouts,0.20,0.485,'Model LMC',charsize=charsize,normal=normal,alignment=align,co=co1,charth=charthk
  ;xyouts,0.32,0.485,'Old LMC',charsize=charsize,normal=normal,alignment=align,co=co3,charth=charthk
  ;xyouts,0.44,0.485,'Model SGR',charsize=charsize,normal=normal,alignment=align,co=co2,charth=charthk
  ;if abs(vcirc-220.) gt 4. then begin
  ;  xyouts,0.56,0.485,'Stitched SGR!d220!n',charsize=charsize,normal=normal,alignment=align,co=co4,charth=charthk
  ;  xyouts,0.69,0.485,'Solitary SGR!d220!n',charsize=charsize,normal=normal,alignment=align,co=co5,charth=charthk
  ;endif else begin
  ;  xyouts,0.56,0.485,'Stitched SGR',charsize=charsize,normal=normal,alignment=align,co=co4,charth=charthk
  ;  xyouts,0.69,0.485,'Solitary SGR',charsize=charsize,normal=normal,alignment=align,co=co5,charth=charthk
  ;endelse
  ;xyouts,0.83,0.485,'SGR-LMC Sep.',charsize=charsize,normal=normal,alignment=align,co=white,charth=charthk
endif else begin
  xyouts,0.30,0.485,'Model LMC',charsize=charsize,normal=normal,alignment=align,co=co1,charth=charthk
  xyouts,0.50,0.485,'Model SGR',charsize=charsize,normal=normal,alignment=align,co=co2,charth=charthk
  xyouts,0.70,0.485,'SGR-LMC Sep.',charsize=charsize,normal=normal,alignment=align,co=white,charth=charthk
endelse

; overplotting the input parameters, mualpha, mudelta, vsini
if not keyword_set(title) then begin
  strvcirc = '!3V!dcirc!n = '+stringize(vcirc,ndec=2)+' km s!u-1!n'
  ;tit = 'Minimum Separation Distance ['+strimin+', '+strimax+']'
  strmualpha = '!4l!da!n!3 cos(!4d!3) = '+stringize(mualpha,ndec=3)+' mas yr!u-1!n'
  strmudelta = '!4l!dd!n!3 = '+stringize(mudelta,ndec=3)+' mas yr!u-1!n'
  strminsep = 'd!dmin!n = '+stringize(min(sep(colind)),ndec=1)+' kpc'
  strdhalo = 'Dhalo = '+stringize(dhalo,ndec=2)+' kpc'
  xyouts,0.2,0.96,strmualpha,charsize=charsize,normal=normal,alignment=align,charth=charthk
  xyouts,0.44,0.96,strmudelta,charsize=charsize,normal=normal,alignment=align,charth=charthk
  xyouts,0.66,0.96,strvcirc,charsize=charsize,normal=normal,alignment=align,charth=charthk
  xyouts,0.85,0.96,strminsep,charsize=charsize,normal=normal,alignment=align,charth=charthk
  ;xyouts,0.85,0.96,strdhalo,charsize=charsize,normal=normal,alignment=align,charth=charthk
endif else begin
  xyouts,0.5,0.96,title,charsize=charsize,normal=normal,alignment=0.5,charth=charthk
endelse

if keyword_set(save) then ps_close

; diagnostics
if keyword_set(dfile) then begin
  diag = importdiag(dfile)
  etot = reform(diag(1,*))
  ekin = reform(diag(2,*))
  epot = reform(diag(3,*))
endif

if keyword_set(stp) then stop

end
