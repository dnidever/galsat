pro plotlmccar,file,dfile,over=over,twoplot=twoplot,save=save

; This program plots the nbody orbits for LMC, SMC, and SGR

over = 1
twoplot = 0

tmin = -7.

; importing nbody data
arr = importnbody(file)
ns = n_elements(arr(*,0,0))
nb = n_elements(arr(0,*,0))
; R for all
r = dblarr(ns,nb)
for i=0,nb-1 do begin
  for j=0,ns-1 do r(j,i) = norm(reform(arr(j,i,2:4)))
end

;; R for LMC
;r = dblarr(ns)
;for i=0,ns-1 do r(i) = norm(reform(arr(i,0,2:4)))
;; R for SGR
;r2 = dblarr(ns)
;for i=0,ns-1 do r2(i) = norm(reform(arr(i,1,2:4)))
; TIME and 1 Gyr steps
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
print,'Minimum Separation Distance LMC-Car = ',stringize(min(sep),ndec=2),' kpc'
; Separation between SMC and SGR
;sep2 = dblarr(ns)
;for i=0,ns-1 do sep2(i) = norm(reform(arr(i,1,2:4)-arr(i,2,2:4)))
;print,'Minimum Separation Distance SMC-Sgr = ',stringize(min(sep2),ndec=2),' kpc'
; Separation between LMC and SMC
;sep3 = dblarr(ns)
;for i=0,ns-1 do sep3(i) = norm(reform(arr(i,0,2:4)-arr(i,2,2:4)))
;print,'Minimum Separation Distance LMC-SMC = ',stringize(min(sep3),ndec=2),' kpc'

; time of closest approach
;tbrack = where(t gt -2.3 and t lt -1.0)
;colind = [minloc(sep(tbrack),/first)]+min(tbrack)
colind = [minloc(sep,/first)]

; printing stats at time of closest approach
lmcpos = reform(arr(colind,0,2:4))
sgrpos = reform(arr(colind,1,2:4))
print,''
print,'Positions at time of closest approach'
print,'LMC: ',stringize(lmcpos(0),ndec=2),' ',stringize(lmcpos(1),ndec=2),' ',stringize(lmcpos(2),ndec=2)
print,'Car: ',stringize(sgrpos(0),ndec=2),' ',stringize(sgrpos(1),ndec=2),' ',stringize(sgrpos(2),ndec=2)

; importing normal SGR orbit
;sgrarr = importnbody('sgr.out')
;ns2 = n_elements(sgrarr(*,0,0))
;rs = dblarr(ns2)
;for i=0,ns2-1 do rs(i) = norm(reform(sgrarr(i,0,2:4)))

; importing maketrail data
;marr = importdat('orbSGR.q9.dat')      ;t,x,y,z,vx,vy,vz,a,b,c
;;marr = importdat('orbLMC_121304.dat')      ;t,x,y,z,vx,vy,vz,a,b,c
;nt = n_elements(marr(*,0))
;mr = dblarr(nt)
;for i=0,nt-1 do mr(i) = norm(reform(marr(i,1:3)))
;gd = where(marr(*,0) le 0 and marr(*,0) ge -2.)

; overplotting the LMC and SGR orbits
;if keyword_set(over) then begin

  cararr = importnbody('car.out')
  gdcar = where(cararr(*,0,0) le 0. and cararr(*,0,0) ge tmin,ngdcar)

  ncar = n_elements(cararr(*,0,0))
  rcar = dblarr(ncar)
  for i=0,ncar-1 do rcar(i) = norm(reform(cararr(i,0,2:4)))

;  ; Loading the data
;  sgrarr1 = importdat2('orb1.0.dat')
;  sgrarr2 = importdat2('orb321.dat')
;  lmcarr = importdat2('orbLMC_121304.dat')
;
;  sgrarr1(0,*) = sgrarr1(0,*)-2.0  ;outer Sgr kludge
;
;  ; selecting the right time ranges
;  sgrind1 = where(sgrarr1(0,*) gt -6.1 and sgrarr1(0,*) lt -2.0)
;  sgrind2 = where(sgrarr2(0,*) ge -2.0 and sgrarr2(0,*) le 0.1)
;  lmcind = where(lmcarr(0,*) gt -6.1 and lmcarr(0,*) le 0.1)
;
;  ; pasting the two Sgr orbits together
;  sgr = [[sgrarr1(*,sgrind1)],[sgrarr2(*,sgrind2)]]
;  lmc = lmcarr(*,lmcind)
;
;  ; LMC stuff
;  nlmc = n_elements(lmc(0,*))
;  rlmc = dblarr(nlmc)
;  for i=0,nlmc-1 do rlmc(i) = norm(reform(lmc(1:3,i)))
;  gdlmc = where(lmc(0,*) le 0. and lmc(0,*) ge tmin,ngdlmc)
;
;  ; SGR stuff
;  nsgr = n_elements(sgr(0,*))
;  rsgr = dblarr(nsgr)
;  for i=0,nsgr-1 do rsgr(i) = norm(reform(sgr(1:3,i)))
;  gdsgr = where(sgr(0,*) le 0. and sgr(0,*) ge tmin,ngdsgr)
;
;endif

;save,lmc,nlmc,rlmc,gdlmc,sgr,nsgr,rsgr,gdsgr,sgrarr,ns2,rs,file='coll_orbits.dat'
;restore,'coll_orbits.dat'

if not keyword_set(save) then erase  ; erase plot window

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


thk = 1.2      ; line thickness
thk2 = 1.2

;linestyle: 0-solid line, 1-dotted, 2-dashed, 3-dot dash, 4-dot dot
;          dot dash, 5-long dashes

if keyword_set(save) then ps_open,'plotlmccar',/color

sym8a = 1.2
sym8b = 1.2

co1 = green
co2 = red
co3 = orange
co4 = blue

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
  red=250    ;300
  backg=white

  co1 = green
  co2 = red
  co3 = orange
  co4 = blue

  sym8a = sym8a*0.8
  sym8b = sym8b*0.7
  thk = 3.0
  thk2 = 5.0
end


; Z vs X  (upper-left)
if not keyword_set(twoplot) then begin
!p.multi=[4,2,2]
psym8,sym8a,/fill,/square
plot,dist(5),/nodata,xtit='XGC',ytit='ZGC',xr=[-70,70],yr=[-130,130],xs=1,ys=1
oplot,[0],[0],ps=1       ;origin
oplot,arr(gdind,0,2),arr(gdind,0,4),linestyle=0,co=co1,thick=thk
oplot,[arr(0,0,2)],[arr(0,0,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,2),arr(gdind,1,4),linestyle=0,co=co2,thick=thk2
oplot,[arr(0,1,2)],[arr(0,1,4)],linestyle=0,co=co2,thick=thk2,ps=8    ;current position
;oplot,arr(gdind,2,2),arr(gdind,2,4),linestyle=0,co=co3,thick=thk
;oplot,[arr(0,2,2)],[arr(0,2,4)],linestyle=0,co=co3,thick=thk2,ps=8    ;current position
;oplot,marr(gd,1),marr(gd,3),linestyle=2
;oplot,sgrarr(*,0,2),sgrarr(*,0,4),co=white,thick=0.2

; equally spaced points in time
psym8,sym8b
sym1 = 0.8
sym2 = 0.5
oplot,arr(t1ind,0,2),arr(t1ind,0,4),co=co1,ps=8,symsize=sym1
oplot,arr(t2ind,0,2),arr(t2ind,0,4),co=co1,ps=8,symsize=sym2
oplot,arr(t1ind,1,2),arr(t1ind,1,4),co=co2,ps=8,symsize=sym1
oplot,arr(t2ind,1,2),arr(t2ind,1,4),co=co2,ps=8,symsize=sym2
;oplot,arr(t1ind,2,2),arr(t1ind,2,4),co=co3,ps=8,symsize=sym1
;oplot,arr(t2ind,2,2),arr(t2ind,2,4),co=co3,ps=8,symsize=sym2

if keyword_set(over) then begin
;  oplot,lmc(1,gdlmc),lmc(3,gdlmc),linestyle=5,co=co3,thick=thk
;  oplot,sgr(1,gdsgr),sgr(3,gdsgr),linestyle=5,co=co4,thick=thk
  oplot,cararr(gdcar,0,2),cararr(gdcar,0,4),linestyle=5,co=white
endif

; point of closest approach
oplot,arr(colind,0,2),arr(colind,0,4),co=black,ps=2,symsize=1.7
oplot,arr(colind,0,2),arr(colind,0,4),co=yellow,ps=2,symsize=0.8
oplot,arr(colind,1,2),arr(colind,1,4),co=black,ps=2,symsize=1.7
oplot,arr(colind,1,2),arr(colind,1,4),co=yellow,ps=2,symsize=0.8

endif ; not twoplot


; R vs T (upper-right)
!p.multi=[3,2,2]
if keyword_set(twoplot) then !p.multi=[2,2,0]
psym8,sym8a,/fill,/square
plot,dist(5),/nodata,xtit='t',ytit='RGC',xr=[tmin,0],yr=[0,150],xs=1,ys=1
oplot,arr(gdind,0,0),r(gdind,0),linestyle=0,co=co1,thick=thk
;oplot,[arr(0,0,0)],[r(0)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,0),r(gdind,1),linestyle=0,co=co2,thick=thk2
;oplot,arr(gdind,2,0),r(gdind,2),linestyle=0,co=co3,thick=thk
;oplot,[arr(0,1,0)],[r2(0)],linestyle=0,co=co2,thick=thk,ps=8    ;current position
;oplot,marr(gd,0),mr(gd),linestyle=2

;psym8,sym8b
;oplot,arr(tind,0,0),r(tind),co=co1,ps=8
;oplot,arr(tind,1,0),r2(tind),co=co2,ps=8

oplot,t,sep,co=white
;oplot,t,sep2,co=white,linestyle=2
;oplot,t,sep3,co=white,linestyle=3

; equally spaced points in time
psym8,sym8b
sym1 = 0.8
sym2 = 0.5
oplot,arr(t1ind,0,0),r(t1ind,0),co=co1,ps=8,symsize=sym1
oplot,arr(t2ind,0,0),r(t2ind,0),co=co1,ps=8,symsize=sym2
oplot,arr(t1ind,1,0),r(t1ind,1),co=co2,ps=8,symsize=sym1
oplot,arr(t2ind,1,0),r(t2ind,1),co=co2,ps=8,symsize=sym2
;oplot,arr(t1ind,2,0),r(t1ind,2),co=co3,ps=8,symsize=sym1
;oplot,arr(t2ind,2,0),r(t2ind,2),co=co3,ps=8,symsize=sym2

if keyword_set(over) then begin
;  oplot,lmc(0,gdlmc),rlmc(gdlmc),linestyle=5,co=co3,thick=thk
;  oplot,sgr(0,gdsgr),rsgr(gdsgr),linestyle=5,co=co4,thick=thk
  oplot,cararr(gdcar,0,0),rcar(gdcar),linestyle=5,co=white
endif

; point of closest approach
oplot,arr(colind,0,0),r(colind,0),co=black,ps=2,symsize=1.7
oplot,arr(colind,0,0),r(colind,0),co=yellow,ps=2,symsize=0.8
oplot,arr(colind,1,0),r(colind,1),co=black,ps=2,symsize=1.7
oplot,arr(colind,1,0),r(colind,1),co=yellow,ps=2,symsize=0.8

; Y vs X (lower-left)
if not keyword_set(twoplot) then begin
!p.multi=[2,2,2]
psym8,sym8a,/fill,/square
plot,dist(5),/nodata,xtit='XGC',ytit='YGC',xr=[-70,70],yr=[-130,150],xs=1,ys=1
oplot,[0],[0],ps=1       ;origin
oplot,arr(gdind,0,2),arr(gdind,0,3),linestyle=0,co=co1,thick=thk
oplot,[arr(0,0,2)],[arr(0,0,3)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,2),arr(gdind,1,3),linestyle=0,co=co2,thick=thk2
oplot,[arr(0,1,2)],[arr(0,1,3)],linestyle=0,co=co2,thick=thk,ps=8    ;current position
;oplot,arr(gdind,2,2),arr(gdind,2,3),linestyle=0,co=co3,thick=thk
;oplot,[arr(0,2,2)],[arr(0,2,3)],linestyle=0,co=co3,thick=thk,ps=8    ;current position
;oplot,marr(gd,1),marr(gd,2),linestyle=2
;oplot,sgrarr(*,0,2),sgrarr(*,0,3),co=white,thick=0.2

psym8,sym8b
;oplot,arr(tind,0,2),arr(tind,0,3),co=co1,ps=8
;oplot,arr(tind,1,2),arr(tind,1,3),co=co2,ps=8
oplot,arr(t1ind,0,2),arr(t1ind,0,3),co=co1,ps=8,symsize=sym1
oplot,arr(t2ind,0,2),arr(t2ind,0,3),co=co1,ps=8,symsize=sym2
oplot,arr(t1ind,1,2),arr(t1ind,1,3),co=co2,ps=8,symsize=sym1
oplot,arr(t2ind,1,2),arr(t2ind,1,3),co=co2,ps=8,symsize=sym2
;oplot,arr(t1ind,2,2),arr(t1ind,2,3),co=co3,ps=8,symsize=sym1
;oplot,arr(t2ind,2,2),arr(t2ind,2,3),co=co3,ps=8,symsize=sym2

if keyword_set(over) then begin
;  oplot,lmc(1,gdlmc),lmc(2,gdlmc),linestyle=5,co=co3,thick=thk
;  oplot,sgr(1,gdsgr),sgr(2,gdsgr),linestyle=5,co=co4,thick=thk
  oplot,cararr(gdcar,0,2),cararr(gdcar,0,3),linestyle=5,co=white
endif

; point of closest approach
oplot,arr(colind,0,2),arr(colind,0,3),co=black,ps=2,symsize=1.7
oplot,arr(colind,0,2),arr(colind,0,3),co=yellow,ps=2,symsize=0.8
oplot,arr(colind,1,2),arr(colind,1,3),co=black,ps=2,symsize=1.7
oplot,arr(colind,1,2),arr(colind,1,3),co=yellow,ps=2,symsize=0.8

end

; Z vs Y  (lower-right)
!p.multi=[1,2,2]
if keyword_set(twoplot) then !p.multi=[1,2,0]
psym8,sym8a,/fill,/square
;plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=[-70,70],yr=[-70,70],xs=1,ys=1
plot,dist(5),/nodata,xtit='YGC',ytit='ZGC',xr=[-130,150],yr=[-130,130],xs=1,ys=1
oplot,[0],[0],ps=1       ;origin
oplot,arr(gdind,0,3),arr(gdind,0,4),linestyle=0,co=co1,thick=thk
oplot,[arr(0,0,3)],[arr(0,0,4)],linestyle=0,co=co1,thick=thk,ps=8    ;current position
oplot,arr(gdind,1,3),arr(gdind,1,4),linestyle=0,co=co2,thick=thk2
oplot,[arr(0,1,3)],[arr(0,1,4)],linestyle=0,co=co2,thick=thk,ps=8    ;current position
;oplot,arr(gdind,2,3),arr(gdind,2,4),linestyle=0,co=co3,thick=thk2
;oplot,[arr(0,2,3)],[arr(0,2,4)],linestyle=0,co=co3,thick=thk,ps=8    ;current position
;oplot,marr(gd,2),marr(gd,3),linestyle=2
;oplot,sgrarr(*,0,3),sgrarr(*,0,4),co=white,thick=0.2

psym8,sym8b
;oplot,arr(tind,0,3),arr(tind,0,4),co=co1,ps=8
;oplot,arr(tind,1,3),arr(tind,1,4),co=co2,ps=8
oplot,arr(t1ind,0,3),arr(t1ind,0,4),co=co1,ps=8,symsize=sym1
oplot,arr(t2ind,0,3),arr(t2ind,0,4),co=co1,ps=8,symsize=sym2
oplot,arr(t1ind,1,3),arr(t1ind,1,4),co=co2,ps=8,symsize=sym1
oplot,arr(t2ind,1,3),arr(t2ind,1,4),co=co2,ps=8,symsize=sym2
;oplot,arr(t1ind,2,3),arr(t1ind,2,4),co=co3,ps=8,symsize=sym1
;oplot,arr(t2ind,2,3),arr(t2ind,2,4),co=co3,ps=8,symsize=sym2

if keyword_set(over) then begin
;  oplot,lmc(2,gdlmc),lmc(3,gdlmc),linestyle=5,co=co3,thick=thk
;  oplot,sgr(2,gdsgr),sgr(3,gdsgr),linestyle=5,co=co4,thick=thk
  oplot,cararr(gdcar,0,3),cararr(gdcar,0,4),linestyle=5,co=white
endif

; point of closest approach
oplot,arr(colind,0,3),arr(colind,0,4),co=black,ps=2,symsize=1.7
oplot,arr(colind,0,3),arr(colind,0,4),co=yellow,ps=2,symsize=0.8
oplot,arr(colind,1,3),arr(colind,1,4),co=black,ps=2,symsize=1.7
oplot,arr(colind,1,3),arr(colind,1,4),co=yellow,ps=2,symsize=0.8

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
