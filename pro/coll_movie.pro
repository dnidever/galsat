pro coll_movie,file=file,save=save,anim=anim

; This makes a movie of the collision
;
; some of this was copied from astr534/plot_orb.pro

deg2rad = !dpi/180.d
rad2deg = (180.d)/!dpi

file = 'coll_movie7'
arr = importnbody('best_coll.out')

; time, mass, x, y, z, vx, vy, vz
si = sort(arr(*,0,0))
lmc = reform(arr(si,0,*))
sgr = reform(arr(si,1,*))
tarr = reform(arr(si,0,0))
npts = n_elements(lmc(*,0))

separr = sqrt( (lmc(*,2)-sgr(*,2))^2. + (lmc(*,3)-sgr(*,3))^2. + (lmc(*,4)-sgr(*,4))^2. )

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
  backgr = 0.
  coarr = [green,orange,yellow,blue,purple,lred,lblue]

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
    lred=240
    ;backg=white
    backg=yellow
    coarr = [green,orange,yellow,blue,purple,lred,aqua]
  endif

  loadct,13
  !p.multi=0
  ;xd = max(lxarr)-min(lxarr)
  ;yd = max(lyarr)-min(lyarr)
  ;zd = max(lzarr)-min(lzarr)
  ;fac=0.5
  ;xr=[min(lxarr)-xd*fac,max(lxarr)+xd*fac]
  ;yr=[min(lyarr)-yd*fac,max(lyarr)+yd*fac]
  ;zr=[min(lzarr)-zd*fac*1.5,max(lzarr)+zd*fac]
  xr = [-100.,100.]
  yr = [-100.,100.]
  zr = [-100.,100.]

  ; creating the disk
  rdisk = 16.   ;25
  diskx = rdisk*sin(findgen(1000)/999.*2.*!dpi)
  disky = rdisk*cos(findgen(1000)/999.*2.*!dpi)
  diskz = disky*0.
  bd = where(abs(diskx) le 5. or abs(disky) le 5.,nbd)
  remove,bd,diskx,disky,diskz

  ; creating the bulge
  sphere,2.,bx,by,bz,/upper

  ; LMC and SGR spheres
  sphere,3.,lx,ly,lz
  sphere,1.5,sx,sy,sz

  sphere,0.5,sunx,suny,sunz

for i=0,npts-1 do begin

  if i/10 eq i/10. then print,i

  time = tarr(i)

  lxarr = lmc(i,2)
  lyarr = lmc(i,3)
  lzarr = lmc(i,4)

  sxarr = sgr(i,2)
  syarr = sgr(i,3)
  szarr = sgr(i,4)

  ;daz = 80./npts
  daz = 360./npts
  az0 = 360.  ;330.  ; 275.  ;-45.
  ax =  20.   ;45
  az = daz*i+az0   ;45.
  charsize = 1d-10   ; 2.5
  surface,dist(10),/nodata,az=az,ax=ax,/save,xr=xr,yr=yr,zr=zr,$
    ; xtit=' ',ytit=' ',ztit=' ',charsize=2.5,xs=1,ys=1,zs=1,$
    ; xtickname=[''],xticks=1,yticks=1,zticks=1,charthick=1d-10
    ; xtit='X',ytit='Y',ztit='Z',charsize=2.5,xs=1,ys=1,zs=1
     xtit='X',ytit='Y',ztit='Z',charsize=2.5,xs=4,ys=4,zs=4

  ; x-axes
  plots,[-100,100],[-100,-100],[-100,-100],/t3d
  plots,[-100,100],[100,100],[-100,-100],/t3d
  plots,[-100,100],[-100,-100],[100,100],/t3d
  plots,[-100,100],[100,100],[100,100],/t3d

  ; y-axes
  plots,[-100,-100],[-100,100],[-100,-100],/t3d
  plots,[100,100],[-100,100],[-100,-100],/t3d
  plots,[-100,-100],[-100,100],[100,100],/t3d
  plots,[100,100],[-100,100],[100,100],/t3d

  ; z-axes
  plots,[100,100],[100,100],[-100,100],/t3d
  plots,[-100,-100],[100,100],[-100,100],/t3d
  plots,[100,100],[-100,-100],[-100,100],/t3d
  plots,[-100,-100],[-100,-100],[-100,100],/t3d

  ;stop

  ;AXIS, /XAxis, 0, 100, -100, /T3D, charsize=charsize,xticks=1,yticks=1,zticks=1
  ;AXIS, /XAxis, 0, -100, 100, /T3D, charsize=charsize,xticks=1,yticks=1,zticks=1
  ;AXIS, /XAxis, 0, 100, 100, /T3D, charsize=charsize,xticks=1,yticks=1,zticks=1
  ;
  ;AXIS, /YAxis, 100, 0, -100, /T3D, charsize=charsize,yticks=1
  ;AXIS, /YAxis, -100, 0, 100, /T3D, charsize=charsize,yticks=1
  ;AXIS, /YAxis, 100, 0, 100, /T3D, charsize=charsize,yticks=1
  ;
  ;AXIS, /ZAxis, -100, -100, 0, /T3D, charsize=charsize,zticks=1
  ;AXIS, /ZAxis, 100, -100, 0, /T3D, charsize=charsize,zticks=1
  ;AXIS, /ZAxis, 100, 100, 0, /T3D, charsize=charsize,zticks=1

  ;stop

  ;overplot the milky way disk and bulge
  plots,diskx,disky,diskz,/t3d,color=blue
  polyfill,diskx,disky,diskz,/t3d,color=blue    ;,/line_fill
  ;plots,diskx,disky,diskz*0.+zr(0),/t3d,color=blue
  ;;polyfill,diskx,disky,diskz*0.+zr(0),/t3d,/line_fill,color=blue
  ;plots,diskx,disky*0.+yr(1),diskz,/t3d,color=blue
  ;plots,diskx*0.+xr(1),disky,diskz,/t3d,color=blue
  plots,bx,by,bz,color=yellow,/t3d  ; bulge

  ; sun's motion
  ; 220 km/s, 1 km/s = 1 pc/Myr = 1 kpc/Gyr
  ; 220 km/s = 220 pc/Myr = 220 kpc/Gyr
  ; period = 2*pi*8.5/220 
  period = 2.*!dpi*8.5/220.
  omega = 2.*!dpi/period
  sunposx = -8.5*cos(omega*time)    ; x=-8.5 at t=0
  sunposy = 8.5*sin(omega*time)    ; y=0    at t=0
  plots,sunx+sunposx,suny+sunposy,sunz,color=orange,/t3d  ; sun

  ; overplot LMC orbit
  plots,lxarr,lyarr,lzarr,/t3d,color=green
  ;plots,lxarr,lyarr,lzarr*0.+zr(0),/t3d,color=green
  ;plots,lxarr,lyarr*0.+yr(1),lzarr,/t3d,color=green
  ;plots,lxarr*0.+xr(1),lyarr,lzarr,/t3d,color=green
  ;plots,lxarr(npts-1)+lx,lyarr(npts-1)+ly,lzarr(npts-1)+lz,/t3d,color=green
  plots,lxarr+lx,lyarr+ly,lzarr+lz,/t3d,color=green
  plots,lmc(0:i,2),lmc(0:i,3),lmc(0:i,4),/t3d,color=green

  ; overplot SGR orbit
  plots,sxarr,syarr,szarr,/t3d,color=red
  ;plots,sxarr,syarr,szarr*0.+zr(0),/t3d,color=red
  ;plots,sxarr,syarr*0.+yr(1),szarr,/t3d,color=red
  ;plots,sxarr*0.+xr(1),syarr,szarr,/t3d,color=red
  ;plots,sxarr(npts-1)+sx,syarr(npts-1)+sy,szarr(npts-1)+sz,/t3d,color=red
  plots,sxarr+sx,syarr+sy,szarr+sz,/t3d,color=red
  plots,sgr(0:i,2),sgr(0:i,3),sgr(0:i,4),/t3d,color=red

  ; colliding, stars and BANG!, min(sep) = 2.6654452
  sep = sqrt( (lxarr-sxarr)^2. + (lyarr-syarr)^2. + (lzarr-szarr)^2. )
  if (sep lt 9.) then begin
  ;if (sep lt 6.) then begin
  ;if (sep lt 2.67) then begin
    cx = mean([lxarr,sxarr])
    cy = mean([lyarr,syarr])
    cz = mean([lzarr,szarr])   ;+20.
    rstar = 10.
    caz = cos(az*deg2rad)
    saz = sin(az*deg2rad)

    plots,[cx,cx],[cy,cy],[cz-rstar,cz+rstar],/t3d,color=yellow,thick=1.5
    plots,[cx-rstar*caz,cx+rstar*caz],[cy+rstar*saz,cy-rstar*saz],[cz,cz],/t3d,color=yellow,thick=1.5
    plots,[cx-0.5*rstar*caz,cx+0.5*rstar*caz],[cy+0.5*rstar*saz,cy-0.5*rstar*saz],$
          [cz-0.5*rstar,cz+0.5*rstar],/t3d,color=yellow,thick=1.5
    plots,[cx-0.5*rstar*caz,cx+0.5*rstar*caz],[cy+0.5*rstar*saz,cy-0.5*rstar*saz],$
          [cz+0.5*rstar,cz-0.5*rstar],/t3d,color=yellow,thick=1.5
    ;plots,[cx,cx],[cy,cy],[cz,cz],/t3d,color=yellow        ; diag 1
    ;plots,[cx,cx],[cy,cy],[cz,cz],/t3d,color=yellow        ; daig 2
    ;plots,cx,cy,cz,/t3d,color=yellow,ps=2,symsize=5

    ;xyouts,21.,40.,'BANG!',color=yellow,charsize=2.
    ;;xyouts,cx,cy,'BANG!',color=yellow,charsize=2.

   ;stop
  endif

  if keyword_set(anim) then begin
    num = strtrim(long(i),2)
    if num lt 10 then num='0'+num
    if num lt 100 then num='0'+num
    filename = file+'_'+num+'.tiff'

    ; creating tiff image
    tiff = tvrd(true=1)
    tiff = reverse(tiff,3)
    write_tiff,file+'/'+filename,tiff,1
  endif

  ;stop

end

stop

end
