pro plotseparr8c,file,psave=psave,color=color,post=post,$
                single=single,sub=sub,top=top,$
                vcirc=vcirc,ds=ds,anim=anim,savefile=savefile,$
                colrange=colrange

;  This program plots the separr grid.
;
;  INPUTS
;
;    FILE       The name of the separation array file
;                  default='collgrid6_separr.dat'
;    SAVEFILE   The filename for the postscript file
;
;  KEYWORDS
;
;    /PSAVE      Creates a postscript file
;    /SINGLE    Only plots the separation and not the time as well
;    /SUB       Only plot a sub-section of the grid, can also set
;               the limits (in pixels/indices) explicitly
;    /TOP       Sets the maximum for the distance scale
;    /VCIRC     The vcirc to be plotted, default vcirc=220 km/s
;    /ANIM      Create an animation with TIFFs, can use savefile
;    DS         Which particular halo softening parameter value to use
;    COLRANGE   Time range to use for minimum, default colrange=[-2.5,-1.0]

; Restoring the data file
; includes separr, tarr, separr3, tarr2
; mualpharr, mudeltarr, vcircarr
if not keyword_set(file) then file = 'collgrid8c_separr.dat
restore,file

if not keyword_set(vcirc) then vcirc=220.
if not keyword_set(ds) then ds=13.
if not keyword_set(colrange) then colrange=[-2.5,-1.0]
orig_colrange = colrange
if max(colrange) gt 0. then colrange=-abs(colrange)
si = sort(colrange)
colrange = colrange(si)

; interpolate the vcirc values
;oldseparr = separr
;separr = congrid(separr,21,n_elements(separr(0,*,0)),n_elements(separr(0,0,*)),/inter)
;stop
;
; We want to go three sigma in all directions
; From van der Marel
;mualph = 1.68  ; +/- 0.16 mas/yr
;mudelt = 0.34  ; +/- 0.16 mas/yr
mnmualph = 1.68
mnmudelt = 0.34
;sigma = 4.
;mualph0 = 1.68-sigma*0.16
;mudelt0 = 0.34-sigma*0.16
;ni = 21.   ;10
;nj = 30.
;nk = 30.
;dmualph = (2.*sigma*0.16)/(nj-1.)
;dmudelt = (2.*sigma*0.16)/(nk-1.)
;
;mualph = dindgen(nj)*dmualph + mualph0
;mudelt = dindgen(nk)*dmudelt + mudelt0
;
;vcirc0 = 190.
;dvcirc = 40./(ni-1)
;
;vcirc = dindgen(ni)*dvcirc + vcirc0

; Setting the time range (from colrange)
; 0 - (-1.0,0.0), 1 - (-2.5,-1.0), 2 - (-4.0,-2.5)
; 3 - (-5.5,-4.0), 4 - (-7.0,-5.5)
;colarr = [-7.0,-5.5,-4.0,-2.5,-1.0,0.0]

;separr = separr(*,0:1,*,*)
;tarr = tarr(*,0:1,*,*)

;colarr = [0.0,-1.0,-2.5,-4.0,-5.5,-7.0]
;locol = closest(colrange(0),colarr,ind=loind,/lower)
;hicol = closest(colrange(1),colarr,ind=hiind,/upper)
;nrange = loind-hiind
;oldseparr = separr
;oldtarr = tarr
;
;ni = n_elements(oldseparr(0,*,0,0))
;nj = n_elements(oldseparr(0,0,*,0))
;nk = n_elements(oldseparr(0,0,0,*))
nj = n_elements(separr(*,0))
nk = n_elements(separr(0,*))
;separr = reform(separr(hiind,*,*,*)) 
;tarr = reform(tarr(hiind,*,*,*)) 
;
;; need to find the minimum of various ranges
;if nrange gt 1 then begin
;  for nr=1,nrange-1 do begin
;    for i=0,ni-1 do begin
;      for j=0,nj-1 do begin
;        for k=0,nk-1 do begin
;          if oldseparr(nr+hiind,i,j,k) lt separr(i,j,k) then begin
;            separr(i,j,k) = oldseparr(nr+hiind,i,j,k)
;            tarr(i,j,k) = oldtarr(nr+hiind,i,j,k)
;          endif
;          ;if oldseparr(nr+hiind,i,j,k) lt separr(i,j,k) then stop
;          ;separr = reform( oldseparr(nr+loind,*,*,*) < oldseparr(nr+loind+1,*,*,*) )
;        end  ; k
;      end  ; j
;    end  ; i
;  end                              ; nr
;endif  ; nrange > 1

;stop

; using separr2, the narrow range
;oldseparr = separr
;oldtarr = tarr
;separr = separr2
;tarr = tarr2

; only plotting a subarray
;sub=1
if keyword_set(sub) then begin

  if n_elements(sub) eq 4 then begin
    xmin = sub(0)
    xmax = sub(1)
    ymin = sub(2)
    ymax = sub(3)
    zmin = sub(4)
    zmax = sub(5)
  endif else begin
    xmin = 0
    ;xmax = 15
    xmax = ni-1
    ymin = 0
    ymax = 31  ;30  ;69
    zmin = 0
    zmax = 45    ;69
  endelse

  separr = separr(xmin:xmax,ymin:ymax,zmin:zmax)
  tarr = tarr(xmin:xmax,ymin:ymax,zmin:zmax)
  dsarr = dsarr(xmin:xmax)
 ; vcircarr = vcircarr(xmin:xmax)
  mualpharr= mualpharr(ymin:ymax)
  mudeltarr= mudeltarr(zmin:zmax)
endif

; minimum values for the arrays, can set this manually
smin=min(separr)
smax=max(separr)
;smax = 70.
if keyword_set(top) then smax=top

tmin = min(tarr)
tmax = max(tarr)


;oldseparr = separr
;separr = dblarr(1,ni,nj)
;separr(0,*,*) = oldseparr

;nx = n_elements(separr(*,0,0))
ny = n_elements(separr(*,0))
nz = n_elements(separr(0,*))

;if not keyword_set(psave) then window,xsize=640*1.2,ysize=512*1.2
fac=1.3
if not keyword_set(psave) then window,xsize=900.,ysize=512    ;default
;if not keyword_set(psave) then window,xsize=1000.,ysize=650    ;default
;window,xsize=640*fac,ysize=512*fac    ;default

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

; setting for postscript
if keyword_set(color) then begin
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
  ;backg=white
  backg=yellow
endif

;white=16777215

; LOOP for animation
;FOR i=0,nx-1 DO BEGIN
i=0
flag=0
WHILE (flag ne 1) do begin

  ;if not keyword_set(anim) and keyword_set(ds) then begin
  ;  dsoft=closest(ds,dsarr,ind=i)
  ;endif

  if not keyword_set(savefile) then savefile='plotseparr8c'
  if keyword_set(psave) then ps_open,savefile,color=color

  imin = min(separr(*,*))
  imax = max(separr(*,*))
  strimin = stringize(imin,ndec=1)
  strimax = stringize(imax,ndec=1)

  ;minl = minloc(separr(i,*,*))   ; first goes through all x, then increment y
  ;zmin = long(minl/nz)
  ;ymin = minl-zmin*nz

  ;alltit = '!3Ds = '+stringize(dsarr(i),ndec=2)+' kpc'
  alltit=''
  tit = 'Minimum Separation Distance ['+strimin+', '+strimax+']'
  ;tit = 'Minimum Separation Distance (min='+strimin+',max='+strimax+' kpc)'
  xt = '!4l!da!n!3 cos(!4d!3) mas yr!u-1!n'
  yt = '!4l!dd!n!3 mas yr!u-1!n'
  charsize=1.5
  if keyword_set(psave) and not keyword_set(single) then charsize=charsize*0.65
  image = reform(separr(*,*))
;  image = max(image)+min(image)-image   ; flip it
  xs = mualpharr
  ys = mudeltarr
  nxs = n_elements(xs)
  nys = n_elements(ys)
  ;display,reform(separr(i,*,*)),vyarr,vzarr,min=smin,max=smax,tit=tit,xtit=xt,ytit=yt
  ;wait,0.1
  ;plot,[vyarr(ymin)],[vzarr(zmin)],ps=1,co=red,/noerase,xstyle=5,ystyle=5,$
  ;  xr=[min(vyarr)-2.5,max(vyarr)+2.5],yr=[min(vzarr)-2.5,max(vzarr)+2.5]

  interp = 1.
  minval = smin
  maxval = smax
  xtitle = xt
  ytitle = yt
  title = tit
  ;Device, Decomposed=0
  loadct,13

  ;  This is the meat of DISPLAY.pro
  ;---------------------------------
  Erase


   position = [0.10,0.12,0.95,0.75]
   if not keyword_set(single) then begin
     position = [0.15,0.12,0.95,0.75]
     position(0) = position(0)/2.
     position(2) = position(2)/2.
   endif

 ;  position = [0.1,0.1,0.9,0.9]
   ; for X window
   xsize = (position(2) - position(0)) * !D.X_VSIZE
   ysize = (position(3) - position(1)) * !D.Y_VSIZE
   xstart = position(0) * !D.X_VSIZE
   ystart = position(1) * !D.Y_VSIZE
   device = 1
   normal = 0

   ; for postscript
   if keyword_set(psave) then begin
     position = [0.12,0.12,0.95,0.75]

     if not keyword_set(single) then begin
       position = [0.15,0.12,0.95,0.75]
       position(0) = position(0)/2.
       position(2) = position(2)/2.
     endif

     xstart = position(0)
     ystart = position(1)
     xsize = position(2)-position(0)
     ysize = position(3)-position(1)
     device = 0
     normal = 1
   endif

   interp = 0
   if not keyword_set(psave) then im = CONGRID(image, xsize, ysize,interp=interp) $
     else im = congrid(image,nxs*20 < 700,nys*20 < 700,interp=interp)
   dev_pos = [xstart,ystart,xstart+xsize,ystart+ysize]
   xrange = [min(xs),max(xs)]
   yrange = [min(ys),max(ys)]

;  im = ImgExp(image, xs, ys, xscale, yscale, xrange, yrange, $
;		Aspect=aspect, Interpolate=Keyword_Set(interp), $
;		MaskValue=maskvalue, Position=dev_pos, PS_Interp_Size=psis, $
;		No_Expand=Keyword_Set(no_expand))

;  sz = Size(im)
;  im_x_width = Float(sz(1))                   ;image width
;  im_y_width = Float(sz(2))                   ;image height

;  dev_x_width = dev_pos(2) - dev_pos(0) + 1
;  dev_y_width = dev_pos(3) - dev_pos(1) + 1


  if keyword_set(hist) then byte_im=hist_equal(im) else $
  byte_im = ImgScl(im, Min=minval, Max=maxval, Top=!D.Table_Size-2, $
			Log=log_scaling, Levels=l, MaskValue=maskvalue)

;  TV, byte_im, /Device, dev_pos(0), dev_pos(1), $
;	XSize=dev_pos(2)-dev_pos(0), YSize=dev_pos(3)-dev_pos(1)


  if not keyword_set(psave) then Device, Decomposed=0   ; turns display color on
  TV, byte_im, device=device, normal=normal, xstart, ystart, $
        XSize=xsize, YSize=ysize
  if not keyword_set(psave) then Device, Decomposed=1   ; turns display color off

  xstyle = 1
  ystyle = 1
  Plot, [0,1], /NoErase, /NoData, XStyle=xstyle, YStyle=ystyle, $
                  device=device, normal=normal, Position=dev_pos, $
                  XRange=xrange, YRange=yrange, $
                  Title=title, XTitle=xtitle, YTitle=ytitle,$
                  charsize=charsize,xticks=xticks,yticks=yticks,$
                  color=framecolor     ;,xgridstyle=1,ygridstyle=1,ticklen=1  ; gridlines

   ; overplotting the colorbar
   colpos = [0.1,0.87,0.95,0.92]
   ;colpos(0) = colpos(0)/2.
   ;colpos(2) = colpos(2)/2.
   colpos(0) = position(0)
   colpos(2) = position(2)
   

   if keyword_set(psave) then begin
     colpos = [0.12,0.87,0.95,0.92]
     colpos(0) = position(0)
     colpos(2) = position(2)
   endif
  ; if keyword_set(psave) then colpos = [0.1,0.85,0.95,0.90]
   colorbar,position=colpos,minrange=smin,maxrange=smax,$    ;/invert
            charsize=charsize

;   ; overplot sigma boxes
;   device=1
;   plot,[mnmualph], [mnmudelt],ps=1,$
;       position=dev_pos,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+2.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+2.*[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+3.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+3.*[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;
;   ; overploting the individual proper motion measurements, from Palma and van der Marel
;   psym8,1.5
;   ;mua = [1.20, 1.30, 1.94, 1.60, 2.00]  ; Jones et al 1994, Kroupa et al 1994, Kroupa & Bastian 1997
;   ;mud = [0.26, 1.10, -0.14, 0.19, 0.4]  ; van Leeuwen & Evans 1998, Pedreros, Anguita & Maza 2002
;   ;err = [0.28, 0.65, 0.33, 0.33, 0.2]
;   mua = [1.30, 1.36, 1.94, 1.60, 1.83, 1.40]  ; from van der Marel, plus Drake
;   mud = [1.10, -0.16, -0.14, 0.19, 0.66, 0.38]  
;   err = [0.65, 0.28, 0.33, 0.2, 0.33, 0.32]
;   errthk = 2.0
;
;   ploterr2, mua, mud, err, position=dev_pos,device=device,normal=normal,/noerase,$
;        xstyle=5,ystyle=5,XRange=xrange, YRange=yrange, thick=errthk
;
;   if max(xrange) ge 2.17 and max(yrange) ge 0.79 then begin  ;make sure the plot is large enough
;   xyouts,1.68+0.16+0.01, 0.34+0.16-0.03,'1!4r!3',charsize=charsize
;   xyouts,1.68+2.*0.16+0.01, 0.34+2.*0.16-0.03,'2!4r!3',charsize=charsize
;   xyouts,1.68+3.*0.16+0.01, 0.34+3.*0.16-0.03,'3!4r!3',charsize=charsize
;   endif

  ; Plotting the tarr on the right hand side
  if NOT keyword_set(single) then begin

   pos2 = position
   pos2(0) = pos2(0)+0.5
   pos2(2) = pos2(2)+0.5

   xsize2 = (pos2(2) - pos2(0)) * !D.X_VSIZE
   ysize2 = (pos2(3) - pos2(1)) * !D.Y_VSIZE
   xstart2 = pos2(0) * !D.X_VSIZE
   ystart2 = pos2(1) * !D.Y_VSIZE
   dev_pos2 = [xstart2,ystart2,xstart2+xsize2,ystart2+ysize2]
   device = 1
   normal = 0

   ; for postscript
   if keyword_set(psave) then begin
     xstart2 = pos2(0)
     ystart2 = pos2(1)
     xsize2 = pos2(2)-pos2(0)
     ysize2 = pos2(3)-pos2(1)
     dev_pos2 = [xstart2,ystart2,xstart2+xsize2,ystart2+ysize2]
     device = 0
     normal = 1
   endif

   itmin = min(tarr(i,*,*))
   itmax = max(tarr(i,*,*))
   stritmin = stringize(itmin,ndec=1)
   stritmax = stringize(itmax,ndec=1)

   tit2='Time of Minimum Separation ['+stritmin+', '+stritmax+']'
   decomp = 1
   if keyword_set(psave) then decomp=0
   display,reform(tarr(*,*)),xs,ys,pos=pos2,/noerase,decomp=decomp,$
           xtit=xtitle,ytit=ytitle,tit=tit2,charsize=charsize,$
           XStyle=xstyle, YStyle=ystyle,min=itmin,max=itmax,$
           interp=interp,xrange=xrange,yrange=yrange
   ; overplotting the colorbar
   colpos2 = colpos
   colpos2(0) = colpos2(0)+0.5
   colpos2(2) = colpos2(2)+0.5
   colorbar,position=colpos2,minrange=tmin,maxrange=tmax,$    ;/invert
            charsize=charsize,format='(F4.1)'

   ; overplot sigma boxes
;   device=1
;   plot,[mnmualph], [mnmudelt],ps=1,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+2.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+2.*[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+3.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+3.*[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;
;   ; overploting the individual proper motion measurements, from Palma
;   psym8,1.5
;   ;mua = [1.20, 1.30, 1.94, 1.60, 2.00]  ; Jones et al 1994, Kroupa et al 1994, Kroupa & Bastian 1997
;   ;mud = [0.40, 1.10, -0.14, 0.19, 0.4]  ; van Leeuwen & Evans 1998, Pedreros, Anguita & Maza 2002
;   ;err = [0.28, 0.65, 0.33, 0.33, 0.2]
;   mua = [1.30, 1.36, 1.94, 1.60, 1.83, 1.40]  ; from van der Marel, plus Drake
;   mud = [1.10, -0.16, -0.14, 0.19, 0.66, 0.38]  
;   err = [0.65, 0.28, 0.33, 0.2, 0.33, 0.32]
;   
;
;   ploterr2, mua, mud, err, position=dev_pos2,device=device,normal=normal,/noerase,$
;        xstyle=5,ystyle=5,XRange=xrange, YRange=yrange, thick=errthk
;
;   plot,[1.20], [0.26],ps=8,$       ;Jones et al. 1994  (+/-0.28)
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,[1.30], [1.10],ps=8,$       ;Kroupa et al. 1994  (+/-0.65)
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,[1.94], [-0.14],ps=8,$      ;Kroupa & Bastian 1997  (+/-0.33)
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,[1.60], [0.19],ps=8,$       ;van Leeuwen & Evans 1998  (+/-0.33)
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   ;plot,[1.70], [2.8],ps=8,$        ;Anguita 1998, Anguita, Loyola & Pedreros 2000  (+/-0.2)
;   ;    position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,[2.00], [0.4],ps=8,$        ;Pedreros, Anguita, & Maza 2002  (+/-0.2)
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   ;plot,[1.96], [0.08],ps=8,$        ;Drake et al. (2002), from average and other points
;   ;    position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;
;   if max(xrange) ge 2.17 and max(yrange) ge 0.79 then begin   ; make sure the plot is large enough
;   xyouts,1.68+0.16+0.01, 0.34+0.16-0.03,'1!4r!3',charsize=charsize
;   xyouts,1.68+2.*0.16+0.01, 0.34+2.*0.16-0.03,'2!4r!3',charsize=charsize
;   xyouts,1.68+3.*0.16+0.01, 0.34+3.*0.16-0.03,'3!4r!3',charsize=charsize
;   endif

  endif  ; not single

   ;plotting the overall title
   xyouts,0.5,0.95,alltit,charsize=charsize*1.5,/normal,alignment=0.5
 ;  xyouts,0.4,0.95,alltit,charsize=charsize*1.5,/normal

  if keyword_set(anim) then begin
    num = strtrim(long(i),2)
    if num lt 10 then num='0'+num
    if num lt 100 then num='0'+num
    ext = 'plotseparr8c'
    if keyword_set(savefile) then ext = savefile
    file = ext+'_'+num+'.tiff'

    dir = '/home/frosty/dln5q/nbody/animsep7/'
    if keyword_set(savefile) then dir = '/home/frosty/dln5q/nbody/'+savefile+'/'
    dum = findfile(dir)
    if dum(0) eq '' then spawn,'mkdir '+dir

    ; creating tiff image
    tiff = tvrd(true=1)
    tiff = reverse(tiff,3)
    write_tiff,dir+file,tiff,1
  endif

  ;stop

  if keyword_set(psave) then begin
    ps_close
 ;   num = strtrim(long(i),2)
 ;   if num lt 10 then num='0'+num
 ;   if num lt 100 then num='0'+num
 ;   file = 'plotseparr6_'+num+'.ps'
 ; ;  spawn,'cp idl.ps animsep3/'+file
 ;   ;file = 'plotseparr6_'+num+'.ps'
 ;   ;spawn,'cp idl.ps animsep2/'+file
  endif

  if not keyword_set(anim) then stop

  ;breaking the loop
  if not keyword_set(anim) then flag=1
  if keyword_set(anim) and i eq nx-1 then flag=1
  
  ;incrementing
  i=i+1

END  ; anim loop

;if keyword_set(anim) then $
;spawn,'/astro8/bin/convert -delay 20 animsep5/*.tiff animsep3/animsept3.gif'

;spawn,'ps2gif animsep3/plotseparr3_*ps'
;spawn,'mogrify -rotate -90 animsep3/plotseparr3_*gif'
;spawn,'gifmerge animsep3/plotseparr3_*gif > animsep3/animsep3.gif'
;;spawn,'ps2gif animsep2/plotseparr2_*ps'
;;spawn,'mogrify -rotate -90 animsep2/plotseparr2_*gif'
;;spawn,'gifmerge animsep2/plotseparr2_*gif > animsep2/animsep2.gif'

stop

colrange = orig_colrange

end
