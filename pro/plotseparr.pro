pro plotseparr,file,psave=psave,color=color,post=post,$
                single=single,sub=sub,top=top,$
                vcirc=vcirc,ds=ds,anim=anim,savefile=savefile,$
                colrange=colrange,rms=rms,stp=stp,$
                interp=interp,seprms=seprms,overobs=overobs,$
                top2=top2,charsize=charsize,charthick=charthick,$
                thick=thick,title=title,tit2=tit2,overbest=overbest

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
;    TOP        Sets the maximum for the first plot
;    TOP2       Sets the maximum for the second plot
;    VCIRC      The vcirc to be plotted, default vcirc=220 km/s
;    /ANIM      Create an animation with TIFFs, can use savefile
;    DS         Which particular halo softening parameter value to use
;    COLRANGE   Time range to use for minimum, default colrange=[-2.5,-1.0]
;    /STP       Stop at the end
;    /INTERP    Interpolates b/w the points (off by default)
;    /RMS       Plot the RMS array instead of SEPARR
;    /SEPRMS    Plot SEPARR on the left and RMS on the right
;    /OVEROBS   Overplot the individual LMC prop.mot. measurements
;    /OVERBEST  Overplot the proper motion value for the best-fit orbit
;    CHARSIZE   The character size
;    CHARTHICK  The thickness of the characters
;    THICK      The thickness of the lines
;    TITLE      The title of the first plot
;    TIT2       The title of the second plot

; Restoring the data file
; includes separr, tarr, separr3, tarr2
; mualpharr, mudeltarr, vcircarr
;if not keyword_set(file) then file = 'collgrid9_separr.dat'
if not keyword_set(file) then begin
  print,'Syntax -  plotseparr,file,psave=psave,color=color,post=post,'
  print,'              single=single,sub=sub,top=top,'
  print,'              vcirc=vcirc,ds=ds,anim=anim,savefile=savefile,'
  print,'              colrange=colrange,rms=rms,stp=stp,interp=interp'
  return
endif
restore,file

if not keyword_set(vcirc) then vcirc=220.
if not keyword_set(ds) then ds=13.
if not keyword_set(colrange) then colrange=[-2.5,-1.0]
orig_colrange = colrange
if max(colrange) gt 0. then colrange=-abs(colrange)
si = sort(colrange)
colrange = colrange(si)

mnmualph = 1.68
mnmudelt = 0.34

; What are we plotting?
if keyword_set(rms) then begin
  array1 = rmsarr
  single = 1
  if not keyword_set(top) then begin
    dum = array1
    b=where(dum gt 4*median(dum),nb)
    dum(b)=0.
    top=median(dum)+2.*stdev(dum)
  endif
endif else begin
  array1 = separr
  array2 = tarr
endelse
if keyword_set(seprms) then begin
  array1 = separr
  array2 = rmsarr
  if not keyword_set(top2) then begin
    dum = array2
    b = where(finite(dum) eq 0,nb)
    if nb gt 0 then dum(b)=0.
    med = median(dum) 
    b=where(dum gt 4*med,nb)
    dum(b)=0.
    top2=med+2.*stdev(dum)
  endif
endif

; Graphics keywords
if not keyword_set(thick) then $
  if keyword_set(psave) then thick=2.5 else thick=1.0
if not keyword_set(charthick) then $
  if keyword_set(psave) then charthick=2.5 else charthick=1.0
if not keyword_set(charsize) then $
  if keyword_set(psave) then charsize=1.7 else charsize=1.0

nj = n_elements(array1(*,0))
nk = n_elements(array1(0,*))

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

  array1 = array1(xmin:xmax,ymin:ymax,zmin:zmax)
  array2 = array2(xmin:xmax,ymin:ymax,zmin:zmax)
  ;tarr = tarr(xmin:xmax,ymin:ymax,zmin:zmax)
  dsarr = dsarr(xmin:xmax)
 ; vcircarr = vcircarr(xmin:xmax)
  mualpharr= mualpharr(ymin:ymax)
  mudeltarr= mudeltarr(zmin:zmax)
endif

; minimum values for the arrays, can set this manually
smin=min(array1)
smax=max(array1)
;smax = 70.
if keyword_set(top) then smax=top

tmin = min(array2)
tmax = max(array2)
if keyword_set(top2) then tmax=top2


;oldseparr = separr
;separr = dblarr(1,ni,nj)
;separr(0,*,*) = oldseparr

;nx = n_elements(separr(*,0,0))
ny = n_elements(array1(*,0))
nz = n_elements(array1(0,*))

;if not keyword_set(psave) then window,xsize=640*1.2,ysize=512*1.2
fac=1.3
if not keyword_set(psave) then begin
  if not keyword_set(single) then window,xsize=900.,ysize=512 else $   ;default
          window,xsize=640,ysize=512    ;default
endif
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

  if not keyword_set(savefile) then savefile=maketemp('plotseparr')
  if keyword_set(psave) then ps_open,savefile,color=color

  imin = min(array1(*,*))
  imax = max(array1(*,*))
  strimin = stringize(imin,ndec=1)
  strimax = stringize(imax,ndec=1)

  ;minl = minloc(separr(i,*,*))   ; first goes through all x, then increment y
  ;zmin = long(minl/nz)
  ;ymin = minl-zmin*nz

  ;alltit = '!3Ds = '+stringize(dsarr(i),ndec=2)+' kpc'
  alltit=''
  if not keyword_set(rms) then tit = 'Minimum Separation Distance ['+strimin+', '+strimax+']' $
    else tit = 'RMS ['+strimin+', '+strimax+']'
  ;tit = 'Minimum Separation Distance (min='+strimin+',max='+strimax+' kpc)'
  if keyword_set(title) then tit=title
  xt = '!4l!da!n!3 cos(!4d!3) mas yr!u-1!n'
  yt = '!4l!dd!n!3 mas yr!u-1!n'
  ;charsize=1.5
  ;if keyword_set(psave) and not keyword_set(single) then charsize=charsize*0.65
  image = reform(array1(*,*))
;  image = max(image)+min(image)-image   ; flip it
  xs = mualpharr
  ys = mudeltarr
  nxs = n_elements(xs)
  nys = n_elements(ys)
  ;display,reform(separr(i,*,*)),vyarr,vzarr,min=smin,max=smax,tit=tit,xtit=xt,ytit=yt
  ;wait,0.1
  ;plot,[vyarr(ymin)],[vzarr(zmin)],ps=1,co=red,/noerase,xstyle=5,ystyle=5,$
  ;  xr=[min(vyarr)-2.5,max(vyarr)+2.5],yr=[min(vzarr)-2.5,max(vzarr)+2.5]

  if n_elements(interp) eq 0 then interp = 0   ; no interpolation by default
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

   decomp = 1
   if keyword_set(psave) then decomp=0
   display,reform(array1(*,*)),xs,ys,pos=position,/noerase,decomp=decomp,$
           xtit=xtitle,ytit=ytitle,tit=tit,charsize=charsize,$
           XStyle=xstyle, YStyle=ystyle,min=smin,max=smax,$
           interp=interp,xrange=xrange,yrange=yrange,thick=thick,$
           charthick=charthick

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
            charsize=charsize,xthick=thick,ythick=thick,charthick=charthick

   ; overplot sigma boxes
   device=0
   normal=1
   ;plot,[mnmualph], [mnmudelt],ps=1,$
   ;    position=position,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
   ;plot,mnmualph+[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+[1.,1.,-1.,-1.,1.]*0.16,$
   ;    position=position,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
   ;plot,mnmualph+2.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+2.*[1.,1.,-1.,-1.,1.]*0.16,$
   ;    position=position,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
   ;plot,mnmualph+3.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+3.*[1.,1.,-1.,-1.,1.]*0.16,$
   ;    position=position,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange

   ; overploting the individual proper motion measurements, from Palma and van der Marel
   if keyword_set(overobs) then begin
     ;psym8,1.5
     ;mua = [1.20, 1.30, 1.94, 1.60, 2.00]  ; Jones et al 1994, Kroupa et al 1994, Kroupa & Bastian 1997
     ;mud = [0.26, 1.10, -0.14, 0.19, 0.4]  ; van Leeuwen & Evans 1998, Pedreros, Anguita & Maza 2002
     ;err = [0.28, 0.65, 0.33, 0.33, 0.2]   ; Kallivayalil et al. 2005
     mua = [1.30, 1.36, 1.94, 1.60, 1.83, 1.40, 2.03]  ; from van der Marel, plus Drake
     mud = [1.10, -0.16, -0.14, 0.19, 0.66, 0.38, 0.44]  ; Kroupa et al. 1994, Jones et al. 1994,
     err = [0.65, 0.28, 0.33, 0.2, 0.33, 0.32, 0.07]     ; Kroupa & Bastian 1997, (van Leeuwn & Evans 1998)
                                                   ; Pedreros et al. 2003, Drake et al. 2001, Kallivayalil et al. 2005
     if keyword_set(thick) then errthk = thick else errthk = 2.0

     psym8,1.0
     ploterr2, mua, mud, err, position=position,device=device,normal=normal,/noerase,$
         xstyle=5,ystyle=5,XRange=xrange, YRange=yrange, thick=errthk, ps=8
   endif ; /overobs
 
   ; overplotting the best fitting point
   if keyword_set(overbest) then begin
     bestind = array_indices(rmsarr,minloc(rmsarr))
     bmualph = mualpharr(bestind(0))
     bmudelt = mudeltarr(bestind(1))
     plot,[bmualph],[bmudelt],position=position,device=device,normal=normal,/noerase,$
           xs=5,ys=5,xrange=xrange,yrange=yrange,thick=thick,ps=1
   endif ;/overbest

 ;  if max(xrange) ge 2.17 and max(yrange) ge 0.79 then begin  ;make sure the plot is large enough
 ;  xyouts,1.68+0.16+0.01, 0.34+0.16-0.03,'1!4r!3',charsize=charsize
 ;  xyouts,1.68+2.*0.16+0.01, 0.34+2.*0.16-0.03,'2!4r!3',charsize=charsize
 ;  xyouts,1.68+3.*0.16+0.01, 0.34+3.*0.16-0.03,'3!4r!3',charsize=charsize
 ;  endif

  ; PLOTTING ARRAY2 ON THE RIGHT SIDE
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

   ;itmin = min(array2(i,*,*))
   ;itmax = max(array2(i,*,*))
   ;stritmin = stringize(itmin,ndec=1)
   ;stritmax = stringize(itmax,ndec=1)
   strtmin = stringize(tmin,ndec=1)
   strtmax = stringize(tmax,ndec=1)

   tit2='Time of Minimum Separation ['+strtmin+', '+strtmax+']'
   if not keyword_set(seprms) then tit2='Time of Minimum Separation ['+strtmin+', '+strtmax+']' $
          else tit2 = 'RMS ['+strtmin+', '+strtmax+']'
   if keyword_set(title2) then tit2=title2
   decomp = 1
   if keyword_set(psave) then decomp=0
   display,reform(array2(*,*)),xs,ys,pos=pos2,/noerase,decomp=decomp,$
           xtit=xtitle,ytit=ytitle,tit=tit2,charsize=charsize,$
           XStyle=xstyle, YStyle=ystyle,min=tmin,max=tmax,$
           interp=interp,xrange=xrange,yrange=yrange,thick=thick,$
           charthick=charthick

   ; overplotting the colorbar
   colpos2 = colpos
   colpos2(0) = colpos2(0)+0.5
   colpos2(2) = colpos2(2)+0.5
   colorbar,position=colpos2,minrange=tmin,maxrange=tmax,$    ;/invert
            charsize=charsize,xthick=thick,ythick=thick,charthick=charthick  ;,format='(F4.1)'

   ; overplot sigma boxes
   device=0
   normal=1
;   plot,[mnmualph], [mnmudelt],ps=1,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+2.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+2.*[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;   plot,mnmualph+3.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+3.*[1.,1.,-1.,-1.,1.]*0.16,$
;       position=dev_pos2,device=device,normal=normal,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
;
   ; overploting the individual proper motion measurements, from Palma
   if keyword_set(overobs) then begin
     ;psym8,1.5
     ;mua = [1.20, 1.30, 1.94, 1.60, 2.00]  ; Jones et al 1994, Kroupa et al 1994, Kroupa & Bastian 1997
     ;mud = [0.40, 1.10, -0.14, 0.19, 0.4]  ; van Leeuwen & Evans 1998, Pedreros, Anguita & Maza 2002
     ;err = [0.28, 0.65, 0.33, 0.33, 0.2]   
     mua = [1.30, 1.36, 1.94, 1.60, 1.83, 1.40, 2.03]  ; from van der Marel, plus Drake
     mud = [1.10, -0.16, -0.14, 0.19, 0.66, 0.38, 0.44]   ; Kallivayalil et al. 2005
     err = [0.65, 0.28, 0.33, 0.2, 0.33, 0.32, 0.07]

     psym8,1.0
     ploterr2, mua, mud, err, position=pos2,device=device,normal=normal,/noerase,$
          xstyle=5,ystyle=5,XRange=xrange, YRange=yrange, thick=errthk, ps=8

   endif; /overobs

   ; overplotting the best fitting point
   if keyword_set(overbest) then begin
     bestind = array_indices(rmsarr,minloc(rmsarr))
     bmualph = mualpharr(bestind(0))
     bmudelt = mudeltarr(bestind(1))
     plot,[bmualph],[bmudelt],position=pos2,device=device,normal=normal,/noerase,$
           xs=5,ys=5,xrange=xrange,yrange=yrange,thick=thick,ps=1
   endif ; /overbest

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
    ext = 'plotseparr'
    if keyword_set(savefile) then ext = savefile
    file = ext+'_'+num+'.tiff'

    base = '/home/frosty/dln5q/'
    dum = findfile(base)
    if dum(0) eq '' then base = '/Users/davidnidever/'

    dir = base+'nbody/plotseparr'+maketemp()
    if keyword_set(savefile) then dir = base+'nbody/'+savefile+'/'
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

  ;if not keyword_set(anim) then stop

  ;breaking the loop
  if not keyword_set(anim) then flag=1
  if keyword_set(anim) and i eq nj-1 then flag=1
  
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

if keyword_set(stp) then stop

colrange = orig_colrange

end
