pro plotseparr4,save=save,color=color

;  This program plots the separr grid.

restore,'collgrid4_separr.dat'
;restore,'collgrid2_separr.dat'

; for collgrid2 and collgrid3
;vxarr = findgen(50)*5.-130.
;vyarr = findgen(40)*5.-320.
;vzarr = findgen(40)*5+60.

; We want to go three sigma in all directions
; From van der Marel
;mualph = 1.68  ; +/- 0.16 mas/yr
;mudelt = 0.34  ; +/- 0.16 mas/yr
mnmualph = 1.68
mnmudelt = 0.34
sigma = 4.
mualph0 = 1.68-sigma*0.16
mudelt0 = 0.34-sigma*0.16
ni = 50.
nj = 50.
dmualph = (2.*sigma*0.16)/(ni-1.)
dmudelt = (2.*sigma*0.16)/(nj-1.)

mualph = dindgen(ni)*dmualph + mualph0
mudelt = dindgen(nj)*dmudelt + mudelt0

smin=min(separr)
smax=max(separr)

oldseparr = separr
separr = dblarr(1,ni,nj)
separr(0,*,*) = oldseparr

nx = n_elements(separr(*,0,0))
ny = n_elements(separr(0,*,0))
nz = n_elements(separr(0,0,*))

window,xsize=640*1.2,ysize=512*1.2
;window,xsize=640,ysize=512    ;default

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

for i=0,nx-1 do begin

  ;if keyword_set(save) then ps_open,color=color

  imin = min(separr(i,*,*))
  imax = max(separr(i,*,*))
  strimin = stringize(imin,ndec=1)
  strimax = stringize(imax,ndec=1)

  minl = minloc(separr(i,*,*))   ; first goes through all x, then increment y
  zmin = long(minl/nz)
  ymin = minl-zmin*nz

  tit = '  (min='+strimin+', max='+strimax+' kpc)'
  xt = '!4l!da!n!3 cos(!4d!3) mas yr!u-1!n'
  yt = '!4l!dd!n!3 mas yr!u-1!n'
  charsize=1.5
  image = reform(separr(i,*,*))
;  image = max(image)+min(image)-image   ; flip it
  xs = mualph
  ys = mudelt
  ;display,reform(separr(i,*,*)),vyarr,vzarr,min=smin,max=smax,tit=tit,xtit=xt,ytit=yt
  ;wait,0.1
  ;plot,[vyarr(ymin)],[vzarr(zmin)],ps=1,co=red,/noerase,xstyle=5,ystyle=5,$
  ;  xr=[min(vyarr)-2.5,max(vyarr)+2.5],yr=[min(vzarr)-2.5,max(vzarr)+2.5]

  interp = 1.
  minval = smin
  maxvel = smax
  xtitle = xt
  ytitle = yt
  title = tit
  ;Device, Decomposed=0
  loadct,13

  ;  This is the meat of DISPLAY.pro
  ;---------------------------------
  Erase

   position = [0.1,0.1,0.95,0.75]
 ;  position = [0.1,0.1,0.9,0.9]
   xsize = (position(2) - position(0)) * !D.X_VSIZE
   ysize = (position(3) - position(1)) * !D.Y_VSIZE
   xstart = position(0) * !D.X_VSIZE
   ystart = position(1) * !D.Y_VSIZE
   im = CONGRID(image, xsize, ysize,/interp)
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

  Device, Decomposed=0   ; turns display color on
  TV, byte_im, /Device, xstart, ystart, $
        XSize=xsize, YSize=ysize
  Device, Decomposed=1   ; turns display color off

  xstyle = 1
  ystyle = 1
  Plot, [0,1], /NoErase, /NoData, XStyle=xstyle, YStyle=ystyle, $
                  /Device, Position=dev_pos, $
                  XRange=xrange, YRange=yrange, $
                  Title=title, XTitle=xtitle, YTitle=ytitle,$
                  charsize=charsize,xticks=xticks,yticks=yticks,$
                  color=framecolor     ;,xgridstyle=1,ygridstyle=1,ticklen=1  ; gridlines

   ; overplotting the colorbar
   colorbar,position=[0.1, 0.85, 0.95, 0.90],minrange=smin,maxrange=smax,$    ;/invert
            charsize=charsize

   ; overplot sigma boxes
   plot,[mnmualph], [mnmudelt],ps=1,$
       position=dev_pos,/device,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
   plot,mnmualph+[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+[1.,1.,-1.,-1.,1.]*0.16,$
       position=dev_pos,/device,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
   plot,mnmualph+2.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+2.*[1.,1.,-1.,-1.,1.]*0.16,$
       position=dev_pos,/device,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
   plot,mnmualph+3.*[-1.,1.,1.,-1.,-1.]*0.16, mnmudelt+3.*[1.,1.,-1.,-1.,1.]*0.16,$
       position=dev_pos,/device,/noerase,xstyle=5,ystyle=5,XRange=xrange, YRange=yrange
   xyouts,1.68+0.16+0.01, 0.34+0.16-0.03,'1!4r!3',charsize=charsize
   xyouts,1.68+2.*0.16+0.01, 0.34+2.*0.16-0.03,'2!4r!3',charsize=charsize
   xyouts,1.68+3.*0.16+0.01, 0.34+3.*0.16-0.03,'3!4r!3',charsize=charsize

  if keyword_set(save) then begin
    num = strtrim(long(i),2)
    if num lt 10 then num='0'+num
    if num lt 100 then num='0'+num
    file = 'plotseparr4.tiff'

    ; creating tiff image
    tiff = tvrd(true=1)
    tiff = reverse(tiff,3)
    write_tiff,file,tiff,1
  endif

  stop

 ; if keyword_set(save) then begin
 ;   ps_close
 ;   num = strtrim(long(i),2)
 ;   if num lt 10 then num='0'+num
 ;   if num lt 100 then num='0'+num
 ;   file = 'plotseparr3_'+num+'.ps'
 ;   spawn,'cp idl.ps animsep3/'+file
 ;   ;file = 'plotseparr2_'+num+'.ps'
 ;   ;spawn,'cp idl.ps animsep2/'+file
 ; endif

  ;stop

end

if keyword_set(save) then $
spawn,'/astro8/bin/convert -delay 20 animsep3/*.tiff animsep3/animsept3.gif'

;spawn,'ps2gif animsep3/plotseparr3_*ps'
;spawn,'mogrify -rotate -90 animsep3/plotseparr3_*gif'
;spawn,'gifmerge animsep3/plotseparr3_*gif > animsep3/animsep3.gif'
;;spawn,'ps2gif animsep2/plotseparr2_*ps'
;;spawn,'mogrify -rotate -90 animsep2/plotseparr2_*gif'
;;spawn,'gifmerge animsep2/plotseparr2_*gif > animsep2/animsep2.gif'

stop

end
