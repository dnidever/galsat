;+
;
; GALSAT
;
; This program runs the c++ galsat program.  It's probably best to
; call it from a batch file. 
; NOTE: All velocities should be input in KM/S.  They will also be
;       output this way.  The velocity internal to the C++ galsat
;       program is kpc/Gyr, but galsat.pro and importnbody.pro 
;       automatically makes the right conversions.
;
; INPUT
;   input   Array of inputs [Nbodies, [mass,soft,x,y,z,Vx,Vy,Vz]]
;   step    Step size control parameter (default: step = 0.01)
;   dstep   Diagnostic interval (default: dstep = tmin)
;   ostep   Output interval (default: ostep = 0.01)
;   minstep Minimum timestep (positive), when using test particles
;   maxstep Maximum timestep (positive), when using test particles
;   fstep   Fixed timestep (positive), when using test particles
;   tmin    Minimum time to integrate to, ending time value (default tmin = 4.0)
;   prog    Name of the program to use (default prog='galsat8')
;   vcirc   Rotation velocity of the Milky Way disk (default vcirc = 220.0)
;   dhalo   Milky Way halo softening parameter (default ds = 13.0)
;   vhalo   Milky Way logarithmic halo Vhalo parameter (default vhalo
;              = 121.0).  Only for version galsat8 and later
;   /forward  Integrate forward instead of backwards (default: forward=0)
;   /mw     Do NOT use the Milky Way potential (default: mw=0)
;   /df     USE dynamical friction (default df=0)
;   /xtra   Output extra debugging information
;   /noprint Don't print anything to the screen
;
; OUTPUT
;   output  Array of galsat output. [Nsnap,Nbodies,[t,mass,soft,x,y,z,Vx,Vy,Vz]]
;   =diag   The diagnostics.  [Nsnap,[time, Etot, Ekin, Epot, Nsteps]]
;
; Created by David Nidever July 2005
;-

pro galsat,input,output,step=step,dstep=dstep,ostep=ostep,prog=prog,$
           vcirc=vcirc,dhalo=dhalo,vhalo=vhalo,tmin=tmin,forward=forward,mw=mw,df=df,$
           minstep=minstep,maxstep=maxstep,fstep=fstep,noprint=noprint,$
           xtra=xtra,diag=diag,stp=stp,nodelete=nodelete

; Bad Input Parameters
if n_params() eq 0 or n_elements(input) eq 0 then begin
  print,'syntax - galsat,input,output,step=step,dstep=dstep,ostep=ostep,'
  print,'                prog=prog,vcirc=vcirc,dhalo=dhalo,vhalo=vhalo,tmin=tmin,foward=forward'
  print,'                mw=mw,df=df,minstep=minstep,maxstep=maxstep,fstep=fstep'
  print,'                noprint=noprint,diag=diag'
  return
endif

; Where am I?
;dir = userdir()
dir = '~/'
bindir = dir+'nbody/'
;bindir = '/net/halo/dln5q/galsat/'
;bindir = '/Volumes/data/net/halo/dln5q/galsat/'

; 1 kpc = 3.0857E16 km, 1 Gyr = 3.1557E16 sec
; 1 kpc/Gyr = 3.0857/3.1557 = 0.977818
kms2kpcGyr = 3.0857d/3.1557d
kpcGyr2kms = 3.1557d/3.0857d

; galsat versions
vers = long(strmid(prog,6))
if strtrim(vers,2) eq '' then vers=1
if strmid(prog,0,6) ne 'galsat' then vers=0
if keyword_set(mw) eq 0 and vers lt 1 then message,'MW potential not available before galsat1'
if n_elements(df) gt 0 and vers lt 6 then message,'/df not available before galsat6'
if n_elements(forward) gt 0 and vers lt 4 then message,'/forward not available before galsat4'
if n_elements(vcirc) gt 0 and vers lt 2 then message,'=vcirc not available before galsat2'
if n_elements(dhalo) gt 0 and vers lt 3 then message,'=dhalo not available before galsat3'
if n_elements(vhalo) gt 0 and vers lt 8 then message,'=vhalo not available before galsat3'
if n_elements(minstep) gt 0 and vers lt 4 then message,'=minstep not available before galsat4'
if n_elements(maxstep) gt 0 and vers lt 4 then message,'=maxstep not available before galsat4'
if n_elements(fstep) gt 0 and vers lt 4 then message,'=fstep not available before galsat4'

; Defaults
if n_elements(prog) eq 0 then prog='galsat8'  ; 'galsat6'
if n_elements(tmin) eq 0 then tmin=4.0
if n_elements(step) eq 0 then step=0.01
if n_elements(dstep) eq 0 then dstep=tmin
if n_elements(ostep) eq 0 then ostep=0.01
if n_elements(vcirc) eq 0 then vcirc=220.
if n_elements(dhalo) eq 0 then dhalo=13.0
if n_elements(vhalo) eq 0 then vhalo=121.0
if n_elements(forward) eq 0 then forward=0
if n_elements(mw) eq 0 then mw=0
if n_elements(df) eq 0 then df=0

; # of Bodies
sz = size(input)
case sz[0] of
  1: nbody = 1
  2: nbody = sz(1)
endcase

; CONVERTING VELOCITIES (km/s -> kpc/Gyr)
input2 = input
; More than one particle
if nbody gt 1 then begin
  case sz[2] of
  ; mass, x, y, z, vx, vy, vz
  7: begin
    input2[*,4] = input2[*,4]*kms2kpcgyr
    input2[*,5] = input2[*,5]*kms2kpcgyr
    input2[*,6] = input2[*,6]*kms2kpcgyr
  end
  ; mass, soft, x, y, z, vx, vy, vz
  8: begin
    input2[*,5] = input2[*,5]*kms2kpcgyr
    input2[*,6] = input2[*,6]*kms2kpcgyr
    input2[*,7] = input2[*,7]*kms2kpcgyr
  end
  else: stop,'Format not understood'
  endcase
; One particle
endif else begin
  case sz[2] of
  ; mass, x, y, z, vx, vy, vz
  7: begin
    input2[4] = input2[4]*kms2kpcgyr
    input2[5] = input2[5]*kms2kpcgyr
    input2[6] = input2[6]*kms2kpcgyr
  end
  ; mass, soft, x, y, z, vx, vy, vz
  8: begin
    input2[5] = input2[5]*kms2kpcgyr
    input2[6] = input2[6]*kms2kpcgyr
    input2[7] = input2[7]*kms2kpcgyr
  end
  else: stop,'Format not understood'
  endcase
endelse
vcirc_kpcgyr = vcirc*kms2kpcgyr
vhalo_kpcgyr = vhalo*kms2kpcgyr


; Getting inputs ready
if tmin lt 0. then strtmin = stringize(-tmin,ndec=2) else strtmin=strtrim(tmin,2)
strstep = strtrim(step,2)
strdstep = strtrim(dstep,2)
strostep = strtrim(ostep,2)
;strvcirc = strtrim(vcirc,2)
strvcirc = strtrim(vcirc_kpcgyr,2)
strdhalo = strtrim(dhalo,2)
;strvhalo = strtrim(vhalo,2)
strvhalo = strtrim(vhalo_kpcgyr,2)
strnbody = strtrim(nbody,2)

; Setting flags
flags = ''
if keyword_set(forward) then flags = flags+' -f'
if keyword_set(mw) and vers ge 1 then flags = flags+' -m'
if keyword_set(df) then flags = flags+' -y'
if keyword_set(xtra) then flags = flags+' -x'

; CREATING INPUT FILE
tfile = maketemp('tgal')
openw,unit,/get_lun,tfile+'.in'
printf,unit,strnbody
printf,unit,'0'
if vers ge 2 then printf,unit,strvcirc    ; Vcirc
if vers ge 3 then printf,unit,strdhalo    ; dhalo
if vers ge 8 then printf,unit,strvhalo
if nbody gt 1 then $
    for i=0.,nbody-1 do printf,unit,strtrim(transpose(input2[i,*]),2)
;    for i=0.,nbody-1 do printf,unit,strtrim(input2(i,*),2)
if nbody eq 1 then printf,unit,strtrim(input2,2)
;printf,unit,'2.0e10 ',strsoft,' ',strxyz,' ',strvxvyvz
;printf,unit,'7.5e8 0.0 16.1556 2.27043 -5.88738  237.832 -42.3568 221.998'
close,unit
free_lun,unit

; Printing info
if not keyword_set(noprint) then begin

  ; Printing the param inputs
  print,'PARAM INPUTS'
  print,'-----------------------------------------------------'
  print,'Prog = ',prog
  print,'Nbodies = ',strnbody
  print,'Tmin   = ',strtmin,'   Gyrs'
  if vers gt 2 then print,'Vcirc  = ',strtrim(vcirc,2),'   km/s'
  if vers gt 3 then print,'Dhalo  = ',strdhalo,'   kpc'
  if vers ge 8 then print,'Vhalo  = ',strtrim(vhalo,2),'   km/s'
  print,'Step   = ',strstep
  print,'Ostep  = ',strostep,' Gyrs'
  print,'Dstep  = ',strdstep,' Gyrs'
  if keyword_set(forward) then print,'Integrating FORWARDS'
  if keyword_set(mw) then print,'NOT Using Milky Way Potential' else $
     print,'Using Milky Way Potential'
  if keyword_set(df) then print,'Using Dynamical Friction'
  print,'-----------------------------------------------------'
endif ; not noprint

; RUNNING GALSAT
print,' '
;print,'Running GALSAT program ... 3 sec.'
print,'Running GALSAT program ...'
;spawn,'\rm '+tfile+'.out'
cmd = '( '+bindir+prog+' -d '+strstep+' -o '+strostep+' -t '+strtmin+' -i -e '+strdstep
if keyword_set(minstep) then cmd = cmd+' -l '+strtrim(minstep,2)
if keyword_set(maxstep) then cmd = cmd+' -c '+strtrim(maxstep,2)
if keyword_set(fstep) then cmd = cmd+' -s '+strtrim(fstep,2)
cmd = cmd + flags+' < '+tfile+'.in > '+tfile+'.out ) > & '+tfile+'.diag'

;print,cmd

spawn,cmd,dum

;spawn,'./'+prog+' -d 0.01 -o 0.01 -t '+strtmin+' -i -e '+strtmin+' < '+tfile+'.in > '+tfile+'.out',dum

; GETTING THE OUTPUT
fileinfo = file_info(tfile+'.out')
fsize = fileinfo.size

; Everything in one file
if (fsize gt 0) then begin
  output = importnbody(tfile+'.out')

; In separate files
endif else begin
  galfiles = file_search('gal*.out')
  nfiles = n_elements(galfiles)

  d = importnbody(galfiles(0))
  sz = size(d)
  output = dblarr(nfiles,sz(2),sz(3))

  ; Restoring the files
  for i=0L,nfiles-1 do begin
    d = importnbody(galfiles(i))
    output(i,*,*) = d
  end

  ; Sort by time
  si = sort(output(*,0,0))
  output = output(si,*,*)

  stop

endelse

; Read in the diagnostic information
;  time  Etot  Ekin  Epot  Nsteps
readline,tfile+'.diag',diaglines,comment='#',count=ndiaglines
if keyword_set(xtra) then begin
  diaglines0 = diaglines
  bd = where(stregex(diaglines,'internal data for',/boolean) eq 1,nbd)
  if nbd gt 0 then remove,bd,diaglines
  brk = where(stregex(diaglines,'for debugging',/boolean),nbrk)
  lo = brk-1
  hi = [brk[1:*]-2,n_elements(diaglines)-1]
  diag = fltarr(nbrk,nbody,18)
  ; loop through outputs
  for i=0,nbrk-1 do begin
    temp = diaglines[lo[i]:hi[i]]
    sumline = temp[0]  ; t, etot, ekin, epot, nsteps
    sumlinearr = strsplit(sumline,' ',/extract)
    time = float(sumlinearr[0])
    temp = temp[2:*]  ; cut out header lines
    ; mass, soft, x, y, z, vx, vy, vz, accx, accy, accz, jerkx, jerky, jerkz
    ;   epot, ekin, etot
    tarr = strsplitter(temp,' ',/extract)
    tarr = float(transpose(tarr))
    diag[i,*,0] = time
    diag[i,*,1:*] = tarr
  endfor
endif

if keyword_set(stp) then stop

; REMOVING TEMPORARY FILES
if not keyword_set(nodelete) then file_delete,tfile+['.in','.out','.diag'],/allow

end
