pro lsmc_orbit,input,output,nbody=nbody,movie=movie,df=df,R0=R0,vcirc=vcirc,$
             noplot=noplot,step=step,ostep=ostep,dstep=dstep,fstep=fstep,$
             dhalo=dhalo

;+
; This program runs an N-test-body simulation of a
; galaxy with a plummer distribution of stars that
; orbits the Milky Way.
;
; INPUTS:
;  input    Input parameters [t,mass,soft,x,y,z,vx,vy,vz]
;  =nbody   Number of bodies to use in N-test-body simulation
;  /movie   Show the movie
;  /df      Use dynamical friction
;  =R0      Galactocentric distance of sun
;  =Vcirc   Circular rotation velocity of the sun
;  /noplot  Don't plot anything
;
; OUTPUTS:
;  output   Output parameters 
;
; CALLS:
;  plummer   To get the distribution of the N-bodies
;  galsat    To run the N-test-body simulation
;  plotnbody To plot the results or show the movie
; 
;
; By David Nidever   March 2006
;-

if n_elements(dhalo) eq 0 then dhalo = 13.0
if n_elements(vcirc) eq 0 then vcirc = 220.
if n_elements(R0) eq 0 then R0 = 8.5
if n_elements(df) eq 0 then df = 1
if n_elements(nbody) eq 0 then nbody=1000

;######################################################
;# STEP 1:  The current parameters
;######################################################

tstart = 1.0

;###################################
;# For the LMC
;###################################
; From Kallivayalil et al.2006
; Vx = -86+/-12, Vy = -268+/-11, Vz = 252+/-16 km/s   in GSR
Vxlmc = -86.0
Vylmc = -268.0
vzlmc = 252.0

; Calculating the galactocentric xyz coordinates for the LMC
; Taken from Connors et al. 2006
;R0 = 8.5
distlmc = 49.43
gllmc = 280.46              ;current l,b
gblmc = -32.89
lbd2xyz,gllmc,gblmc,distlmc,xlmc,ylmc,zlmc,R0=R0,/noprin

lmass = 2d10
lsoft = 0.9*(lmass/1d9)^(1./3.)  ; rscale
;lsoft = 1.0

inplmc = [lmass,lsoft,xlmc,ylmc,zlmc,vxlmc,vylmc,vzlmc]


;###################################
;# For the SMC
;###################################
; From Kallivayalil et al.2007
; Vx = -87+/-48, Vy = -247+/-42, Vz = 149+/-37 km/s   in GSR
Vxsmc = -87.0
Vysmc = -247.0
Vzsmc = 149.0

; Calculating the galactocentric xyz coordinates for the LMC
; Taken from Connors et al. 2006
;R0 = 8.5
distsmc = 57.02
glsmc = 302.79             ;current l,b
gbsmc = -44.30
lbd2xyz,glsmc,gbsmc,distsmc,xsmc,ysmc,zsmc,R0=R0,/noprin

smass = 1e9   ; SMC
ssoft = 0.9*(smass/1d9)^(1./3.)  ; rscale

inpsmc = [smass,ssoft,xsmc,ysmc,zsmc,vxsmc,vysmc,vzsmc]


;######################################################
;# STEP 2:  Run single particles backwards to get
;#          initial condition
;######################################################

print,'Running the two single particles backwards'

; Making the intial array
input1 = fltarr(2,8)
input1(0,*) = inplmc
input1(1,*) = inpsmc

; Running GALSAT
GALSAT,input1,output1,step=0.01,ostep=0.01,dstep=tstart,tmin=tstart,prog='galsat6',$
       dhalo=dhalo,vcirc=vcirc,df=df

; Getting the initial conditions
nlast = n_elements(output1[*,0,0])-1
inplmc2 = reform(output1[nlast,0,1:8])
inpsmc2 = reform(output1[nlast,1,1:8])



;######################################################
;# STEP 3:  Make the galaxies
;######################################################

;nbody = 1000

; Making the LMC parameters
print,'Making the LMC'
MAKEGAL,inplmc2,lmcarr,nbody=nbody

; Making the SMC parameters
print,'Making the SMC'
MAKEGAL,inpsmc2,smcarr,nbody=nbody


;######################################################
;# STEP 4:  Run the full N-body simulation
;######################################################


; Combining the arrays
input2 = [lmcarr,smcarr]


t0 = systime(1)

; SETTING GALSAT PARAMETERS
if not keyword_set(step) then step = 0.001               ; step size control parameter 
if keyword_set(movie) then ostep=0.01 else ostep=tstart  ; output interval
dstep = tstart                                           ; diagnostic interval
if not keyword_set(fstep) then fstep = 0.001             ; fixed timestep

; RUNNING GALSAT
galsat,input2,output2,step=step,ostep=ostep,dstep=dstep,tmin=tstart,df=df,dhalo=dhalo,$
           minstep=minstep,maxstep=maxstep,fstep=fstep,/forward,vcirc=vcirc,prog='galsat6'

; Correct the time
output2[*,*,0] = output2[*,*,0]-tstart

print,systime(1)-t0,' seconds'


; PLOTTING THE LAST SNAP (CURRENT)
nsnap = n_elements(output2[*,0,0])
if not keyword_set(noplot) then plotnbody,arr=output2,/last,R0=R0

; SHOWING THE MOVIE
if keyword_set(movie) then begin

  ;; Use the napo at the end
  ;colarr = reform(output2[nsnap-1,*,9])

  ;; Making reasonable colors, bound stars are white
  ;colarr2 = (colarr gt 0)*(colarr*30.+80) + (colarr eq 0)*255.
  ;colarr2(0) = 250   ; satellite is red

  colarr2=[fltarr(nbody)+150,fltarr(nbody)+250]
  plotnbody,arr=output2,/movie,/notrail,ps=3,R0=R0,colarr=colarr2 ;,/dotfirst,colarr=colarr2
endif


stop

end
