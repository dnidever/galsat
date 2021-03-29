pro drive_collgrid,mlmc,soft,file

; This drives collgrid for a given Mlmc and soft.
; it will run three grids: 
;   1. The whole area
;   2. First minimum (based on grid of whole area)
;   3. Second minimum (based on grid of whole area)

nmlmc = n_elements(mlmc)
nsoft = n_elements(soft)
nfile = n_elements(file)

; Not enough parameters entered
if nmlmc eq 0 or nsoft eq 0 or nfile eq 0 then begin
  print,'Syntax - driver_collgrid,mlmc,soft,file'
  return
endif

prog = 'galsat6'

; Getting the base filename
orig_file = file
fil = file
lo = strpos(file,'.dat')
if lo ne -1 then begin
  fil = strmid(file,0,lo)
end

; Running it on the whole area
mualphr = [0.90, 2.5]
mudeltr = [-0.5, 1.15]
n=100
collgrid,mualphr,mudeltr,n,/df,soft=soft,prog=prog,savefile=fil+'a.dat',mlmc=mlmc

restore,fil+'a.dat'


; FIRST MINIMUM

; Looking for the minimum
;1.4803878     -0.30699628
dm = 0.3
mualphr = [1.48-dm,  1.48+dm]
mudeltr = [-0.31-dm, -0.31+dm]

dum = closest(mualphr(0),mualpharr,ind=ind1)
dum = closest(mualphr(1),mualpharr,ind=ind2)
dum = closest(mudeltr(0),mudeltarr,ind=ind3)
dum = closest(mudeltr(1),mudeltarr,ind=ind4)

x = mualpharr(ind1:ind2)
y = mudeltarr(ind3:ind4)
sep1 = separr(ind1:ind2,ind3:ind4)
mini=first_el(minloc(sep1))
minind = array_indices(sep1,mini)
minmualph = x(minind(0))
minmudelt = y(minind(1))

; Running galsat
dm2 = 0.2
mualphr = [minmualph-dm2, minmualph+dm2]
mudeltr = [minmudelt-dm2,minmudelt+dm2]
n=150
collgrid,mualphr,mudeltr,n,/df,soft=soft,prog=prog,savefile=fil+'b.dat',mlmc=mlmc


; SECOND MINIMUM

; Running it on the second minimum
;    1.1674898     0.098867463
dm = 0.3
mualphr = [1.17-dm, 1.17+dm]
mudeltr = [0.10-dm, 0.10+dm]

dum = closest(mualphr(0),mualpharr,ind=ind1)
dum = closest(mualphr(1),mualpharr,ind=ind2)
dum = closest(mudeltr(0),mudeltarr,ind=ind3)
dum = closest(mudeltr(1),mudeltarr,ind=ind4)

x = mualpharr(ind1:ind2)
y = mudeltarr(ind3:ind4)
sep1 = separr(ind1:ind2,ind3:ind4)
mini=first_el(minloc(sep1))
minind = array_indices(sep1,mini)
minmualph = x(minind(0))
minmudelt = y(minind(1))

; Running galsat
dm2 = 0.2
mualphr = [minmualph-dm2, minmualph+dm2]
mudeltr = [minmudelt-dm2,minmudelt+dm2]
n=150
collgrid,mualphr,mudeltr,n,/df,soft=soft,prog=prog,savefile=fil+'c.dat',mlmc=mlmc

;stop

end