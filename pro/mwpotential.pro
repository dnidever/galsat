pro mwpotential

; Plots David Law's MW potential

G = 4.3e-6
Mdisk = 1e11
Msph = 3.4e10
a = 6.5   ;kpc
b = 0.26  ;kpc
c = 0.7   ;kpc
vhalo = 121.
q = 0.9
d = 13.
z = 0.
rho = scale_vector(findgen(1000),0.1,100.)
r = sqrt(rho^2.+z^2.)

phi_disk = -G*Mdisk/sqrt((rho^2.) + (a+sqrt((z^2.) + (b^2.)))^2.)
phi_sph = -G*Msph/(r+c)
phi_halo = (vhalo^2.)*alog((rho^2.)+((z^2.)/((q^2.))+(d^2.)))

phi_tot = phi_disk + phi_sph + phi_halo

dr = r(1)-r(0)
acc = slope(phi_tot)/dr

xr=minmax(r)
yr=[-1e5,2e5]
xtit='Radius (kpc)'
ytit='Potential'
tit='Milky Way Gravitational Potential'

plot,r,phi_tot,xr=xr,yr=yr,xs=1,ys=1,xtit=xtit,ytit=ytit,tit=tit,/ylog
oplot,r,r*0
oplot,r,phi_disk,co=50
oplot,r,phi_sph,co=100
oplot,r,phi_halo,co=200


stop

end