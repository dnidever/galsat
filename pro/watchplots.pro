pro watchplots

; Used to look at the plots from collgrid.pro

; lower limits
psi0 = 70.
vperp0 = 230.
; range
rpsi = 40.
rvperp = 20
; number of steps
npsi = 10.
nvperp = 10.
; step sizes
dpsi = 1.+rpsi/npsi
dvperp = 1.+rvperp/nvperp

; run run_coll
for i=0,npsi-1 do begin
  psi = psi0 + i*dpsi

  for j=0,nvperp-1 do begin
    vperp = vperp0 + j*dvperp

    print,''
    print,'##########################################'
    print,'#  COLLGRID    PSI=',stringize(psi,ndec=2),'   VPERP=',stringize(vperp,ndec=2),'  #'
    print,'##########################################'
    print,''

    ; copying the files to correct location
    name = 'collgrid_'+stringize(psi,ndec=2)+'_'+stringize(vperp,ndec=2)
    spawn,'gv -swap collgrid/'+name+'.ps &'

    stop

  end  ;for j
end  ;for i

stop

end
