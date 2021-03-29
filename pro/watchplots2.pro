pro watchplots2

; This program runs run_coll.pro with a number of
; different parameter settings in order to get the
; "correct" orbit.

; lower limits
psi0 = 70.
vperp0 = 236.  ;230.
; range
rpsi = 10.     ;40.
rvperp = 10.   ;20
; number of steps
npsi = 11.     ;10.
nvperp = 11.   ;10.
; step sizes
dpsi = rpsi/(npsi-1.)         ;1.+rpsi/nspi
dvperp = rvperp/(nvperp-1.)   ;1.+rvperp/nvperp

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


; run plotcoll,/save

; save tcoll.out and tcoll.diag to other file with parameter
; values in name


stop

end
