FUNCTION rx2nh, in_arr, $
                IN_RX = in_rx, $
                RX_OUT = rx_out, $
                SCAT = scat
                

if (n_elements(model) eq 0) then begin
    print, '======================================'
    print, 'NO MODEL SET. RUNNING CARROLL+21 MODEL'
    print, '======================================'
    model = 'POWER'
endif
out_arr = dblarr(n_elements(in_arr))-9999.

nh = [21.:25.:0.25];5:0.25]
;; choose XSPEC model
if (n_elements(in_rx) eq 0) then begin
    rx = [0.0000000,-0.0033829937,-0.0093284138,-0.019683618,-0.037441140, $
          -0.065567324,-0.11028862,-0.17726416,-0.27116512,-0.39944674, $
          -0.58061745,-0.85145044,-1.2666158,-1.8337063,-2.2191357, $
          -2.2841807,-2.2952262];,-2.2993812,-2.3002889]
;    rx = [0.0000000,-0.0033660974,-0.0092815038,-0.019583447,-0.037246642, $
;          -0.065215346,-0.10966440,-0.17617641,-0.26929564,-0.39619911, $
;          -0.57459371,-0.83845970,-1.2304071,-1.7084625,-1.9592897, $
;          -1.9936579,-1.9992526];,-2.0013390,-2.0017935]
endif else rx = in_rx

;; Calculate either NH->L/L -OR- L/L->NH
if keyword_set(rx_out) then begin
    ;; Calculate Lx/Lx(Lir)
    ihi = where(in_arr gt nh[-1],nhi)
    if (nhi gt 0.) then out_arr[ihi] = rx[-1]
    ilo = where(in_arr lt nh[0] and in_arr ne -9999.,nlo)
    if (nlo gt 0.) then out_arr[ilo] = rx[0]
    igd = where(in_arr ge nh[0] and in_arr le nh[-1],ngd)
    if (ngd gt 0.) then out_arr[igd] = interpol(rx,nh,in_arr[igd])
endif else begin
    ;; Calculate NH
    ihi = where(in_arr gt rx[0],nhi)
    if (nhi gt 0.) then out_arr[ihi] = nh[0]
    ilo = where(in_arr lt rx[-1] and in_arr ne -9999.,nlo)
    if (nlo gt 0.) then out_arr[ilo] = nh[-1]
    igd = where(in_arr ge rx[-1] and in_arr le rx[0],ngd)
    if (ngd gt 0.) then out_arr[igd] = interpol(nh,rx,in_arr[igd])
endelse

;; add scatter
if (n_elements(scat) gt 0) then out_arr += randomn(seed,n_elements(out_arr))*scat

return, out_arr


END





