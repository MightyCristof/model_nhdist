FUNCTION rx2nh, in_arr, $
                RX_OUT = rx_out, $
                SCAT = scat
                

common _rx


out_arr = dblarr(n_elements(in_arr))-9999.

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





