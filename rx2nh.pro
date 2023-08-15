FUNCTION rx2nh, in_arr, $
                RX_OUT = rx_out, $
                SCAT = scat, $
                MCMC = mcmc
                

common _rxnh
;seed = 1001

out_arr = dblarr(n_elements(in_arr))-9999.

;; Calculate either NH->L/L -OR- L/L->NH
if keyword_set(rx_out) then begin
    ;; Calculate Lx/Lx(Lir)
    if keyword_set(mcmc) then begin
        ;; random index for FSCAT
        ifs = fix(randomn(seed,10000)*0.8+2.)
        ifs = (ifs[where(ifs gt 0 and ifs lt 3)])[0:n_elements(in_arr)-1]
        
        ihi = where(in_arr gt nh[-1],nhi)
        if (nhi gt 0.) then out_arr[ihi] = rx[-1,ifs[ihi]]
        ilo = where(in_arr lt nh[0] and in_arr ne -9999.,nlo)
        if (nlo gt 0.) then out_arr[ilo] = rx[0,ifs[ilo]]
        igd = where(in_arr ge nh[0] and in_arr le nh[-1],ngd)
        if (ngd gt 0.) then begin
            for i = 0,ngd-1 do out_arr[igd[i]] = interpol(rx[*,ifs[igd[i]]],nh,in_arr[igd[i]])
        endif
    endif else begin
        ihi = where(in_arr gt nh[-1],nhi)
        if (nhi gt 0.) then out_arr[ihi] = rx[-1]
        ilo = where(in_arr lt nh[0] and in_arr ne -9999.,nlo)
        if (nlo gt 0.) then out_arr[ilo] = rx[0]
        igd = where(in_arr ge nh[0] and in_arr le nh[-1],ngd)
        if (ngd gt 0.) then out_arr[igd] = interpol(rx,nh,in_arr[igd])
    endelse
endif else begin
    ;; Calculate NH
    if keyword_set(mcmc) then begin
        ;; random index for FSCAT
        ifs = fix(randomn(seed,10000)*0.8+2.)
        ifs = (ifs[where(ifs gt 0 and ifs lt 3)])[0:n_elements(in_arr)-1]
        
        ihi = where(in_arr gt rx[0,ifs],nhi)
        if (nhi gt 0.) then out_arr[ihi] = nh[0]
        ilo = where(in_arr lt rx[-1,ifs] and in_arr ne -9999.,nlo)
        if (nlo gt 0.) then out_arr[ilo] = nh[-1]
        igd = where(in_arr ge rx[-1,ifs] and in_arr le rx[0,ifs],ngd)
        if (ngd gt 0.) then begin
            for i = 0,ngd-1 do out_arr[igd[i]] = interpol(nh,rx[*,ifs[igd[i]]],in_arr[igd[i]])
        endif        
    endif else begin
        ihi = where(in_arr gt rx[0],nhi)
        if (nhi gt 0.) then out_arr[ihi] = nh[0]
        ilo = where(in_arr lt rx[-1] and in_arr ne -9999.,nlo)
        if (nlo gt 0.) then out_arr[ilo] = nh[-1]
        igd = where(in_arr ge rx[-1] and in_arr le rx[0],ngd)
        if (ngd gt 0.) then out_arr[igd] = interpol(nh,rx,in_arr[igd])
    endelse
endelse

;; add scatter
if (n_elements(scat) gt 0) then out_arr += randomn(seed,n_elements(out_arr))*scat

return, out_arr


END





