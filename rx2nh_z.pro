FUNCTION rx2nh_z, in_arr, $
                  in_z, $
                  RX_OUT = rx_out, $
                  SCAT = scat
                

dir = '/Users/ccarroll/Research/projects/model_nhdist/workspace/data_prep/rx_grid.sav'
restore,dir
loc = value_locate(row_z,in_z)
out_arr = dblarr(n_elements(in_arr))

;; cut down to NH = 25
inh = where(col_nh le 25.)
col_nh = col_nh[inh]
rx_grid = rx_grid[inh,*]

;; Calculate either NH->L/L -OR- L/L->NH
if keyword_set(rx_out) then begin
    ;; Calculate Lx/Lx(Lir)
    ihi = where(in_arr gt col_nh[-1],nhi)
    if (nhi gt 0) then out_arr[ihi] = rx_grid[-1,loc[ihi]]
    ilo = where(in_arr lt col_nh[0] and in_arr ne -9999.,nlo)
    if (nlo gt 0) then out_arr[ilo] = rx_grid[0,loc[ilo]]
    igd = where(in_arr ge col_nh[0] and in_arr le col_nh[-1],ngd)
    if (ngd gt 0) then for i = 0,ngd-1 do out_arr[igd[i]] = interpol(rx_grid[*,loc[igd[i]]],col_nh,in_arr[igd[i]])
endif else begin
    ;; Calculate NH
    ihi = where(in_arr gt 0.,nhi)
    if (nhi gt 0) then out_arr[ihi] = col_nh[0]
    ilo = where(in_arr lt min(rx_grid),nlo)
    if (nlo gt 0) then out_arr[ilo] = col_nh[-1]
    igd = where(in_arr lt 0. and in_arr ge min(rx_grid),ngd)
    if (ngd gt 0) then for i = 0,ngd-1 do out_arr[igd[i]] = interpol(col_nh,rx_grid[*,loc[igd[i]]],in_arr[igd[i]])
endelse

;; add scatter
if (n_elements(scat) gt 0) then out_arr += randomn(seed,n_elements(out_arr))*scat

return, out_arr


END





