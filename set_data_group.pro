PRO set_data_group, GROUP = group
                    

common _data
common _nhdist
common _nhobs
common _rxnh


group = strupcase(group)

;; set variables for modeling
;ii = sdst_det gt 60. or sdst_non gt 60.
case group of 
    'WAC': begin 
        rxd = rxdet[where(iiwd,ndet)]
        e_rxd = e_rxdet[where(iiwd)]
        rxdn = dblarr(n_elements(where(iiwn,nnon)))-9999.
        rxl = rxlim[where(iiwac)]
        frac_hix = total(iiwac and iihix)/total(iiwac)
        frac_lox = total(iiwac and iilox)/total(iiwac)
        end
    'SEC': begin
        rxd = rxdet[where(iisd,ndet)]
        e_rxd = e_rxdet[where(iisd)]
        rxdn = dblarr(n_elements(where(iisn,nnon)))-9999.
        rxl = rxlim[where(~iiwac)]
        frac_hix = total(~iiwac and iihix)/total(~iiwac)
        frac_lox = total(~iiwac and iilox)/total(~iiwac)
        end
    'ALL': begin
        rxd = rxdet[where(iiad,ndet)]
        e_rxd = e_rxdet[where(iiad)]
        rxdn = dblarr(n_elements(where(iian,nnon)))-9999.
        rxl = rxlim
        frac_hix = total(iihix)/total(iihix or iilox)
        frac_lox = total(iilox)/total(iihix or iilox)
        end
    'WAC_HIX': begin
        rxd = rxdet[where(iiwd and iihix,ndet)]
        e_rxd = e_rxdet[where(iiwd and iihix)]
        rxdn = dblarr(n_elements(where(iiwn and iihix,nnon)))-9999.
        rxl = rxlim[where(iiwac and iihix)]
        frac_hix = 1.
        frac_lox = 0.
        end
    'WAC_LOX': begin
        rxd = rxdet[where(iiwd and iilox,ndet)]
        e_rxd = e_rxdet[where(iiwd and iilox)]
        rxdn = dblarr(n_elements(where(iiwn and iilox,nnon)))-9999.
        rxl = rxlim[where(iiwac and iilox)]
        frac_hix = 0.
        frac_lox = 1.
        end
    'OFFSET': begin
        rxd = rxdet[where(iiwd,ndet)]+0.3
        e_rxd = e_rxdet[where(iiwd)]
        rxdn = dblarr(n_elements(where(iiwn,nnon)))-9999.
        rxl = rxlim[where(iiwac)]
        frac_hix = total(iiwac and iihix)/total(iiwac)
        frac_lox = total(iiwac and iilox)/total(iiwac)
        end
    'OFFAXIS': begin 
        rxd = rxdet[where(iiwd and sdst_det gt 60.,ndet)]
        e_rxd = e_rxdet[where(iiwd and sdst_det gt 60.)]
        rxdn = dblarr(n_elements(where(iiwn and sdst_non gt 60.,nnon)))-9999.
        rxl = rxlim[where(iiwac and (sdst_det gt 60. or sdst_non gt 60.))]
        frac_hix = total(iiwac and iihix and (sdst_det gt 60. or sdst_non gt 60.))/total(iiwac and (sdst_det gt 60. or sdst_non gt 60.))
        frac_lox = total(iiwac and iilox and (sdst_det gt 60. or sdst_non gt 60.))/total(iiwac and (sdst_det gt 60. or sdst_non gt 60.))
        end
    'WAC2': begin 
        rxd = rxdet[where(iiwd and iitype2,ndet)]
        e_rxd = e_rxdet[where(iiwd and iitype2)]
        rxdn = dblarr(n_elements(where(iiwn and iitype2,nnon)))-9999.
        rxl = rxlim[where(iiwac and iitype2)]
        frac_hix = total(iiwac  and iitype2 and iihix)/total(iiwac)
        frac_lox = total(iiwac  and iitype2 and iilox)/total(iiwac)
        end
endcase
;; full number of sources, detected and non-detected
nsrc = n_elements(rxl)
;; data detection fraction
ddetf = 1.*ndet/nsrc

sav_vars = ['GROUP','RXD','E_RXD','RXDN','RXL','FRAC_HIX','FRAC_LOX','NSRC','NDET','NNON','DDETF']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="select_group.sav"')


END

