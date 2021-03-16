PRO set_data_group, GROUP = group


common _data
common _nhobs
common _rxnh

group = strupcase(group)

;; set variables for modeling
case group of 
    'WAC': begin 
        rxd = rxdet[where(iiwd,ndet)]
        e_rxd = e_rxdet[where(iiwd)]
        rxl = rxlim[where(iiwac)]
        frac_hix = total(iiwac and iixh)/total(iiwac)
        frac_lox = total(iiwac and iixl)/total(iiwac)
        end
    'SEC': begin
        rxd = rxdet[where(iisd,ndet)]
        e_rxd = e_rxdet[where(iisd)]
        rxl = rxlim[where(~iiwac)]
        frac_hix = total(~iiwac and iixh)/total(~iiwac)
        frac_lox = total(~iiwac and iixl)/total(~iiwac)
        end
    'ALL': begin
        rxd = rxdet[where(iiad,ndet)]
        e_rxd = e_rxdet[where(iiad)]
        rxl = rxlim
        frac_hix = total(iixh)/total(iixh or iixl)
        frac_lox = total(iixl)/total(iixh or iixl)
        end
    'WAC_HIX': begin
        rxd = rxdet[where(iiwd and iixh,ndet)]
        e_rxd = e_rxdet[where(iiwd and iixh)]
        rxl = rxlim[where(iiwac and iixh)]
        frac_hix = 1.
        frac_lox = 0.
        end
    'WAC_LOX': begin
        rxd = rxdet[where(iiwd and iixl,ndet)]
        e_rxd = e_rxdet[where(iiwd and iixl)]
        rxl = rxlim[where(iiwac and iixl)]
        frac_hix = 0.
        frac_lox = 1.
        end
endcase
;; full number of sources, detected and non-detected
nsrc = n_elements(rxl)

sav_vars = ['GROUP','RXD','E_RXD','RXL','FRAC_HIX','FRAC_LOX','NSRC']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="select_group.sav"')


END

