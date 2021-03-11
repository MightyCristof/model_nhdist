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
        end
    'SEC': begin
        rxd = rxdet[where(iisd,ndet)]
        e_rxd = e_rxdet[where(iisd)]
        rxl = rxlim[where(~iiwac)]
        end
    'ALL': begin
        rxd = rxdet[where(iiad,ndet)]
        e_rxd = e_rxdet[where(iiad)]
        rxl = rxlim
        end
    'WAC_HIX': begin
        rxd = rxdet[where(iiwd and iihix,ndet)]
        e_rxd = e_rxdet[where(iiwd and iihix)]
        rxl = rxlim[where(iiwac and iihix)]
        end
    'WAC_LOX': begin
        rxd = rxdet[where(iiwd and iilox,ndet)]
        e_rxd = e_rxdet[where(iiwd and iilox)]
        rxl = rxlim[where(iiwac and iilox)]
        end
endcase
;; full number of sources, detected and non-detected
nsrc = n_elements(rxl)

sav_vars = ['GROUP','RXD','E_RXD','RXL','NSRC']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="select_group.sav"')


END

