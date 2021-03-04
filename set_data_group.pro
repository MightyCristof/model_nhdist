PRO set_data_group, MODE = mode


common _data
common _nhobs
common _rxnh

mode = strupcase(mode)

;; set variables for modeling
case mode of 
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
    'XHI': begin
        rxd = rxdet[where(iiwd and iixh,ndet)]
        e_rxd = e_rxdet[where(iiwd and iixh)]
        rxl = rxlim[where(iiwac and iixh)]
        end
    'XLO': begin
        rxd = rxdet[where(iiwd and iixl,ndet)]
        e_rxd = e_rxdet[where(iiwd and iixl)]
        rxl = rxlim[where(iiwac and iixl)]
        end
endcase
;; full number of sources, detected and non-detected
nsrc = n_elements(rxl)

sav_vars = ['MODE','RXD','E_RXD','RXL','NSRC']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="select_group.sav"')


END

