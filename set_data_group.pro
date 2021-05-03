PRO set_data_group, GROUP = group, $
                    TEST = test
                    

common _data
common _nhdist
common _nhobs
common _rxnh


;; AD test or JOINT (Fisher method)
if (n_elements(test) eq 0) then message, 'NO TEST STAT SPECIFIED.' else $
                                test = strupcase(test)

group = strupcase(group)

;; set variables for modeling
case group of 
    'WAC': begin 
        rxd = rxdet[where(iiwd,ndet)]
        e_rxd = e_rxdet[where(iiwd)]
        rxl = rxlim[where(iiwac)]
        frac_hix = total(iiwac and iihix)/total(iiwac)
        frac_lox = total(iiwac and iilox)/total(iiwac)
        end
    'SEC': begin
        rxd = rxdet[where(iisd,ndet)]
        e_rxd = e_rxdet[where(iisd)]
        rxl = rxlim[where(~iiwac)]
        frac_hix = total(~iiwac and iihix)/total(~iiwac)
        frac_lox = total(~iiwac and iilox)/total(~iiwac)
        end
    'ALL': begin
        rxd = rxdet[where(iiad,ndet)]
        e_rxd = e_rxdet[where(iiad)]
        rxl = rxlim
        frac_hix = total(iihix)/total(iihix or iilox)
        frac_lox = total(iilox)/total(iihix or iilox)
        end
    'WAC_HIX': begin
        rxd = rxdet[where(iiwd and iihix,ndet)]
        e_rxd = e_rxdet[where(iiwd and iihix)]
        rxl = rxlim[where(iiwac and iihix)]
        frac_hix = 1.
        frac_lox = 0.
        end
    'WAC_LOX': begin
        rxd = rxdet[where(iiwd and iilox,ndet)]
        e_rxd = e_rxdet[where(iiwd and iilox)]
        rxl = rxlim[where(iiwac and iilox)]
        frac_hix = 0.
        frac_lox = 1.
        end
    'OFFSET': begin
        rxd = rxdet[where(iiwd,ndet)]+0.3
        e_rxd = e_rxdet[where(iiwd)]
        rxl = rxlim[where(iiwac)]
        frac_hix = total(iiwac and iihix)/total(iiwac)
        frac_lox = total(iiwac and iilox)/total(iiwac)
        end
endcase
;; full number of sources, detected and non-detected
nsrc = n_elements(rxl)

sav_vars = ['GROUP','RXD','E_RXD','RXL','FRAC_HIX','FRAC_LOX','NSRC','TEST']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="select_group.sav"')


END

