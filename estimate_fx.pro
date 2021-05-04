FUNCTION estimate_fx, rx_model, $
                      iidet, $
                      CHA = cha, $
                      ITERATE = iterate


common _data
common _nhdist
common _nhobs
common _rxnh
common _group


;; use only Chandra sources to match Chandra X-ray stacking
if keyword_set(cha) then iicha = xnon eq 'CHA' else $
                         iicha = xnon ne ''
                         
;; iterate over each object, or run single case
if keyword_set(iterate) then niter = 100 else $
                             niter = 1

;; data for non-detected sample        
case group of 
    'WAC': begin 
        loglxir_non = loglxir[where(iiwn and iicha,nnon)]
        z_non = z[where(iiwn and iicha,nnon)]
        dl2_non = dl2[where(iiwn and iicha,nnon)]
        end
    'SEC': begin
        loglxir_non = loglxir[where(iisn and iicha,nnon)]
        z_non = z[where(iisn and iicha,nnon)]
        dl2_non = dl2[where(iisn and iicha,nnon)]
        end
    'ALL': begin
        loglxir_non = loglxir[where(iicha,nnon)]
        z_non = z[where(iicha,nnon)]
        dl2_non = dl2[where(iicha,nnon)]
        end
    'WAC_HIX': begin 
        loglxir_non = loglxir[where(iiwn and iixh and iicha,nnon)]
        z_non = z[where(iiwn and iixh and iicha,nnon)]
        dl2_non = dl2[where(iiwn and iixh and iicha,nnon)]
        end
    'WAC_LOX': begin 
        loglxir_non = loglxir[where(iiwn and iixl and iicha,nnon)]
        z_non = z[where(iiwn and iixl and iicha,nnon)]
        dl2_non = dl2[where(iiwn and iixl and iicha,nnon)]
        end
    'OFFSET': begin
        loglxir_non = loglxir[where(iiwn and iicha,nnon)]
        z_non = z[where(iiwn and iicha,nnon)]
        dl2_non = dl2[where(iiwn and iicha,nnon)]
        end
    else: message, message, 'NO INPUT MODE SET FOR FXEST: WAC/SEC/ALL'
endcase

;; prep hard and soft conversion arrays
iz = value_locate(zv,z_non)
;c_hard_non = c_hard[*,iz]
;c_soft_non = c_soft[*,iz]
c_hard_non = c_hard_fine[*,iz]
c_soft_non = c_soft_fine[*,iz]

;; reform to 2D array if needed
sz = size(rx_model)
nmod = product(sz[2:sz[0]])
rx_modd = reform(rx_model,sz[1],nmod)
ii_modd = reform(iidet,sz[1],nmod)
lin_vars = 'FX_'+['FULL','HARD','SOFT']
clin_vars = 'C'+lin_vars
lin_vars = [lin_vars,'E_'+lin_vars]
clin_vars = [clin_vars,'E_'+clin_vars]
for i = 0,n_elements(lin_vars)-1 do re = execute(lin_vars[i]+' = dblarr(nmod)')
for i = 0,n_elements(clin_vars)-1 do re = execute(clin_vars[i]+' = dblarr(nmod)')

;; sample from model RL for each non-detected observation (NNON)
for i = 0,nmod-1 do begin  
    ;; "non-detected" model sources
    inon = where(ii_modd[*,i] eq 0,nonct)
    if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
    rx_modn = (rx_modd[*,i])[inon]
    fullv = dblarr(nnon,niter)
    hardv = dblarr(nnon,niter)
    softv = dblarr(nnon,niter)    
    cfullv = dblarr(nnon,niter)
    chardv = dblarr(nnon,niter)
    csoftv = dblarr(nnon,niter)    
    for n = 0,niter-1 do begin
        ;; random draw from model to match number of non-detected sources
        rx_samp = mc_samp(rx_modn,nnon)
        ;; estimate flux
        lx_non = 10.^(rx_samp+loglxir_non)
        fullv[*,n] = lx_non/(4.*!const.pi*dl2_non)
        hardv[*,n] = fullv[*,n]
        softv[*,n] = fullv[*,n]
        ;; bounded indices
        ibnd = where(rx_samp ge min(rx),nbnd,complement=iunb,ncomplement=nunb)
        hardv[ibnd,n] *= c_hard_non[value_locate(rx_fine,rx_samp[ibnd]),ibnd]
        softv[ibnd,n] *= c_soft_non[value_locate(rx_fine,rx_samp[ibnd]),ibnd]
        hardv[iunb,n] *= c_hard_non[-1,iunb]
        softv[iunb,n] *= c_soft_non[-1,iunb]
        ;; simulate contamination
        cfullv[*,n] = fullv[*,n]
        chardv[*,n] = hardv[*,n]
        csoftv[*,n] = softv[*,n]
        irand = randomi(nnon*0.10,nnon,/nodup)
        cfullv[irand,n] = 0.
        chardv[irand,n] = 0.
        csoftv[irand,n] = 0.                        
    endfor
    ;; find the mean flux over all sources
    full = mean(fullv,dim=1,/nan)
    hard = mean(hardv,dim=1,/nan)
    soft = mean(softv,dim=1,/nan)
    e_full = medabsdev(fullv,dim=1)
    e_hard = medabsdev(hardv,dim=1)
    e_soft = medabsdev(softv,dim=1)
    cfull = mean(cfullv,dim=1,/nan)
    chard = mean(chardv,dim=1,/nan)
    csoft = mean(csoftv,dim=1,/nan)
    e_cfull = medabsdev(cfullv,dim=1)
    e_chard = medabsdev(chardv,dim=1)
    e_csoft = medabsdev(csoftv,dim=1)
    ;; and the mode of the 
    if (niter gt 1) then begin
        full = median(mode(full,kde=kde_bandwidth(full)),/even)
        hard = median(mode(hard,kde=kde_bandwidth(hard)),/even)
        soft = median(mode(soft,kde=kde_bandwidth(soft)),/even)
        e_full = median(mode(e_full,kde=kde_bandwidth(e_full)),/even)
        e_hard = median(mode(e_hard,kde=kde_bandwidth(e_hard)),/even)
        e_soft = median(mode(e_soft,kde=kde_bandwidth(e_soft)),/even)
        cfull = median(mode(cfull,kde=kde_bandwidth(cfull)),/even)
        chard = median(mode(chard,kde=kde_bandwidth(chard)),/even)
        csoft = median(mode(csoft,kde=kde_bandwidth(csoft)),/even)
        e_cfull = median(mode(e_cfull,kde=kde_bandwidth(e_cfull)),/even)
        e_chard = median(mode(e_chard,kde=kde_bandwidth(e_chard)),/even)
        e_csoft = median(mode(e_csoft,kde=kde_bandwidth(e_csoft)),/even)
    endif
    fx_full[i] = full
    fx_hard[i] = hard
    fx_soft[i] = soft
    e_fx_full[i] = e_full
    e_fx_hard[i] = e_hard
    e_fx_soft[i] = e_soft
    cfx_full[i] = cfull
    cfx_hard[i] = chard
    cfx_soft[i] = csoft
    e_cfx_full[i] = e_cfull
    e_cfx_hard[i] = e_chard
    e_cfx_soft[i] = e_csoft
endfor
;; reconstruct original array if greater than 2D
dim = strjoin(strtrim(sz[2:sz[0]],2),',')
for i = 0,n_elements(lin_vars)-1 do re = execute(lin_vars[i]+' = reform('+lin_vars[i]+','+dim+')')
for i = 0,n_elements(clin_vars)-1 do re = execute(clin_vars[i]+' = reform('+clin_vars[i]+','+dim+')')

;; save to structure
eng = ['FULL','HARD','SOFT']
lin = [eng,'E_'+eng]
ceng = ['CFULL','CHARD','CSOFT']
clin = [ceng,'E_'+ceng]

str = strjoin([lin+':'+lin_vars,clin+':'+clin_vars],',')
re = execute('fxest = soa2aos({'+str+'})')

return, fxest


END

