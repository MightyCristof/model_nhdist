FUNCTION estimate_fx, rx_model, $
                      iidet, $
                      ITERATE = iterate, $
                      DETECTIONS = detections


common _data
common _nhdist
common _nhobs
common _rxnh
common _group

                         
;; iterate over each object, or run single case
if keyword_set(iterate) then niter = 100 else $
                             niter = 1

;; use only Chandra sources to match Chandra X-ray stacking
;; set detections vs non-detections
if keyword_set(detections) then begin
    iicha = xdet eq 'CHA'
    iiflg = iidet eq 1
endif else begin
    iicha = xnon eq 'CHA'
    iiflg = iidet eq 0
endelse

;; data for non-detected sample        
case group of 
    'WAC': if keyword_set(detections) then iigrp = iiwd else iigrp = iiwn
    'SEC': if keyword_set(detections) then iigrp = iisd else iigrp = iisn
    'ALL': if keyword_set(detections) then iigrp = iiad else iigrp = iian
    'WAC_HIX': if keyword_set(detections) then iigrp = iiwd and iixh else iigrp = iiwn and iixh
    'WAC_LOX': if keyword_set(detections) then iigrp = iiwd and iixl else iigrp = iiwn and iixl
    'OFFSET': if keyword_set(detections) then iigrp = iiwd else iigrp = iiwn
    else: message, message, 'NO INPUT MODE SET FOR FXEST: WAC/SEC/ALL'
endcase
;; set variables for chosen group
igrp = where(iigrp and iicha,ngrp)
loglxir_grp = loglxir[igrp]
z_grp = z[igrp]
dl2_grp = dl2[igrp]

;; prep hard and soft conversion arrays
iz = value_locate(zv,z_grp)
;c_hard_grp = c_hard[*,iz]
;c_soft_grp = c_soft[*,iz]
c_hard_grp = c_hard_fine[*,iz]
c_soft_grp = c_soft_fine[*,iz]

;; reform to 2D array if needed
sz = size(rx_model)
nmod = product(sz[2:sz[0]])
rx_modd = reform(rx_model,sz[1],nmod)
ii_flag = reform(iiflg,sz[1],nmod)
lin_vars = 'FX_'+['FULL','HARD','SOFT']
clin_vars = 'C'+lin_vars
lin_vars = [lin_vars,'E_'+lin_vars]
clin_vars = [clin_vars,'E_'+clin_vars]
for i = 0,n_elements(lin_vars)-1 do re = execute(lin_vars[i]+' = dblarr(nmod)')
for i = 0,n_elements(clin_vars)-1 do re = execute(clin_vars[i]+' = dblarr(nmod)')

;; sample from model RL for each group observation (ngrp)
for i = 0,nmod-1 do begin  
    ;; "flagged" model sources
    iflg = where(ii_flag[*,i],flgct)
    if (flgct eq 0) then message, 'NO FLAGGED SOURCES IN MODEL.'
    rx_modn = (rx_modd[*,i])[iflg]
    fullv = dblarr(ngrp,niter)
    hardv = dblarr(ngrp,niter)
    softv = dblarr(ngrp,niter)    
    cfullv = dblarr(ngrp,niter)
    chardv = dblarr(ngrp,niter)
    csoftv = dblarr(ngrp,niter)    
    for n = 0,niter-1 do begin
        ;; random draw from model to match number of sources
        rx_samp = mc_samp(rx_modn,ngrp)
        ;; estimate flux
        lx_grp = 10.^(rx_samp+loglxir_grp)
        fullv[*,n] = lx_grp/(4.*!const.pi*dl2_grp)
        hardv[*,n] = fullv[*,n]
        softv[*,n] = fullv[*,n]
        ;; bounded indices
        ibnd = where(rx_samp ge min(rx),nbnd,complement=iunb,ncomplement=nunb)
        hardv[ibnd,n] *= c_hard_grp[value_locate(rx_fine,rx_samp[ibnd]),ibnd]
        softv[ibnd,n] *= c_soft_grp[value_locate(rx_fine,rx_samp[ibnd]),ibnd]
        hardv[iunb,n] *= c_hard_grp[-1,iunb]
        softv[iunb,n] *= c_soft_grp[-1,iunb]
        ;; simulate contamination
        cfullv[*,n] = fullv[*,n]
        chardv[*,n] = hardv[*,n]
        csoftv[*,n] = softv[*,n]
        irand = randomi(ngrp*0.10,ngrp,/nodup)  ;; reliability
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

