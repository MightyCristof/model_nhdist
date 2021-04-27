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
if keyword_set(iterate) then niter = 1000 else $
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
log_vars = 'LOG'+lin_vars
lin_vars = [lin_vars,'E_'+lin_vars]
for i = 0,n_elements(lin_vars)-1 do re = execute(lin_vars[i]+' = dblarr(nmod)')

;; sample from model RL for each non-detected observation (NNON)
for i = 0,nmod-1 do begin  
    ;; "non-detected" model sources
    inon = where(ii_modd[*,i] eq 0,nonct)
    if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
    rx_modn = (rx_modd[*,i])[inon]
    fullv = dblarr(nnon,niter)
    hardv = dblarr(nnon,niter)
    softv = dblarr(nnon,niter)    
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
    endfor
    ;; find the mean flux over all sources
    full = mean(fullv,dim=1)
    hard = mean(hardv,dim=1)
    soft = mean(softv,dim=1)
    e_full = medabsdev(fullv,dim=1)
    e_hard = medabsdev(hardv,dim=1)
    e_soft = medabsdev(softv,dim=1)
    ;; and the mode of the 
    if (niter gt 1) then begin
        full = mode(full)
        hard = mode(hard)
        soft = mode(soft)
        e_full = mode(e_full)
        e_hard = mode(e_hard)
        e_soft = mode(e_soft)
    endif
    fx_full[i] = full
    fx_hard[i] = hard
    fx_soft[i] = soft
    e_fx_full[i] = e_full
    e_fx_hard[i] = e_hard
    e_fx_soft[i] = e_soft
endfor
;; reconstruct original array if greater than 2D
dim = strjoin(strtrim(sz[2:sz[0]],2),',')
for i = 0,n_elements(lin_vars)-1 do re = execute(lin_vars[i]+' = reform('+lin_vars[i]+','+dim+')')
;; log space
for i = 0,n_elements(log_vars)-1 do re = execute(log_vars[i]+' = alog10('+lin_vars[i]+')')
log_vars = [log_vars,'E_'+log_vars]
for i = i,n_elements(log_vars)-1 do re = execute(log_vars[i]+' = '+lin_vars[i]+'/(alog(10.)*'+lin_vars[i-3]+')')

;; save to structure
eng = ['FULL','HARD','SOFT']
lin = [eng,'E_'+eng]
log = ['LOG'+eng,'E_LOG'+eng]

str = strjoin([lin+':'+lin_vars,log+':'+log_vars],',')
re = execute('fxest = soa2aos({'+str+'})')

return, fxest


END

