PRO estimate_fx, CHA = cha


common _data
common _nhobs
common _rxnh
common _group
common _ctfest
common _split
common _rxmod

;; Follow procedure of Carroll+20 (XSTACK_OUTPUT, STACK_NONDET_FLUX, MC_NONDET_DIST)

;; use only Chandra sources to match Chandra X-ray stacking
if keyword_set(cha) then iicha = xnon eq 'CHA'
                         iicha = xnon ne ''

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
    else: message, message, 'NO INPUT MODE SET FOR FXEST: WAC/SEC/ALL'
endcase

;; prep observed frame hard and soft conversion arrays
iz = value_locate(zv,z_non)
c_hard_non = c_hard[*,iz]
c_soft_non = c_soft[*,iz]

;; run for KS test

;; capture mean and median of the modeled log FX distribution 
;niter = 10000
niter = n_elements(rx_mod2_ksv[0,*])
logfx_full_ksv = dblarr(nnon,niter)
logfx_hard_ksv = dblarr(nnon,niter)
logfx_soft_ksv = dblarr(nnon,niter)
;; in Carroll+20, we take the model RL and subtract from it the observations
;; this doesn't work here as the observations are greater than the simulated observed
;; for now, subtracting the model detected 
;inon = where(iimod2_ks[*,iks2] eq 0,nonct)
;if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
;rx_modn_ks = (rx_mod2_ks[*,iks2])[inon]

;; sample from model RL for each non-detected observation (NNON)
for j = 0,niter-1 do begin    
    ;; test new method
    inon = where(iimod2_ksv[*,j] eq 0,nonct)
    if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
    rx_modn_ks = (rx_mod2_ksv[*,j])[inon]    
    ;;

    rx_samp_ks = mc_samp(rx_modn_ks,nnon)
    loglx_full_ks = rx_samp_ks + loglxir_non
    logfx_full_ks = loglx_full_ks - alog10(4.*!const.pi*dl2_non)
    logfx_full_ksv[*,j] = logfx_full_ks
    ;; convert 2-10keV to soft and hard fluxes
    logfx_hard_ks = logfx_full_ks
    logfx_soft_ks = logfx_full_ks
    ibnd = where(rx_samp_ks ge min(rx),nbnd,complement=iunb,ncomplement=nunb)
    for i = 0,nbnd-1 do logfx_hard_ks[ibnd[i]] += interpol(c_hard_non[*,ibnd[i]],rx,rx_samp_ks[ibnd[i]])
    for i = 0,nbnd-1 do logfx_soft_ks[ibnd[i]] += interpol(c_soft_non[*,ibnd[i]],rx,rx_samp_ks[ibnd[i]])
    logfx_hard_ks[iunb] += c_hard_non[-1,iunb]
    logfx_soft_ks[iunb] += c_soft_non[-1,iunb]
    logfx_hard_ksv[*,j] = logfx_hard_ks
    logfx_soft_ksv[*,j] = logfx_soft_ks
endfor
logfx_full_ks = median(logfx_full_ksv,dim=2)
logfx_hard_ks = median(logfx_hard_ksv,dim=2)
logfx_soft_ks = median(logfx_soft_ksv,dim=2)
logfx_full_ksx = median(logfx_full_ks)
logfx_hard_ksx = median(logfx_hard_ks)
logfx_soft_ksx = median(logfx_soft_ks)
e_logfx_full_ksx = stddev(logfx_full_ks)
e_logfx_hard_ksx = stddev(logfx_hard_ks)
e_logfx_soft_ksx = stddev(logfx_soft_ks)
fx_full_ksx = 10.^logfx_full_ksx
fx_hard_ksx = 10.^logfx_hard_ksx
fx_soft_ksx = 10.^logfx_soft_ksx
e_fx_full_ksx = e_logfx_full_ksx * alog(10.) * fx_full_ksx
e_fx_hard_ksx = e_logfx_hard_ksx * alog(10.) * fx_hard_ksx
e_fx_soft_ksx = e_logfx_soft_ksx * alog(10.) * fx_soft_ksx

sav_vars = ['LOGFX_FULL_KSV','LOGFX_HARD_KSV','LOGFX_SOFT_KSV', $
            'LOGFX_FULL_KS','LOGFX_HARD_KS','LOGFX_SOFT_KS', $
            'LOGFX_FULL_KSX','LOGFX_HARD_KSX','LOGFX_SOFT_KSX', $
            'E_LOGFX_FULL_KSX','E_LOGFX_HARD_KSX','E_LOGFX_SOFT_KSX', $
            'FX_FULL_KSX','FX_HARD_KSX','FX_SOFT_KSX', $
            'E_FX_FULL_KSX','E_FX_HARD_KSX','E_FX_SOFT_KSX']
sav_inds = []


;; run for AD test

;; capture mean and median of the modeled log FX distribution 
;niter = 10000
niter = n_elements(rx_mod2_adv[0,*])
logfx_full_adv = dblarr(nnon,niter)
logfx_hard_adv = dblarr(nnon,niter)
logfx_soft_adv = dblarr(nnon,niter)

;; in Carroll+20, we take the model RL and subtract from it the observations
;; this doesn't work here as the observations are greater than the simulated observed
;; for now, subtracting the model detected 
;inon = where(iimod2_ad[*,iad2] eq 0,nonct)
;if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
;rx_modn_ad = (rx_mod2_ad[*,iad2])[inon]

;; sample from model RL for each non-detected observation (NNON)
for j = 0,niter-1 do begin  
    ;; test new method
    inon = where(iimod2_adv[*,j] eq 0,nonct)
    if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
    rx_modn_ad = (rx_mod2_adv[*,j])[inon]
    ;;
  
    rx_samp_ad = mc_samp(rx_modn_ad,nnon)
    loglx_full_ad = rx_samp_ad + loglxir_non
    logfx_full_ad = loglx_full_ad - alog10(4.*!const.pi*dl2_non)
    logfx_full_adv[*,j] = logfx_full_ad
    ;; convert 2-10keV to soft and hard fluxes
    logfx_hard_ad = logfx_full_ad
    logfx_soft_ad = logfx_full_ad
    ibnd = where(rx_samp_ad ge min(rx),nbnd,complement=iunb,ncomplement=nunb)
    for i = 0,nbnd-1 do logfx_hard_ad[ibnd[i]] += interpol(c_hard_non[*,ibnd[i]],rx,rx_samp_ad[ibnd[i]])
    for i = 0,nbnd-1 do logfx_soft_ad[ibnd[i]] += interpol(c_soft_non[*,ibnd[i]],rx,rx_samp_ad[ibnd[i]])
    logfx_hard_ad[iunb] += c_hard_non[-1,iunb]
    logfx_soft_ad[iunb] += c_soft_non[-1,iunb]
    logfx_hard_adv[*,j] = logfx_hard_ad
    logfx_soft_adv[*,j] = logfx_soft_ad
endfor
logfx_full_ad = median(logfx_full_adv,dim=2)
logfx_hard_ad = median(logfx_hard_adv,dim=2)
logfx_soft_ad = median(logfx_soft_adv,dim=2)
logfx_full_adx = median(logfx_full_ad)
logfx_hard_adx = median(logfx_hard_ad)
logfx_soft_adx = median(logfx_soft_ad)
e_logfx_full_adx = stddev(logfx_full_ad)
e_logfx_hard_adx = stddev(logfx_hard_ad)
e_logfx_soft_adx = stddev(logfx_soft_ad)
fx_full_adx = 10.^logfx_full_adx
fx_hard_adx = 10.^logfx_hard_adx
fx_soft_adx = 10.^logfx_soft_adx
e_fx_full_adx = e_logfx_full_adx * alog(10.) * fx_full_adx
e_fx_hard_adx = e_logfx_hard_adx * alog(10.) * fx_hard_adx
e_fx_soft_adx = e_logfx_soft_adx * alog(10.) * fx_soft_adx

sav_vars = [sav_vars, $
            'LOGFX_FULL_ADV','LOGFX_HARD_ADV','LOGFX_SOFT_ADV', $
            'LOGFX_FULL_AD','LOGFX_HARD_AD','LOGFX_SOFT_AD', $
            'LOGFX_FULL_ADX','LOGFX_HARD_ADX','LOGFX_SOFT_ADX', $
            'E_LOGFX_FULL_ADX','E_LOGFX_HARD_ADX','E_LOGFX_SOFT_ADX', $
            'FX_FULL_ADX','FX_HARD_ADX','FX_SOFT_ADX', $
            'E_FX_FULL_ADX','E_FX_HARD_ADX','E_FX_SOFT_ADX']
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',file="fx_estimate.sav"')


END

