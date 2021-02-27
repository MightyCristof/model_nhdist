PRO estimate_fx2, CHA = cha


common _data
common _nhobs
common _rxnh
common _ctfest
common _split
common _rxmod

;; Follow procedure of Carroll+20 (XSTACK_OUTPUT, STACK_NONDET_FLUX, MC_NONDET_DIST)

;; use only Chandra sources to match Chandra X-ray stacking
if keyword_set(cha) then iwn = where(iiwn and xnon eq 'CHA',nnon) else $
                         iwn = where(iiwn,nnon)
;; data for non-detected sample
loglxir_non = loglxir[iwn]
z_non = z[iwn]
dl2_non = dl2[iwn]
;; prep observed frame hard and soft conversion arrays
iz = value_locate(zv,z_non)
c_hard_non = c_hard[*,iz]
c_soft_non = c_soft[*,iz]

;; run for KS test

;; capture mean and median of the modeled log FX distribution 
niter = 10000
logfx_full_ksv = dblarr(nnon,niter)
logfx_hard_ksv = dblarr(nnon,niter)
logfx_soft_ksv = dblarr(nnon,niter)
;; in Carroll+20, we take the model RL and subtract from it the observations
;; this doesn't work here as the observations are greater than the simulated observed
;; for now, subtracting the model detected 

;inon = where(iimod2_ks[*,iks2] eq 0,nonct)
;if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
;rx_modn_ks = (rx_mod2_ks[*,iks2])[inon]
rx_modn_ks = hist2d_avg(rx_mod2_ksv[where(iimod2_ksv eq 0)],/hist)

;; sample from model RL for each non-detected observation (NNON)
for j = 0,niter-1 do begin    
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
stop
sav_vars = ['LOGFX_FULL_KSV','LOGFX_HARD_KSV','LOGFX_SOFT_KSV']
sav_inds = []


;; run for AD test

;; capture mean and median of the modeled log FX distribution 
niter = 10000
logfx_full_adv = dblarr(nnon,niter)
logfx_hard_adv = dblarr(nnon,niter)
logfx_soft_adv = dblarr(nnon,niter)

;; in Carroll+20, we take the model RL and subtract from it the observations
;; this doesn't work here as the observations are greater than the simulated observed
;; for now, subtracting the model detected 
inon = where(iimod2_ad[*,iad2] eq 0,nonct)
if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
rx_modn_ad = (rx_mod2_ad[*,iad2])[inon]

;; sample from model RL for each non-detected observation (NNON)
for j = 0,niter-1 do begin    
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

sav_vars = [sav_vars,'LOGFX_FULL_ADV','LOGFX_HARD_ADV','LOGFX_SOFT_ADV']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',file="fx_estimate.sav"')


END

