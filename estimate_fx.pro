PRO estimate_fx


common _data
common _nhobs
common _ctfest
common _split
common _rxmod

;; Follow procedure of Carroll+20 (XSTACK_OUTPUT, STACK_NONDET_FLUX, MC_NONDET_DIST)

;; LIR and luminosity distance for non-detected sample
ixn = where(iixn,nnon)
loglir_non = loglir[ixn]
dlumsq = dl2[ixn]

;; run for KS test

;; capture mean and median of the modeled log FX distribution 
niter = 10000
logfx_avg_ksv = dblarr(niter)
logfx_med_ksv = dblarr(niter)

;; in Carroll+20, we take the model RL and subtract from it the observations
;; this doesn't work here as the observations are greater than the simulated observed
;; for now, subtracting the model detected 
rx_modn_ks = (rx_mod2_ks.(iks2))[where(iimod2_ks.(iks2) eq 0,nmodn_ks)]

;; sample from model RL for each non-detected observation (NNON)
print, 'KS: ESTIMATE_FX - 0%'
for j = 0,niter-1 do begin    
    rx_samp_ks = mc_rldist(rx_modn_ks,nnon)
    loglx_non_ks = rx_samp_ks + loglir_non
    logfx_non_ks = loglx_non_ks - alog10(4.*!const.pi*dlumsq)
    logfx_avg_ksv[j] = mean(logfx_non_ks)
    logfx_med_ksv[j] = median(logfx_non_ks)
    if (j gt 0 and j mod (niter/2.) eq 0) then print, 'KS: ESTIMATE_FX - 50% COMPLETE'
endfor
print, 'KS: ESTIMATE_FX - 100% COMPLETE'


sav_vars = ['LOGFX_AVG_KSV','LOGFX_MED_KSV']
sav_inds = []


;; run for AD test

;; capture mean and median of the modeled log FX distribution 
niter = 10000
logfx_avg_adv = dblarr(niter)
logfx_med_adv = dblarr(niter)

;; in Carroll+20, we take the model RL and subtract from it the observations
;; this doesn't work here as the observations are greater than the simulated observed
;; for now, subtracting the model detected 
rx_modn_ad = (rx_mod2_ad.(iad2))[where(iimod2_ad.(iad2) eq 0,nmodn_ad)]

;; sample from model RL for each non-detected observation (NNON)
print, 'AD: ESTIMATE_FX - 0%'
for j = 0,niter-1 do begin    
    rx_samp_ad = mc_rldist(rx_modn_ad,nnon)
    loglx_non_ad = rx_samp_ad + loglir_non
    logfx_non_ad = loglx_non_ad - alog10(4.*!const.pi*dlumsq)
    logfx_avg_adv[j] = mean(logfx_non_ad)
    logfx_med_adv[j] = median(logfx_non_ad)
    if (j gt 0 and j mod (niter/2.) eq 0) then print, 'AD: ESTIMATE_FX - 50% COMPLETE'
endfor
print, 'AD: ESTIMATE_FX - 100% COMPLETE'

sav_vars = [sav_vars,'LOGFX_AVG_ADV','LOGFX_MED_ADV']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',file="fx_estimate.sav"')


END

