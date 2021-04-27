PRO estimate_fx, rx_model, $
                 iidet, $
                 CHA = cha


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
;common _fixed
;common _free
;common _model

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

;; capture mean and median of the modeled log FX distribution 
;niter = 10000
niter = n_elements(rx_model[0,*])
logfx_fullv = dblarr(nnon,niter)
logfx_hardv = dblarr(nnon,niter)
logfx_softv = dblarr(nnon,niter)

;; in Carroll+20, we take the model RL and subtract from it the observations
;; this doesn't work here as the observations are greater than the simulated observed
;; for now, subtracting the model detected 
;inon = where(iimod2_ad[*,iad2] eq 0,nonct)
;if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
;rx_modn_ad = (rx_mod2_ad[*,iad2])[inon]

;; sample from model RL for each non-detected observation (NNON)
for n = 0,niter-1 do begin  
    ;; test new method
    inon = where(iidet[*,n] eq 0,nonct)
    if (nonct eq 0) then message, 'NO NON-DETECTIONS IN MODEL.'
    rx_modn = (rx_model[*,n])[inon]
    ;;
  
    rx_samp = mc_samp(rx_modn,nnon)
    loglx_full = rx_samp + loglxir_non
    
    logfx_full = loglx_full - alog10(4.*!const.pi*dl2_non)
    logfx_fullv[*,n] = logfx_full
    ;; convert 2-10keV to soft and hard fluxes
    logfx_hard = logfx_full
    logfx_soft = logfx_full
    ibnd = where(rx_samp ge min(rx),nbnd,complement=iunb,ncomplement=nunb)
    for i = 0,nbnd-1 do logfx_hard[ibnd[i]] += interpol(c_hard_non[*,ibnd[i]],rx,rx_samp[ibnd[i]])
    for i = 0,nbnd-1 do logfx_soft[ibnd[i]] += interpol(c_soft_non[*,ibnd[i]],rx,rx_samp[ibnd[i]])
    logfx_hard[iunb] += c_hard_non[-1,iunb]
    logfx_soft[iunb] += c_soft_non[-1,iunb]
    logfx_hardv[*,n] = logfx_hard
    logfx_softv[*,n] = logfx_soft
endfor
logfx_full = median(logfx_fullv,dim=2)
logfx_hard = median(logfx_hardv,dim=2)
logfx_soft = median(logfx_softv,dim=2)
logfx_fullx = median(logfx_full)
logfx_hardx = median(logfx_hard)
logfx_softx = median(logfx_soft)
e_logfx_fullx = stddev(logfx_full)
e_logfx_hardx = stddev(logfx_hard)
e_logfx_softx = stddev(logfx_soft)
fx_fullx = 10.^logfx_fullx
fx_hardx = 10.^logfx_hardx
fx_softx = 10.^logfx_softx
e_fx_fullx = e_logfx_fullx * alog(10.) * fx_fullx
e_fx_hardx = e_logfx_hardx * alog(10.) * fx_hardx
e_fx_softx = e_logfx_softx * alog(10.) * fx_softx

sav_vars = ['LOGFX_FULLV','LOGFX_HARDV','LOGFX_SOFTV', $
            'LOGFX_FULL','LOGFX_HARD','LOGFX_SOFT', $
            'LOGFX_FULLX','LOGFX_HARDX','LOGFX_SOFTX', $
            'E_LOGFX_FULLX','E_LOGFX_HARDX','E_LOGFX_SOFTX', $
            'FX_FULLX','FX_HARDX','FX_SOFTX', $
            'E_FX_FULLX','E_FX_HARDX','E_FX_SOFTX']
sav_inds = []
stop
sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',file="fx_estimate.sav"')


END

