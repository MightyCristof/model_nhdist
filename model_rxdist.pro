PRO model_rxdist, POSTMOD = postmod


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
;common _uniform
;common _variable


;; use observed NH dist with added unobscured sources after modeling without it to find
if keyword_set(postmod) then nh_obs = nh_ric_int else $
                             postmod = 0
cvm = 0

;; number of iterations for each test
niter = 10000

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RUN FOR FIXED CTF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; resample observed NH distribution to increase data density
nsamp = nsrc*100.
nh_samp = nh_mc(nh_obs,nsamp)
;; number of Compton-thin sources
ithin = where(nh_samp lt 24.,nthin);,complement=ithick,ncomplement=nthick)

;; fraction and source numbers
;fct = mode(fctv,kde=kde_bandwidth(fctv))
;fct = 0.562
;e_fct = stddev(fctv)

;; mcmc results
fct = 0.555857
e_fct = 0.0345998

fcn = 1.-fct
ncn = nthin
nsr = round(ncn/fcn)
nct = round(nsr*fct)
;; NH_RESAMP: structure of increased CT sources
;; NH_MOD: draw N sources from NH_RESAMP
;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where RX model > RX limits from Carroll+21
;; AD: Anderson-Darling test statistic    
nh_resamp = [nh_samp[ithin],24.+2.*randomu(seed,nct)]
nresamp = n_elements(nh_resamp)
;; modeling variables
nh_modv = dblarr(nsrc,niter)
rx_modv = dblarr(nsrc,niter)
iimodv = bytarr(nsrc,niter)
a2 = dblarr(niter)
p_a2 = dblarr(niter)
;; data variables
rx_detv = dblarr(ndet,niter)
;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
for n = 0,niter-1 do begin
    nh_modv[*,n] = nh_resamp[randomi(nsrc,nresamp)]
    rx_modv[*,n] = rx2nh(nh_modv[*,n],/rx_out,scat=rx_scat,/mcmc)
    iimodv[*,n] = rx_modv[*,n] gt rxl
    idet = where(iimodv[*,n] eq 1,moddet)
    rx_detv[*,n] = rxd+randomn(seed,ndet)*0.23;rx_scat
    ;rx_detv[*,n] = rxd+randomn(seed,ndet)*e_rxd
    if (moddet ge 5) then begin
        a2[n] = ad_test(rx_detv[*,n],rx_modv[idet,n],permute=0,prob=p,cvm=cvm)
        p_a2[n] = p
    endif else if (moddet gt 0) then begin
        a2[n] = -1.
        p_a2[n] = -1.
    endif else message, 'NO MODELED DETECTIONS.'
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FIXED FCT MODELING'
    endif else if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FIXED FCT MODELING'
print, '=============================================='
nh_mod = hist2d_avg(nh_modv,1.d,iidet=iimodv,normalize='where(xh lt 24.,/null)')
rx_mod = hist2d_avg(rx_modv,0.2d,iidet=iimodv,/normalize)
fx_non = estimate_fx(rx_modv,iimodv)
fx_det = estimate_fx(rx_modv,iimodv,/detections)
rx_det = hist2d_avg(rx_detv,0.2d,/normalize)

sav_vars = ['NSR','NCN','FCN', $
            'NCT','FCT','E_FCT', $
            'NH_MODV','RX_MODV','A2','P_A2', $
            'NH_MOD','RX_MOD','FX_NON','FX_DET', $
            'RX_DET']
sav_inds = ['IIMODV','POSTMOD']


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RUN FOR SPLIT CTF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fraction and source numbers
;fct_ = mode(fctv,kde=kde_bandwidth(fctv))
;e_fct_ = stddev(fctv)
;fcn_ = 1.-fct_
;ncn_ = nthin
;nsr_ = round(ncn_/fcn_)
;nct_ = round(nsr_*fct_)
;f24_ = mode(f24v_,kde=kde_bandwidth(f24v_))*fct_
;;f24_ = mean(f24v2)*fct_;mode(f24v2,kde=kde_bandwidth(f24v2))*fct_
;n24_ = round(nsr_*f24_)
;f25_ = mode(f25v_,kde=kde_bandwidth(f25v_))*fct_
;;f25_ = mean(f25v2)*fct_;mode(f25v2,kde=kde_bandwidth(f25v2))*fct_
;n25_ = round(nsr_*f25_)
;;n25_ = nct_-n24_

fct_ = fct
e_fct_ = e_fct
fcn_ = 1.-fct_
ncn_ = nthin
nsr_ = round(ncn_/fcn_)
nct_ = round(nsr_*fct_)
n24_ = round(nct_ * 0.45)
n25_ = round(nct_ * 0.55)

;; NH_RESAMP: structure of increased CT sources
;; NH_MOD: draw N sources from NH_RESAMP
;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where RX model > RX limits from Carroll+21
;; AD: Anderson-Darling test statistic    
nh_resamp = [nh_samp[ithin],24.+randomu(seed,n24_),25.+randomu(seed,n25_)]
nresamp = n_elements(nh_resamp)
;; modeling variables
nh_modv_ = dblarr(nsrc,niter)
rx_modv_ = dblarr(nsrc,niter)
iimodv_ = bytarr(nsrc,niter)
a2_ = dblarr(niter)
p_a2_ = dblarr(niter)
;; data variables
rx_detv_ = dblarr(ndet,niter)
;; fill arrays
for n = 0,niter-1 do begin
    nh_modv_[*,n] = nh_resamp[randomi(nsrc,nresamp)]
    rx_modv_[*,n] = rx2nh(nh_modv_[*,n],/rx_out,scat=rx_scat,/mcmc)
    iimodv_[*,n] = rx_modv_[*,n] gt rxl
    idet = where(iimodv_[*,n] eq 1,moddet)
    rx_detv_[*,n] = rxd+randomn(seed,ndet)*0.23;rx_scat
    ;rx_detv_[*,n] = rxd+randomn(seed,ndet)*e_rxd
    if (moddet ge 5) then begin
        a2_[n] = ad_test(rx_detv_[*,n],rx_modv_[idet,n],permute=0,prob=p,cvm=cvm)
        p_a2_[n] = p
    endif else if (moddet gt 0) then begin
        a2_[n] = -1.
        p_a2_[n] = -1.
    endif else message, 'NO MODELED DETECTIONS.'
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FREE FCT MODELING'
    endif else if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FREE FCT MODELING'
print, '=============================================='

nh_mod_ = hist2d_avg(nh_modv_,1.d,iidet=iimodv_,normalize='where(xh lt 24.,/null)')
rx_mod_ = hist2d_avg(rx_modv_,0.2d,iidet=iimodv_,/normalize)
fx_non_ = estimate_fx(rx_modv_,iimodv_)
fx_det_ = estimate_fx(rx_modv_,iimodv_,/detections)
rx_det_ = hist2d_avg(rx_detv_,0.2d,/normalize)

sav_vars = [sav_vars,'NSR_','NCN_','FCN_', $
                     'NCT_','FCT_','E_FCT_', $
                     'N24_','F24_','N25_','F25_', $
                     'NH_MODV_','RX_MODV_','A2_','P_A2_', $
                     'NH_MOD_','RX_MOD_','FX_NON_','FX_DET_', $
                     'RX_DET_']
sav_inds = [sav_inds,'IIMODV_']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="rx_model.sav"')


END









