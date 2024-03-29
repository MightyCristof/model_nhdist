PRO model_rxdist, POSTMOD = postmod


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
common _fixed
common _free
common _split

;; use observed NH dist with added unobscured sources after modeling without it to find
if keyword_set(postmod) then nh_obs = nh_lan_cor else $
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
fct = mode(fctv_fine,kde=kde_bandwidth(fctv_fine));mean([mode(fctv,bin=scott(fctv)),mode(fctv,bin=freedman(fctv)),mode(fctv,kde=kde_bandwidth(fctv))])
e_fct = stddev(fctv_fine)
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
;; create arrays
nh_modv = dblarr(nsrc,niter)
rx_modv = dblarr(nsrc,niter)
iimodv = bytarr(nsrc,niter)
a2 = dblarr(niter)
p_a2 = dblarr(niter)
;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
for n = 0,niter-1 do begin
    nh_modv[*,n] = nh_resamp[randomi(nsrc,nresamp)]
    rx_modv[*,n] = rx2nh(nh_modv[*,n],/rx_out,scat=rx_scat)
    iimodv[*,n] = rx_modv[*,n] gt rxl
    idet = where(iimodv[*,n] eq 1,detct)
    rxd_scat = randomn(seed,ndet)*e_rxd
    if (detct ge 5) then begin
        a2[n] = ad_test(rxd+rxd_scat,rx_modv[idet,n],permute=(test eq 'JOINT'),prob=p,cvm=cvm)
        p_a2[n] = p
    endif else if (detct gt 0) then begin
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

sav_vars = ['NSR','NCN','FCN', $
            'NCT','FCT','E_FCT', $
            'NH_MODV','RX_MODV','A2','P_A2', $
            'NH_MOD','RX_MOD','FX_NON','FX_DET']
sav_inds = ['IIMODV','POSTMOD']


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RUN FOR SPLIT CTF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fraction and source numbers
fct_ = mode(fctv1_fine,kde=kde_bandwidth(fctv1_fine));mean([mode(fctv1_fine,bin=scott(fctv1_fine)),mode(fctv1_fine,bin=freedman(fctv1_fine)),mode(fctv1_fine,kde=kde_bandwidth(fctv1_fine))])
e_fct_ = stddev(fctv1_fine)
fcn_ = 1.-fct_
ncn_ = nthin
nsr_ = round(ncn_/fcn_)
nct_ = round(nsr_*fct_)
f24_ = mean(f24v2)*fct_;mode(f24v2,kde=kde_bandwidth(f24v2))*fct_
n24_ = round(nsr_*f24_)
f25_ = mean(f25v2)*fct_;mode(f25v2,kde=kde_bandwidth(f25v2))*fct_
n25_ = round(nsr_*f25_)
;n25_ = nct_-n24_
;; NH_RESAMP: structure of increased CT sources
;; NH_MOD: draw N sources from NH_RESAMP
;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where RX model > RX limits from Carroll+21
;; AD: Anderson-Darling test statistic    
nh_resamp = [nh_samp[ithin],24.+randomu(seed,n24_),25.+randomu(seed,n25_)]
nresamp = n_elements(nh_resamp)
;; create arrays
nh_modv_ = dblarr(nsrc,niter)
rx_modv_ = dblarr(nsrc,niter)
iimodv_ = bytarr(nsrc,niter)
a2_ = dblarr(niter)
p_a2_ = dblarr(niter)
;; fill arrays
for n = 0,niter-1 do begin
    nh_modv_[*,n] = nh_resamp[randomi(nsrc,nresamp)]
    rx_modv_[*,n] = rx2nh(nh_modv_[*,n],/rx_out,scat=rx_scat)
    iimodv_[*,n] = rx_modv_[*,n] gt rxl
    idet = where(iimodv_[*,n] eq 1,detct)
    rxd_scat = randomn(seed,ndet)*e_rxd
    if (detct ge 5) then begin
        a2_[n] = ad_test(rxd+rxd_scat,rx_modv[idet,n],permute=(test eq 'JOINT'),prob=p,cvm=cvm)
        p_a2_[n] = p
    endif else if (detct gt 0) then begin
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

sav_vars = [sav_vars,'NSR_','NCN_','FCN_', $
                     'NCT_','FCT_','E_FCT_', $
                     'N24_','F24_','N25_','F25_', $
                     'NH_MODV_','RX_MODV_','A2_','P_A2_', $
                     'NH_MOD_','RX_MOD_','FX_NON_','FX_DET_']
sav_inds = [sav_inds,'IIMODV_']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="rx_model.sav"')


END









