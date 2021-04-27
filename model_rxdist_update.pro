PRO model_rxdist_update, POSTMOD = postmod


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
common _fixed
common _free


;; use observed NH dist with added unobscured sources after modeling without it to find
if keyword_set(postmod) then nh_obs = nh_lan_cor

;; resample observed NH distribution to increase data density
nsamp = nsrc*100.
nh_samp = nh_mc(nh_obs,nsamp)
;; number of Compton-thin sources
ithin = where(nh_samp lt 23.,nthin);,complement=ithick,ncomplement=nthick)
fthin = nthin/nsamp

sav_vars = ['NSAMP','NH_SAMP','NTHIN']
sav_inds = ['ITHIN']

;; number of iterations for each test
niter = 10000

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RUN FOR FIXED CTF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fraction and source numbers
fob = median(fobv)
nsr = round(nthin/(1.-fob))
nob = round(nsr*fob)
fct = median(fctv)
nct = round(nsr*fct)
f23 = median(f23v)
n23 = nob-nct;round(num*f23)
;; NH_RESAMP: structure of increased CT sources
;; NH_MOD: draw N sources from NH_RESAMP
;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where RX model > RX limits from Carroll+21
;; AD: Anderson-Darling test statistic    
nh_resamp = [nh_samp[ithin],23.+randomu(seed,n23),24.+2.*randomu(seed,nct)]
nresamp = n_elements(nh_resamp)
;; create arrays
nh_modv = dblarr(nsrc,niter)
rx_modv = dblarr(nsrc,niter)
iimodv = bytarr(nsrc,niter)
ad = dblarr(2,niter)
;; fill arrays
for n = 0,niter-1 do begin
    nh_modv[*,n] = nh_resamp[randomi(nsrc,nresamp)]
    rx_modv[*,n] = rx2nh(nh_modv[*,n],/rx_out,scat=rx_scat)
    iimodv[*,n] = rx_modv[*,n] gt rxl
    idet = where(iimodv[*,n] eq 1,detct)
    if (detct ge 5) then begin
        adtwo,rxd,rx_modv[idet,n],ad_stat,ad_crit
        ad[*,n] = [ad_stat,ad_crit]
    endif else if (detct gt 0) then begin
        ad[*,n] = [99.,99.]
    endif else message, 'NO MODELED DETECTIONS.'
endfor
rx_mod = hist2d_avg(rx_modv,0.2d,iidet=iimodv)
nh_mod = hist2d_avg(nh_modv,1d,iidet=iimodv)

sav_vars = [sav_vars,'NSR', $
                     'NOB','NCT','N23', $
                     'FOB','FCT','F23', $
                     'NH_MODV','RX_MODV','AD', $
                     'NH_MOD','RX_MOD']
sav_inds = [sav_inds,'IIMODV']


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RUN FOR SPLIT CTF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fraction and source numbers
fob_ = median(fobv1)
nsr_ = round(nthin/(1.-fob_))
nob_ = round(nsr_*fob_)
fct_ = median(fctv2)
nct_ = round(nsr_*fct_)
f23_ = median(f23v2)
n23_ = nob_-nct_;round(nsr_*f23_)
f25_ = median(f25v3)
n25_ = round(nsr_*f25_)
f24_ = median(f24v3)
n24_ = nct_-n25_;round(nsr_*f24_)
;; NH_RESAMP: structure of increased CT sources
;; NH_MOD: draw N sources from NH_RESAMP
;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where RX model > RX limits from Carroll+21
;; AD: Anderson-Darling test statistic    
nh_resamp = [nh_samp[ithin],23.+randomu(seed,n23_),24.+randomu(seed,n24_),25.+randomu(seed,n25_)]
nresamp = n_elements(nh_resamp)
;; create arrays
nh_modv_ = dblarr(nsrc,niter)
rx_modv_ = dblarr(nsrc,niter)
iimodv_ = bytarr(nsrc,niter)
ad_ = dblarr(2,niter)
;; fill arrays
for n = 0,niter-1 do begin
    nh_modv_[*,n] = nh_resamp[randomi(nsrc,nresamp)]
    rx_modv_[*,n] = rx2nh(nh_modv_[*,n],/rx_out,scat=rx_scat)
    iimodv_[*,n] = rx_modv_[*,n] gt rxl
    idet = where(iimodv_[*,n] eq 1,detct)
    if (detct ge 5) then begin
        adtwo,rxd,rx_modv_[idet,n],ad_stat,ad_crit
        ad_[*,n] = [ad_stat,ad_crit]
    endif else if (detct gt 0) then begin
        ad_[*,n] = [99.,99.]
    endif else message, 'NO MODELED DETECTIONS.'
endfor
nh_mod_ = hist2d_avg(nh_modv_,1d,iidet=iimodv_)
rx_mod_ = hist2d_avg(rx_modv_,0.2d,iidet=iimodv_)

sav_vars = [sav_vars,'NSR_', $
                     'NOB','NCT_','N23_','N24_','N25_', $
                     'FOB_','FCT_','F23_','F24_','F25_',$
                     'NH_MODV_','RX_MODV_','AD_', $
                     'NH_MOD_','RX_MOD_']
sav_inds = [sav_inds,'IIMODV_']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="rx_model.sav"')


END









