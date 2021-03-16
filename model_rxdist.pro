PRO model_rxdist


common _data
common _nhobs
common _rxnh
common _group
common _ctfest
common _split


;; resample observed NH distribution to increase data density
nsamp = nsrc*100.
nh_samp = nh_mc(nh_lan,nsamp)
;; number of CT sources
ithin = where(nh_samp le 24.,nthin);,complement=ithick,ncomplement=nthick)
fthin = nthin/nsamp
;; range from NH=24 to end of nh_samp
;nh_diff = width([max(nh_samp),24.])
;; range from NH=24-26
nh_diff = 2.

sav_vars = ['NSAMP','NH_SAMP','NTHIN','NH_DIFF']
sav_inds = ['ITHIN']

;; number of iterations for each test
niter = 10000
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RUN FOR FIXED CTF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; RUN FOR KS TEST

;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
;; reduce KS/AD to 2sigma separately (run was combined in previous step for computational time)
ictf_ksx = value_locate(ctf_ks,median(ctf_ksv))
ctf_ksx = ctf_ks[ictf_ksx]
nct_ks = (nthin/(1.-ctf_ksx))*ctf_ksx
;; NH_RESAMP: structure of increased CT sources. scale1_0 == nh_samp
;; NH_RESAMP has different number of sources per tag, cannot use SOA2AOS
;; NH_MOD: draw (from NH_RESAMP) a number of sources == X-ray non-detections
;; RX_MOD: convert the NH model to RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where rx_mod > RX non-detection from Carroll+20
;; KS: KS test between RX model detections and RX detections from Carroll+20
nh_resamp_ks = [nh_samp[ithin],24.+randomu(seed,nct_ks)*nh_diff]
nresamp = n_elements(nh_resamp_ks)
;; create arrays
nh_mod_ksv = dblarr(nsrc,niter)
rx_mod_ksv = dblarr(nsrc,niter)
iimod_ksv = bytarr(nsrc,niter)
ks = dblarr(2,niter)
;; fill arrays
for i = 0,niter-1 do begin
    nh_mod_ksv[*,i] = nh_resamp_ks[randomi(nsrc,nresamp)]
    rx_mod_ksv[*,i] = rx2nh(nh_mod_ksv[*,i],/rx_out,scat=rx_scat)
    iimod_ksv[*,i] = rx_mod_ksv[*,i] gt rxl
    idet = where(iimod_ksv[*,i] eq 1,detct)
    if (detct ge 5) then begin
        kstwo,rxd,rx_mod_ksv[idet,i],ks_stat,ks_prob
        ks[*,i] = [ks_stat,ks_prob]
    endif else if (detct gt 0) then begin
        ks[*,i] = [99.,99.]
    endif else message, 'NO MODELED DETECTIONS.'
endfor
rx_mod_ks = hist2d_avg(rx_mod_ksv,0.2d,iidet=iimod_ksv)
nh_mod_ks = hist2d_avg(nh_mod_ksv,1d,iidet=iimod_ksv)

sav_vars = [sav_vars,'CTF_KSX','NH_MOD_KSV','RX_MOD_KSV','KS','NH_MOD_KS','RX_MOD_KS']
sav_inds = [sav_inds,'ICTF_KSX','IIMOD_KSV']

;; RUN FOR AD TEST
;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
;; reduce KS/AD to 2sigma separately (run was combined in previous step for computational time)
ictf_adx = value_locate(ctf_ad,median(ctf_adv))
ctf_adx = ctf_ad[ictf_adx]
nct_ad = (nthin/(1.-ctf_adx))*ctf_adx
;; NH_RESAMP: structure of increased CT sources. scale1_0 == nh_samp
;; NH_RESAMP has different number of sources per tag, cannot use SOA2AOS
;; NH_MOD: draw (from NH_RESAMP) a number of sources == X-ray non-detections
;; RX_MOD: convert the NH model to RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where rx_mod > RX non-detection from Carroll+20
;; KS: KS test between RX model detections and RX detections from Carroll+20
nh_resamp_ad = [nh_samp[ithin],24.+randomu(seed,nct_ad)*nh_diff]
nresamp = n_elements(nh_resamp_ad)
;; create arrays
nh_mod_adv = dblarr(nsrc,niter)
rx_mod_adv = dblarr(nsrc,niter)
iimod_adv = bytarr(nsrc,niter)
ad = dblarr(2,niter)
;; fill arrays
for i = 0,niter-1 do begin
    nh_mod_adv[*,i] = nh_resamp_ad[randomi(nsrc,nresamp)]
    rx_mod_adv[*,i] = rx2nh(nh_mod_adv[*,i],/rx_out,scat=rx_scat)
    iimod_adv[*,i] = rx_mod_adv[*,i] gt rxl
    idet = where(iimod_adv[*,i] eq 1,detct)
    if (detct ge 5) then begin
        kstwo,rxd,rx_mod_adv[idet,i],ks_stat,ks_prob
        ks[*,i] = [ks_stat,ks_prob]
    endif else if (detct gt 0) then begin
        ks[*,i] = [99.,99.]
    endif else message, 'NO MODELED DETECTIONS.'
endfor
rx_mod_ad = hist2d_avg(rx_mod_adv,0.2d,iidet=iimod_adv)
nh_mod_ad = hist2d_avg(nh_mod_adv,1d,iidet=iimod_adv)

sav_vars = [sav_vars,'CTF_ADX','NH_MOD_ADV','RX_MOD_ADV','AD','NH_MOD_AD','RX_MOD_AD']
sav_inds = [sav_inds,'ICTF_ADX','IIMOD_ADV']



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RUN FOR SPLIT CTF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; RUN FOR KS TEST

min_ks2 = min(median(ks2v,dim=2),ictf2_ksx)
ctf2_ks = ctf_ks[ictf2_ksx]

;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
ct24_ks = median(ct24_ksv[ictf2_ksx,*])
ct25_ks = 1.-ct24_ks

nh_resamp2_ks = [nh_samp[ithin],24.+randomu(seed,nct_ks*ct24_ks),25.+randomu(seed,nct_ks*ct25_ks)]
nresamp = n_elements(nh_resamp2_ks)
;; create arrays
nh_mod2_ksv = dblarr(nsrc,niter)
rx_mod2_ksv = dblarr(nsrc,niter)
iimod2_ksv = bytarr(nsrc,niter)
ks2 = dblarr(2,niter)
;; fill arrays
for i = 0,niter-1 do begin
    nh_mod2_ksv[*,i] = nh_resamp2_ks[randomi(nsrc,nresamp)]
    rx_mod2_ksv[*,i] = rx2nh(nh_mod2_ksv[*,i],/rx_out,scat=rx_scat)
    iimod2_ksv[*,i] = rx_mod2_ksv[*,i] gt rxl
    idet = where(iimod2_ksv[*,i] eq 1,detct)
    if (detct ge 5) then begin
        kstwo,rxd,rx_mod2_ksv[idet,i],ks_stat,ks_prob
        ks2[*,i] = [ks_stat,ks_prob]
    endif else if (detct gt 0) then begin
        ks2[*,i] = [99.,99.]
    endif else message, 'NO MODELED DETECTIONS.'
endfor
nh_mod2_ks = hist2d_avg(nh_mod2_ksv,1d,iidet=iimod2_ksv)
rx_mod2_ks = hist2d_avg(rx_mod2_ksv,0.2d,iidet=iimod2_ksv)

sav_vars = [sav_vars,'CTF2_KS','CT24_KS','CT25_KS','NH_MOD2_KSV','RX_MOD2_KSV','KS2', $
                     'NH_MOD2_KS','RX_MOD2_KS']
sav_inds = [sav_inds,'IIMOD2_KSV']


;; RUN FOR AD TEST

min_ad2 = min(median(ad2v,dim=2),ictf2_adx)
ctf2_ad = ctf_ad[ictf2_adx]

;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
ct24_ad = median(ct24_adv[ictf2_adx,*])
ct25_ad = 1.-ct24_ad

nh_resamp2_ad = [nh_samp[ithin],24.+randomu(seed,nct_ad*ct24_ad),25.+randomu(seed,nct_ad*ct25_ad)]
nresamp = n_elements(nh_resamp2_ad)
;; create arrays
nh_mod2_adv = dblarr(nsrc,niter)
rx_mod2_adv = dblarr(nsrc,niter)
iimod2_adv = bytarr(nsrc,niter)
ad2 = dblarr(2,niter)
;; fill arrays
for i = 0,niter-1 do begin
    nh_mod2_adv[*,i] = nh_resamp2_ad[randomi(nsrc,nresamp)]
    rx_mod2_adv[*,i] = rx2nh(nh_mod2_adv[*,i],/rx_out,scat=rx_scat)
    iimod2_adv[*,i] = rx_mod2_adv[*,i] gt rxl
    idet = where(iimod2_adv[*,i] eq 1,detct)
    if (detct ge 5) then begin
        adtwo,rxd,rx_mod2_adv[idet,i],ad_stat,ad_crit
        ad2[*,i] = [ad_stat,ad_crit]
    endif else if (detct gt 0) then begin
        ad2[*,i] = [99.,99.]
    endif else message, 'NO MODELED DETECTIONS.'
endfor
nh_mod2_ad = hist2d_avg(nh_mod2_adv,1d,iidet=iimod2_adv)
rx_mod2_ad = hist2d_avg(rx_mod2_adv,0.2d,iidet=iimod2_adv)

sav_vars = [sav_vars,'CTF2_AD','CT24_AD','CT25_AD','NH_MOD2_ADV','RX_MOD2_ADV','AD2', $
                     'NH_MOD2_AD','RX_MOD2_AD']
sav_inds = [sav_inds,'IIMOD2_ADV']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="rx_model.sav"')


END









