PRO model_rxdist, RXZ = rxz


common _data
common _nhobs
common _ctfest
common _split

;; STDDEV observed in LX-LMIR relation of Chen+17
rx_scat = 0.2

;; separate WISE AGN, detections and non-detections
iixd = xdet ne '' and iiwac
iixn = xnon ne '' and iiwac
if keyword_set(rxz) then rxd = rldet[where(iixd,ndet)]-alog10((1+z[where(iixd)])^(1.8-2.0)) else $
                         rxd = rldet[where(iixd,ndet)]
e_rxd = e_rldet[where(iixd)]
iwagn = where(iiwac,nsrc)
rxl = rxlim[iwagn]

;; resample observed NH distribution to increase data density
nsamp = nsrc*100.
nh_samp = nh_mc(nh_lan_nst,nsamp)
;; number of CT sources
ithin = where(nh_samp le 24.,nthin);,complement=ithick,ncomplement=nthick)
fthin = nthin/nsamp
;; range from NH=24 to end of nh_samp
;nh_diff = width([max(nh_samp),24.])
;; range from NH=24-26
nh_diff = 2.

;; run for KS test

;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
;; reduce KS/AD to 2sigma separately (run was combined in previous step for computational time)
ifrac_ks = value_locate(ctf2sig_ks,median(ctf_ksv))
ctf_ks = ctf2sig_ks[ifrac_ks]
;; NH_RESAMP: structure of increased CT sources. scale1_0 == nh_samp
;; NH_RESAMP has different number of sources per tag, cannot use SOA2AOS
;; NH_MOD: draw (from NH_RESAMP) a number of sources == X-ray non-detections
;; RX_MOD: convert the NH model to RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where rx_mod > RX non-detection from Carroll+20
;; KS: KS test between RX model detections and RX detections from Carroll+20
nct_ks = (nthin/(1.-ctf_ks))*ctf_ks
nh_resamp_ks = [nh_samp[ithin],24.+randomu(seed,nct_ks)*nh_diff]
nh_mod_ks = (nh_resamp_ks)[randomi(nsrc,n_elements(nh_resamp_ks))]
if keyword_set(rxz) then rx_mod_ks = rx2nh_z(nh_mod_ks,z[iwagn],/rx_out,scat=rx_scat) else $
                         rx_mod_ks = rx2nh(nh_mod_ks,/rx_out,scat=rx_scat)
iimod_ks = rx_mod_ks gt rxl
kstwo,rxd,rx_mod_ks[where(iimod_ks eq 1)],ks_stat,ks_prob
ks = [ks_stat,ks_prob]
iks = 0

sav_vars = ['RX_SCAT','RXD','E_RXD','RXL', $
            'NSAMP','NH_SAMP','NTHIN','NH_DIFF', $
            'CTF_KS','NH_RESAMP_KS','NH_MOD_KS','RX_MOD_KS','KS']
sav_inds = ['IIXD','IIXN','ITHIN','IIMOD_KS','IKS']


;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
ctf24_ks = median(ctf24_ksv[ifrac_ks,*])
ctf25_ks = 1.-ctf24_ks;rnd(median(ctf25_ksv[ifrac_ks,*]),2)

niter = 1000
fmt2 = 'i0'+strtrim(strlen(strtrim(niter,2)),2)

print, 'KS: MODEL_RXDIST FREE - 0% COMPLETE'
re = execute('tag2 = "iter"+string(0,format="('+fmt2+')")')
re = execute('nh_resamp2_ks = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct_ks*ctf24_ks[0]),25.+randomu(seed,nct_ks*ctf25_ks[0])]}')
re = execute('nh_mod2_ks = {'+tag2+':(nh_resamp2_ks.(0))[randomi(nsrc,n_elements(nh_resamp2_ks.(0)))]}')
if keyword_set(rxz) then re = execute('rx_mod2_ks = {'+tag2+':rx2nh_z(nh_mod2_ks.(0),z[iwagn],/rx_out,scat=rx_scat)}') else $
                         re = execute('rx_mod2_ks = {'+tag2+':rx2nh(nh_mod2_ks.(0),/rx_out,scat=rx_scat)}')
re = execute('iimod2_ks = {'+tag2+':rx_mod2_ks.(0) gt rxl}')
ks2 = dblarr(2,niter)
kstwo,rxd,(rx_mod2_ks.(0))[where(iimod2_ks.(0) eq 1)],ks_stat,ks_prob
ks2[*,0] = [ks_stat,ks_prob]
for j = 1,niter-1 do begin
    re = execute('tag2 = "iter"+string(j,format="('+fmt2+')")')
    nh_resamp2_ks = create_struct(nh_resamp2_ks,tag2,[nh_samp[ithin],24.+randomu(seed,nct_ks*ctf24_ks),25.+randomu(seed,nct_ks*ctf25_ks)])
    nh_mod2_ks = create_struct(nh_mod2_ks,tag2,(nh_resamp2_ks.(j))[randomi(nsrc,n_elements(nh_resamp2_ks.(j)))])
    if keyword_set(rxz) then rx_mod2_ks = create_struct(rx_mod2_ks,tag2,rx2nh_z(nh_mod2_ks.(j),z[iwagn],/rx_out,scat=rx_scat)) else $
                             rx_mod2_ks = create_struct(rx_mod2_ks,tag2,rx2nh(nh_mod2_ks.(j),/rx_out,scat=rx_scat))
    iimod2_ks = create_struct(iimod2_ks,tag2,rx_mod2_ks.(j) gt rxl)
    kstwo,rxd,(rx_mod2_ks.(j))[where(iimod2_ks.(j) eq 1)],ks_stat,ks_prob
    ks2[*,j] = [ks_stat,ks_prob]
    if (j mod (niter/2.) eq 0) then print, 'KS: MODEL_RXDIST FREE - 50% COMPLETE'
endfor
print, 'KS: MODEL_RXDIST FREE - 100% COMPLETE'

ks2med = median(ks2[0,*])
isort = sort(ks2[0,*])
iloc = value_locate(ks2[0,isort],ks2med)
iks2 = isort[iloc]


sav_vars = [sav_vars,'CTF24_KS','CTF25_KS', $
                     'NH_RESAMP2_KS','NH_MOD2_KS','RX_MOD2_KS','KS2']
sav_inds = [sav_inds,'IIMOD2_KS','IKS2']


;; run for AD test

;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
;; reduce KS/AD to 2sigma separately (run was combined in previous step for computational time)
ifrac_ad = value_locate(ctf2sig_ad,median(ctf_adv))
ctf_ad = ctf2sig_ad[ifrac_ad]
;; NH_RESAMP: structure of increased CT sources. scale1_0 == nh_samp
;; NH_RESAMP has different number of sources per tag, cannot use SOA2AOS
;; NH_MOD: draw (from NH_RESAMP) a number of sources == X-ray non-detections
;; RX_MOD: convert the NH model to RX with observed scatter in Chen+17 LX-LMIR
;; IIMOD: where rx_mod > RX non-detection from Carroll+20
;; KS: KS test between RX model detections and RX detections from Carroll+20
nct_ad = (nthin/(1.-ctf_ad))*ctf_ad
nh_resamp_ad = [nh_samp[ithin],24.+randomu(seed,nct_ad)*nh_diff]
nh_mod_ad = (nh_resamp_ad)[randomi(nsrc,n_elements(nh_resamp_ad))]
if keyword_set(rxz) then rx_mod_ad = rx2nh_z(nh_mod_ad,z[iwagn],/rx_out,scat=rx_scat) else $
                         rx_mod_ad = rx2nh(nh_mod_ad,/rx_out,scat=rx_scat)
iimod_ad = rx_mod_ad gt rxl
adtwo,rxd,rx_mod_ad[where(iimod_ad eq 1)],ad_stat,ad_crit
ad = [ad_stat,ad_crit]
iad = 0

sav_vars = [sav_vars,'CTF_AD','NH_RESAMP_AD','NH_MOD_AD','RX_MOD_AD','AD']
sav_inds = [sav_inds,'IIMOD_AD','IAD']


;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
ctf24_ad = median(ctf24_adv[ifrac_ad,*])
ctf25_ad = 1.-ctf24_ad;rnd(median(ctf25_adv[ifrac_ad,*]),2)

niter = 1000
fmt2 = 'i0'+strtrim(strlen(strtrim(niter,2)),2)

print, 'AD: MODEL_RXDIST FREE - 0% COMPLETE'
re = execute('tag2 = "iter"+string(0,format="('+fmt2+')")')
re = execute('nh_resamp2_ad = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct_ad*ctf24_ad[0]),25.+randomu(seed,nct_ad*ctf25_ad[0])]}')
re = execute('nh_mod2_ad = {'+tag2+':(nh_resamp2_ad.(0))[randomi(nsrc,n_elements(nh_resamp2_ad.(0)))]}')
if keyword_set(rxz) then re = execute('rx_mod2_ad = {'+tag2+':rx2nh_z(nh_mod2_ad.(0),z[iwagn],/rx_out,scat=rx_scat)}') else $
                         re = execute('rx_mod2_ad = {'+tag2+':rx2nh(nh_mod2_ad.(0),/rx_out,scat=rx_scat)}')
re = execute('iimod2_ad = {'+tag2+':rx_mod2_ad.(0) gt rxl}')
ad2 = dblarr(2,niter)
adtwo,rxd,(rx_mod2_ad.(0))[where(iimod2_ad.(0) eq 1)],ad_stat,ad_crit
ad2[*,0] = [ad_stat,ad_crit]
for j = 1,niter-1 do begin
    re = execute('tag2 = "iter"+string(j,format="('+fmt2+')")')
    nh_resamp2_ad = create_struct(nh_resamp2_ad,tag2,[nh_samp[ithin],24.+randomu(seed,nct_ad*ctf24_ad),25.+randomu(seed,nct_ad*ctf25_ad)])
    nh_mod2_ad = create_struct(nh_mod2_ad,tag2,(nh_resamp2_ad.(j))[randomi(nsrc,n_elements(nh_resamp2_ad.(j)))])
    if keyword_set(rxz) then rx_mod2_ad = create_struct(rx_mod2_ad,tag2,rx2nh_z(nh_mod2_ad.(j),z[iwagn],/rx_out,scat=rx_scat)) else $
                             rx_mod2_ad = create_struct(rx_mod2_ad,tag2,rx2nh(nh_mod2_ad.(j),/rx_out,scat=rx_scat))
    iimod2_ad = create_struct(iimod2_ad,tag2,rx_mod2_ad.(j) gt rxl)
    adtwo,rxd,(rx_mod2_ad.(j))[where(iimod2_ad.(j) eq 1)],ad_stat,ad_crit
    ad2[*,j] = [ad_stat,ad_crit]
    if (j mod (niter/2.) eq 0) then print, 'AD: MODEL_RXDIST FREE - 50% COMPLETE'
endfor
print, 'AD: MODEL_RXDIST FREE - 100% COMPLETE'

ad2med = median(ad2[0,*])
isort = sort(ad2[0,*])
iloc = value_locate(ad2[0,isort],ad2med)
iad2 = isort[iloc]


sav_vars = [sav_vars,'CTF24_AD','CTF25_AD', $
                     'NH_RESAMP2_AD','NH_MOD2_AD','RX_MOD2_AD','AD2']
sav_inds = [sav_inds,'IIMOD2_AD','IAD2']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="rx_model.sav"')


END









