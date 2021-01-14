PRO estimate_split, RXZ = rxz, $
                    KCORR = kcorr


common _data
common _nhobs
common _ctfest

;; STDDEV observed in LX-LMIR relation of Chen+17
rx_scat = 0.2

;; separate WISE AGN, detections and non-detections
iixd = xdet ne '' and iiwac
iixn = xnon ne '' and iiwac
if keyword_set(kcorr) then rxd = rldet[where(iixd,ndet)]-alog10((1+z[where(iixd)])^(1.8-2.0)) else $
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
;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
yhks = cdf(ctf_ksv,bin=freedman(ctf_ksv),xloc=xhks)
yhad = cdf(ctf_adv,bin=freedman(ctf_adv),xloc=xhad)
div = 100.
ctf2sig_ks = [rnd(xhks[value_locate(yhks,0.05)],2):rnd(xhks[value_locate(yhks,0.95)],2):1./div]
ctf2sig_ad = [rnd(xhad[value_locate(yhad,0.05)],2):rnd(xhad[value_locate(yhad,0.95)],2):1./div]
ctf2sig = [ctf2sig_ks[0]<ctf2sig_ad[0]:ctf2sig_ks[-1]>ctf2sig_ad[-1]:1./div]
nfrac = n_elements(ctf2sig)

;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
div2 = 100.
ctf25 = dindgen(div2-1,start=1)/div2
ctf24 = 1.-ctf25
nsplit = n_elements(ctf24)



niter = 100
ctf24_ksv = dblarr(nfrac,niter)
ctf25_ksv = dblarr(nfrac,niter)
ctf24_adv = dblarr(nfrac,niter)
ctf25_adv = dblarr(nfrac,niter)

for n = 0,niter-1 do begin
    if n mod (niter/10) eq 0 then print, strtrim(n/(niter/100),2)+'% complete'

    ;print, 'MODEL_RXDIST FREE - 0% COMPLETE'
    for i = 0,nfrac-1 do begin
        nct = (nthin/(1.-ctf2sig[i]))*ctf2sig[i]
        tag2 = 'split'+string(rnd(ctf24[0]*100,0),format='(i02)')+'_'+string(rnd(ctf25[0]*100,0),format='(i02)')
        re = execute('nh_resamp2 = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct*ctf24[0]),25.+randomu(seed,nct*ctf25[0])]}')
        re = execute('nh_mod2 = {'+tag2+':(nh_resamp2.(0))[randomi(nsrc,n_elements(nh_resamp2.(0)))]}')
        if keyword_set(rxz) then re = execute('rx_mod2 = {'+tag2+':rx2nh_z(nh_mod2.(0),z[iwagn],/rx_out,scat=rx_scat)}') else $
                                 re = execute('rx_mod2 = {'+tag2+':rx2nh(nh_mod2.(0),/rx_out,scat=rx_scat)}')
        re = execute('iimod2 = {'+tag2+':rx_mod2.(0) gt rxl}')
        ks2 = dblarr(2,nsplit)
        kstwo,rxd,(rx_mod2.(0))[where(iimod2.(0) eq 1)],ks_stat,ks_prob
        ks2[*,0] = [ks_stat,ks_prob]
        ad2 = dblarr(2,nsplit)
        adtwo,rxd,(rx_mod2.(0))[where(iimod2.(0) eq 1)],ad_stat,ad_crit
        ad2[*,0] = [ad_stat,ad_crit]
        for j = 1,nsplit-1 do begin
            tag2 = 'split'+string(rnd(ctf24[j]*100,0),format='(i02)')+'_'+string(rnd(ctf25[j]*100,0),format='(i02)')
            nh_resamp2 = create_struct(nh_resamp2,tag2,[nh_samp[ithin],24.+randomu(seed,nct*ctf24[j]),25.+randomu(seed,nct*ctf25[j])])
            nh_mod2 = create_struct(nh_mod2,tag2,(nh_resamp2.(j))[randomi(nsrc,n_elements(nh_resamp2.(j)))])
            if keyword_set(rxz) then rx_mod2 = create_struct(rx_mod2,tag2,rx2nh_z(nh_mod2.(j),z[iwagn],/rx_out,scat=rx_scat)) else $
                                     rx_mod2 = create_struct(rx_mod2,tag2,rx2nh(nh_mod2.(j),/rx_out,scat=rx_scat))
            iimod2 = create_struct(iimod2,tag2,rx_mod2.(j) gt rxl)
            kstwo,rxd,(rx_mod2.(j))[where(iimod2.(j) eq 1)],ks_stat,ks_prob
            ks2[*,j] = [ks_stat,ks_prob]
            adtwo,rxd,(rx_mod2.(j))[where(iimod2.(j) eq 1)],ad_stat,ad_crit
            ad2[*,j] = [ad_stat,ad_crit]
        endfor
        ksm = min(ks2[0,*],imin)
        ctf24_ksv[i,n] = ctf24[imin]
        ctf25_ksv[i,n] = ctf25[imin]
        adm = min(ad2[0,*],imin)
        ctf24_adv[i,n] = ctf24[imin]
        ctf25_adv[i,n] = ctf25[imin]
    endfor
endfor

;; save full arrays and reduce to 2sigma per run
ctf24_ksv_full = ctf24_ksv
ctf25_ksv_full = ctf25_ksv
ictf_ks = [value_locate(ctf2sig,ctf2sig_ks[0]):value_locate(ctf2sig,ctf2sig_ks[-1]):1]
ctf24_ksv = ctf24_ksv[ictf_ks,*]
ctf25_ksv = ctf25_ksv[ictf_ks,*]

ctf24_adv_full = ctf24_adv
ctf25_adv_full = ctf25_adv
ictf_ad = [value_locate(ctf2sig,ctf2sig_ad[0]):value_locate(ctf2sig,ctf2sig_ad[-1]):1]
ctf24_adv = ctf24_adv[ictf_ad,*]
ctf25_adv = ctf25_adv[ictf_ad,*]

sav_vars = ['CTF2SIG','CTF24_KSV_FULL','CTF25_KSV_FULL','CTF24_ADV_FULL','CTF25_ADV_FULL', $
            'XHKS','YHKS','CTF2SIG_KS','CTF24_KSV','CTF25_KSV', $
            'XHAD','YHAD','CTF2SIG_AD','CTF24_ADV','CTF25_ADV']            
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="split_estimate.sav"')


END










