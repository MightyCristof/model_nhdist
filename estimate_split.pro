PRO estimate_split


common _data
common _nhobs
common _rxnh
common _ctfest


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


;; RUN FOR KS TEST

;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
yhks = edf(ctf_ksv,xloc=xhks)
div = 100.d
;; all, 3-,2-,1-sigma
;ctf_ks = minmax(ctf_ksv)
ctf_ks = [xhks[value_locate(yhks,0.01)]:xhks[value_locate(yhks,0.99)]:1./div]
;ctf_ks = [xhks[value_locate(yhks,0.05)]:xhks[value_locate(yhks,0.95)]:1./div]
;ctf_ks = [xhks[value_locate(yhks,0.32)]:xhks[value_locate(yhks,0.68)]:1./div]
;; don't loops twice, just do it once
nfrac_ks = n_elements(ctf_ks)

;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
div2 = 100.
ctf25 = dindgen(div2-1,start=1)/div2
ctf24 = 1.-ctf25
nsplit = n_elements(ctf24)

niter = 1000
ct24_ksv = dblarr(nfrac_ks,niter)
ct24_ksv = dblarr(nfrac_ks,niter)
ks2v = dblarr(nfrac_ks,niter)

for n = 0,niter-1 do begin
    if n mod (niter/10) eq 0 then print, strtrim(n/(niter/100),2)+'% complete'

    ;print, 'MODEL_RXDIST FREE - 0% COMPLETE'
    for i = 0,nfrac_ks-1 do begin
        nct = (nthin/(1.-ctf_ks[i]))*ctf_ks[i]
        tag2 = 'split'+string(rnd(ctf24[0]*100,0),format='(i02)')+'_'+string(rnd(ctf25[0]*100,0),format='(i02)')
        re = execute('nh_resamp2 = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct*ctf24[0]),25.+randomu(seed,nct*ctf25[0])]}')
        re = execute('nh_mod2 = {'+tag2+':(nh_resamp2.(0))[randomi(nsrc,n_elements(nh_resamp2.(0)))]}')
        re = execute('rx_mod2 = {'+tag2+':rx2nh(nh_mod2.(0),/rx_out,scat=rx_scat)}')
        re = execute('iimod2 = {'+tag2+':rx_mod2.(0) gt rxwl}')
        ks2 = dblarr(2,nsplit)
        kstwo,rxwd,(rx_mod2.(0))[where(iimod2.(0) eq 1)],ks_stat,ks_prob
        ks2[*,0] = [ks_stat,ks_prob]
        for j = 1,nsplit-1 do begin
            tag2 = 'split'+string(rnd(ctf24[j]*100,0),format='(i02)')+'_'+string(rnd(ctf25[j]*100,0),format='(i02)')
            nh_resamp2 = create_struct(nh_resamp2,tag2,[nh_samp[ithin],24.+randomu(seed,nct*ctf24[j]),25.+randomu(seed,nct*ctf25[j])])
            nh_mod2 = create_struct(nh_mod2,tag2,(nh_resamp2.(j))[randomi(nsrc,n_elements(nh_resamp2.(j)))])
            rx_mod2 = create_struct(rx_mod2,tag2,rx2nh(nh_mod2.(j),/rx_out,scat=rx_scat))
            iimod2 = create_struct(iimod2,tag2,rx_mod2.(j) gt rxwl)
            kstwo,rxwd,(rx_mod2.(j))[where(iimod2.(j) eq 1)],ks_stat,ks_prob
            ks2[*,j] = [ks_stat,ks_prob]
        endfor
        ksm = min(ks2[0,*],imin)
        ks2v[i,n] = ks2[0,imin]
        ct24_ksv[i,n] = ctf24[imin]
        ct24_ksv[i,n] = ctf25[imin]
    endfor
endfor

sav_vars = ['XHKS','YHKS','CTF_KS','CT24_KSV','CT24_KSV','KS2V']            
sav_inds = []


;; RUN FOR KS TEST

;; range from NH=24 to end of nh_samp
;nh_diff = width([max(nh_samp),24.])
;; range from NH=24-26
nh_diff = 2.
;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
yhad = edf(ctf_adv,xloc=xhad)
div = 100.d
;; all, 3-,2-,1-sigma
;ctf_ad = minmax(ctf_adv)
ctf_ad = [xhad[value_locate(yhad,0.01)]:xhad[value_locate(yhad,0.99)]:1./div]
;ctf_ad = [xhad[value_locate(yhad,0.05)]:xhad[value_locate(yhad,0.95)]:1./div]
;ctf_ad = [xhad[value_locate(yhad,0.32)]:xhad[value_locate(yhad,0.68)]:1./div]
;; don't loops twice, just do it once
nfrac_ad = n_elements(ctf_ad)

;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
div2 = 100.
ctf25 = dindgen(div2-1,start=1)/div2
ctf24 = 1.-ctf25
nsplit = n_elements(ctf24)

niter = 1000
ct24_adv = dblarr(nfrac_ad,niter)
ct24_adv = dblarr(nfrac_ad,niter)
ad2v = dblarr(nfrac_ad,niter)

for n = 0,niter-1 do begin
    if n mod (niter/10) eq 0 then print, strtrim(n/(niter/100),2)+'% complete'

    ;print, 'MODEL_RXDIST FREE - 0% COMPLETE'
    for i = 0,nfrac_ad-1 do begin
        nct = (nthin/(1.-ctf_ad[i]))*ctf_ad[i]
        tag2 = 'split'+string(rnd(ctf24[0]*100,0),format='(i02)')+'_'+string(rnd(ctf25[0]*100,0),format='(i02)')
        re = execute('nh_resamp2 = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct*ctf24[0]),25.+randomu(seed,nct*ctf25[0])]}')
        re = execute('nh_mod2 = {'+tag2+':(nh_resamp2.(0))[randomi(nsrc,n_elements(nh_resamp2.(0)))]}')
        re = execute('rx_mod2 = {'+tag2+':rx2nh(nh_mod2.(0),/rx_out,scat=rx_scat)}')
        re = execute('iimod2 = {'+tag2+':rx_mod2.(0) gt rxwl}')
        ad2 = dblarr(2,nsplit)
        adtwo,rxwd,(rx_mod2.(0))[where(iimod2.(0) eq 1)],ad_stat,ad_crit
        ad2[*,0] = [ad_stat,ad_crit]
        for j = 1,nsplit-1 do begin
            tag2 = 'split'+string(rnd(ctf24[j]*100,0),format='(i02)')+'_'+string(rnd(ctf25[j]*100,0),format='(i02)')
            nh_resamp2 = create_struct(nh_resamp2,tag2,[nh_samp[ithin],24.+randomu(seed,nct*ctf24[j]),25.+randomu(seed,nct*ctf25[j])])
            nh_mod2 = create_struct(nh_mod2,tag2,(nh_resamp2.(j))[randomi(nsrc,n_elements(nh_resamp2.(j)))])
            rx_mod2 = create_struct(rx_mod2,tag2,rx2nh(nh_mod2.(j),/rx_out,scat=rx_scat))
            iimod2 = create_struct(iimod2,tag2,rx_mod2.(j) gt rxwl)
            adtwo,rxwd,(rx_mod2.(j))[where(iimod2.(j) eq 1)],ad_stat,ad_crit
            ad2[*,j] = [ad_stat,ad_crit]
        endfor
        adm = min(ad2[0,*],imin)
        ad2v[i,n] = ad2[0,imin]
        ct24_adv[i,n] = ctf24[imin]
        ct24_adv[i,n] = ctf25[imin]
    endfor
endfor

sav_vars = [sav_vars,'XHAD','YHAD','CTF_AD','CT24_ADV','CT24_ADV','AD2V']            
sav_inds = []


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="split_estimate.sav"')


END










