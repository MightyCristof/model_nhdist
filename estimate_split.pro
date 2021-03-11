PRO estimate_split


common _data
common _nhobs
common _rxnh
common _group
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

;; scale the number of CT sources
;; consider fraction of CT sources rather than arbitrary number
yhks = edf(ctf_ksv,xloc=xhks)
yhad = edf(ctf_adv,xloc=xhad)
div = 100.d
;; all, 3-,2-,1-sigma
;ctf_ks = minmax(ctf_ksv)
ctf_ks = [xhks[value_locate(yhks,0.01)]:xhks[value_locate(yhks,0.99)]:1./div]
;ctf_ks = [xhks[value_locate(yhks,0.05)]:xhks[value_locate(yhks,0.95)]:1./div]
;ctf_ks = [xhks[value_locate(yhks,0.32)]:xhks[value_locate(yhks,0.68)]:1./div]
;ctf_ad = minmax(ctf_adv)
ctf_ad = [xhad[value_locate(yhad,0.01)]:xhad[value_locate(yhad,0.99)]:1./div]
;ctf_ad = [xhad[value_locate(yhad,0.05)]:xhad[value_locate(yhad,0.95)]:1./div]
;ctf_ad = [xhad[value_locate(yhad,0.32)]:xhad[value_locate(yhad,0.68)]:1./div]
;; don't do it twice, just loop through once
nfrac_ks = n_elements(ctf_ks)
nfrac_ad = n_elements(ctf_ad)
ctf = [ctf_ks,ctf_ad]
ctf = ctf[uniq(ctf,sort(ctf))]
nfrac = n_elements(ctf)

;; now do the same as before, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
div2 = 100.d
ct25 = dindgen(div2-1,start=1)/div2
ct24 = 1.-ct25
nsplit = n_elements(ct24)

niter = 1000
ct24_ksv = dblarr(nfrac_ks,niter)
ct25_ksv = dblarr(nfrac_ks,niter)
ks2v = dblarr(nfrac_ks,niter)
ct24_adv = dblarr(nfrac_ad,niter)
ct25_adv = dblarr(nfrac_ad,niter)
ad2v = dblarr(nfrac_ad,niter)

for n = 0,niter-1 do begin
    if n mod (niter/10) eq 0 then print, strtrim(n/(niter/100),2)+'% complete'

    for i = 0,nfrac-1 do begin
        nh_mod2 = dblarr(nsrc,nsplit)
        rx_mod2 = dblarr(nsrc,nsplit)
        iimod2 = bytarr(nsrc,nsplit)
        ks2 = dblarr(2,nsplit)
        ad2 = dblarr(2,nsplit)
        nct = round((nthin/(1.-ctf[i]))*ctf[i])
        for j = 0,nsplit-1 do begin
            nh_resamp2 = [nh_samp[ithin],24.+randomu(seed,nct*ct24[j]),25.+randomu(seed,nct*ct25[j])]
            nresamp2 = n_elements(nh_resamp2)
            nh_mod2[*,j] = nh_resamp2[randomi(nsrc,nresamp2)]
            rx_mod2[*,j] = rx2nh(nh_mod2[*,j],/rx_out,scat=rx_scat)
            iimod2[*,j] = rx_mod2[*,j] gt rxl
            idet = where(iimod2[*,j] eq 1,detct)
            if (detct ge 5) then begin
                kstwo,rxd,rx_mod2[idet,j],ks_stat,ks_prob
                ks2[*,j] = [ks_stat,ks_prob]
                adtwo,rxd,rx_mod2[idet,j],ad_stat,ad_crit
                ad2[*,j] = [ad_stat,ad_crit]
            endif else if (detct gt 0) then begin
                ks2[*,i] = [99.,99.]
                ad2[*,i] = [99.,99.]
            endif else message, 'NO MODELED DETECTIONS.'
        endfor
        ctf_str = strtrim(ctf[i],2)
        iks = where(ctf_str eq strtrim(ctf_ks,2),ks_match)
        ;if (ks_match eq 0) then print, ctf[i]
        if (ks_match eq 1) then begin
            ksm = min(ks2[0,*],imin)
            ks2v[iks,n] = ks2[0,imin]
            ct24_ksv[iks,n] = ct24[imin]
            ct25_ksv[iks,n] = ct25[imin]
            ;print, ctf[i], ' ', ctf_ks[iks]
        endif
        iad = where(ctf_str eq strtrim(ctf_ad,2),ad_match)
        ;if (ad_match eq 0) then print, ctf[i]
        if (ad_match eq 1) then begin
            adm = min(ad2[0,*],imin)
            ad2v[iad,n] = ad2[0,imin]
            ct24_adv[iad,n] = ct24[imin]
            ct25_adv[iad,n] = ct25[imin]
            ;print, ctf[i], ' ', ctf_ad[iad]
        endif
    endfor
endfor

sav_vars = ['XHKS','YHKS','CTF_KS','CT24_KSV','CT25_KSV','KS2V', $
            'XHAD','YHAD','CTF_AD','CT24_ADV','CT25_ADV','AD2V']            
sav_inds = []


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="split_estimate.sav"')


END










