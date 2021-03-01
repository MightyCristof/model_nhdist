PRO estimate_ctf


common _data
common _nhobs
common _rxnh
common _group


;; run this script NITER times and look at the distribution in CTF
niter = 10000
ctf_ksv = dblarr(niter)
ctf_adv = dblarr(niter)
ksv = dblarr(niter)
adv = dblarr(niter)

for n = 0,niter-1 do begin
    if n mod (niter/10) eq 0 then print, strtrim(n/(niter/100),2)+'% complete'

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
    div = 100.
    ctf = dindgen(div-1,start=1)/div
    nfrac = n_elements(ctf)
    ;; NH_RESAMP: structure of increased CT sources. scale1_0 == nh_samp
    ;; NH_RESAMP has different number of sources per tag, cannot use SOA2AOS
    ;; NH_MOD: draw (from NH_RESAMP) a number of sources == X-ray non-detections
    ;; RX_MOD: convert the NH model to RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where rx_mod > RX non-detection from Carroll+20
    ;; KS: KS test between RX model detections and RX detections from Carroll+20
    ;; AD: Anderson-Darling test statistic    
    nh_mod = dblarr(nsrc,nfrac)
    rx_mod = dblarr(nsrc,nfrac)
    iimod = bytarr(nsrc,nfrac)
    ks = dblarr(2,nfrac)
    ad = dblarr(2,nfrac)
    
    for i = 0,nfrac-1 do begin
        nct = round((nthin/(1.-ctf[i]))*ctf[i])
        nh_resamp = [nh_samp[ithin],24.+randomu(seed,nct)*nh_diff]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,i] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,i] = rx2nh(nh_mod[*,i],/rx_out,scat=rx_scat)
        iimod[*,i] = rx_mod[*,i] gt rxl
        idet = where(iimod[*,i] eq 1,detct)
        if (detct gt 0) then begin
            kstwo,rxd,rx_mod[idet,i],ks_stat,ks_prob
            ks[*,i] = [ks_stat,ks_prob]
            adtwo,rxd,rx_mod[idet,i],ad_stat,ad_crit
            ad[*,i] = [ad_stat,ad_crit]
        endif else message, 'NO MODELED DETECTIONS.'
    endfor
    ksv[n] = min(ks[0,*],iks)
    ctf_ksv[n] = ctf[iks]
    adv[n] = min(ad[0,*],iad)
    ctf_adv[n] = ctf[iad]
endfor
print, 'ESTIMATE_CTF - COMPLETE'

sav_vars = ['CTF_KSV','KSV','CTF_ADV','ADV']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_estimate.sav"')


END







