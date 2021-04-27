PRO estimate_free_ctf_update


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
common _fixed


;; run this script NITER times and look at the distribution in CTF
niter = 100
ad_fobv1 = dblarr(niter)
fobv1 = dblarr(niter)
;; binning
step = 0.05d
;; vary obscured fraction
fob = [step:1.-step:step]
nfrac = n_elements(fob)
;; further split obscured sources by NH 23, 24, and 25
fct = [step:1.-step:step]
f23 = 1.-fct
nfixed = n_elements(fct)
;; free CT fraction split between NH=24-25 and 25-26
f25 = [step:1.-step:step]
f24 = 1.-f25
nfree = n_elements(f24)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
print, '=============================================='
print, 'BEGIN - FREE CTF, ROUND 1'
for n = 0,niter-1 do begin
    if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
    ;; resample observed NH distribution to increase data density
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp lt 23.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    ;; NH_RESAMP: structure of increased CT sources
    ;; NH_MOD: draw N sources from NH_RESAMP
    ;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where RX model > RX limits from Carroll+21
    ;; AD: Anderson-Darling test statistic    
    nh_mod = dblarr(nsrc,nfrac,nfixed,nfree)
    rx_mod = dblarr(nsrc,nfrac,nfixed,nfree)
    iimod = bytarr(nsrc,nfrac,nfixed,nfree)
    ad = dblarr(2,nfrac,nfixed,nfree)
    for i = 0,nfrac-1 do begin
        ;; vary obscured fraction
        nob = round((nthin/(1.-fob[i]))*fob[i])
        for j = 0,nfixed-1 do begin
            ;; vary CT fraction
            nct = round(nob*fct[j]);(nob/(1.-fct[j]))*fct[j])
            for k = 0,nfree-1 do begin
                ;; vary NH=24 and NH=25
                n25 = round(nct*f25[k])>1
                n24 = round(nct*f24[k])>1
                n23 = nob-(n24+n25)
                nh_resamp = [nh_samp[ithin],23.+randomu(seed,n23),24.+randomu(seed,n24),25.+randomu(seed,n25)]
                nresamp = n_elements(nh_resamp)
                nh_mod[*,i,j,k] = nh_resamp[randomi(nsrc,nresamp)]
                rx_mod[*,i,j,k] = rx2nh(nh_mod[*,i,j,k],/rx_out,scat=rx_scat)
                iimod[*,i,j,k] = rx_mod[*,i,j,k] gt rxl
                idet = where(iimod[*,i,j,k] eq 1,detct)
                if (detct ge 5) then begin
                    adtwo,rxd,rx_mod[idet,i,j,k],ad_stat,ad_crit
                    ad[*,i,j,k] = [ad_stat,ad_crit]
                endif else if (detct gt 0) then begin
                    ad[*,i,j,k] = 99.
                endif else message, 'NO MODELED DETECTIONS.'
            endfor
        endfor
    endfor
    ad_fobv1[n] = min(ad[0,*,*,*],imin)
    ind = array_indices(ad[0,*,*,*],imin)
    fobv1[n] = fob[ind[1]]
endfor
print, 'END   - FREE CTF, ROUND 1'
print, '=============================================='

sav_vars = ['AD_FOBV1','FOBV1']            
sav_inds = []


;; now rerun this script knowing the obscured fraction
ad_fctv2 = dblarr(niter)
f23v2 = dblarr(niter)
fctv2 = dblarr(niter)
;; fixed obscured fraction
fob = median(fobv1)
print, '=============================================='
print, 'BEGIN - FREE CTF, ROUND 2'
for n = 0,niter-1 do begin
    if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
    ;; resample observed NH distribution to increase data density
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp le 23.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    ;; fix CT fraction
    nob = round((nthin/(1.-fob))*fob)
    ;; NH_RESAMP: structure of increased CT sources
    ;; NH_MOD: draw N sources from NH_RESAMP
    ;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where RX model > RX limits from Carroll+21
    ;; AD: Anderson-Darling test statistic    
    nh_mod = dblarr(nsrc,nfixed,nfree)
    rx_mod = dblarr(nsrc,nfixed,nfree)
    iimod = bytarr(nsrc,nfixed,nfree)
    ad = dblarr(2,nfixed,nfree)
    for j = 0,nfixed-1 do begin
        ;; vary CT fraction
        nct = round(nob*fct[j])
        for k = 0,nfree-1 do begin
            ;; vary NH=24 and NH=25
            n25 = round(nct*f25[k])>1
            n24 = round(nct*f24[k])>1
            n23 = nob-(n24+n25)
            nh_resamp = [nh_samp[ithin],23.+randomu(seed,n23),24.+randomu(seed,n24),25.+randomu(seed,n25)]
            nresamp = n_elements(nh_resamp)
            nh_mod[*,j,k] = nh_resamp[randomi(nsrc,nresamp)]
            rx_mod[*,j,k] = rx2nh(nh_mod[*,j,k],/rx_out,scat=rx_scat)
            iimod[*,j,k] = rx_mod[*,j,k] gt rxl
            idet = where(iimod[*,j,k] eq 1,detct)
            if (detct ge 5) then begin
                adtwo,rxd,rx_mod[idet,j,k],ad_stat,ad_crit
                ad[*,j,k] = [ad_stat,ad_crit]
            endif else if (detct gt 0) then begin
                ad[*,j,k] = 99.
            endif else message, 'NO MODELED DETECTIONS.'
        endfor
    endfor
    ad_fctv2[n] = min(ad[0,*,*],imin)
    ind = array_indices(ad[0,*,*],imin)
    f23v2[n] = total(nh_mod[*,ind[1],ind[2]] ge 23. and nh_mod[*,ind[1],ind[2]] lt 24.)/nsrc
    fctv2[n] = total(nh_mod[*,ind[1],ind[2]] ge 24.)/nsrc
endfor
print, 'END   - FREE CTF, ROUND 2'
print, '=============================================='

sav_vars = [sav_vars,'AD_FCTV2','F23V2','FCTV2']
sav_inds = [sav_inds]


;; and finally, rerun this script knowing both the obscured and CT fraction
ad_fctv3 = dblarr(niter)
f24v3 = dblarr(niter)
f25v3 = dblarr(niter)
;; fixed obscured fraction
f23 = median(f23v2)
fct = median(fctv2)
print, '=============================================='
print, 'BEGIN - FREE CTF, ROUND 3'
for n = 0,niter-1 do begin
    if n mod (ncount/10) eq 0 then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
    ;; resample observed NH distribution to increase data density
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp le 23.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    ;; fix obscured fraction
    tot = round(nthin/(1.-fob))
    nob = round(tot*fob)
    ;; NH_RESAMP: structure of increased CT sources
    ;; NH_MOD: draw N sources from NH_RESAMP
    ;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where RX model > RX limits from Carroll+21
    ;; AD: Anderson-Darling test statistic    
    nh_mod = dblarr(nsrc,nfree)
    rx_mod = dblarr(nsrc,nfree)
    iimod = bytarr(nsrc,nfree)
    ad = dblarr(2,nfree)
    ;; fix CT fraction
    nct = round(tot*fct)
    n23 = nob-nct
    for k = 0,nfree-1 do begin
        ;; vary NH=24 and NH=25
        n25 = round(nct*f25[k])>1
        n24 = round(nct-n25)>1
        nh_resamp = [nh_samp[ithin],23.+randomu(seed,n23),24.+randomu(seed,n24),25.+randomu(seed,n25)]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,k] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,k] = rx2nh(nh_mod[*,k],/rx_out,scat=rx_scat)
        iimod[*,k] = rx_mod[*,k] gt rxl
        idet = where(iimod[*,k] eq 1,detct)
        if (detct ge 5) then begin
            adtwo,rxd,rx_mod[idet,k],ad_stat,ad_crit
            ad[*,k] = [ad_stat,ad_crit]
        endif else if (detct gt 0) then begin
            ad[*,k] = 99.
        endif else message, 'NO MODELED DETECTIONS.'
    endfor
    ad_fctv3[n] = min(ad[0,*],imin)
    f24v3[n] = total(nh_mod[*,imin] ge 24. and nh_mod[*,imin] lt 25.)/nsrc
    f25v3[n] = total(nh_mod[*,imin] ge 25.)/nsrc
endfor
print, 'END   - FREE CTF, ROUND 3'
print, '=============================================='

sav_vars = [sav_vars,'AD_FCTV3','F24V3','F25V3']
sav_inds = [sav_inds]


sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_free.sav"')


END










