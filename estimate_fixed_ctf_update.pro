PRO estimate_fixed_ctf_update


common _data
common _nhdist
common _nhobs
common _rxnh
common _group


;; run this script NITER times and look at the distribution in CTF
niter = 100
ad_fobv = dblarr(niter)
fobv = dblarr(niter)
;; step size
step = 0.01d
;; vary obscured fraction (NH>23)
fob = [step:1.-step:step]
nfrac = n_elements(fob)
;; further split obscured sources by 23<NH<24 and NH>24
fct = [step:1.-step:step]
f23 = 1.-fct
nfixed = n_elements(fct)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
print, '=============================================='
print, 'BEGIN - FIXED CTF, ROUND 1'
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
    nh_mod = dblarr(nsrc,nfrac,nfixed)
    rx_mod = dblarr(nsrc,nfrac,nfixed)
    iimod = bytarr(nsrc,nfrac,nfixed)
    ad = dblarr(2,nfrac,nfixed)
    for i = 0,nfrac-1 do begin
        nob = round((nthin/(1.-fob[i]))*fob[i])
        for j = 0,nfixed-1 do begin
            nct = round(nob*fct[j])>1
            n23 = nob-nct
            nh_resamp = [nh_samp[ithin],23.+randomu(seed,n23),24.+2.*randomu(seed,nct)]
            nresamp = n_elements(nh_resamp)
            nh_mod[*,i,j] = nh_resamp[randomi(nsrc,nresamp)]
            rx_mod[*,i,j] = rx2nh(nh_mod[*,i,j],/rx_out,scat=rx_scat)
            iimod[*,i,j] = rx_mod[*,i,j] gt rxl
            idet = where(iimod[*,i,j] eq 1,detct)
            if (detct ge 5) then begin
                adtwo,rxd,rx_mod[idet,i,j],ad_stat,ad_crit
                ad[*,i,j] = [ad_stat,ad_crit]
            endif else if (detct gt 0) then begin
                ad[*,i,j] = 99.
            endif else message, 'NO MODELED DETECTIONS.'
        endfor
    endfor
    ad_fobv[n] = min(ad[0,*,*],imin)
    ind = array_indices(ad[0,*,*],imin)
    fobv[n] = fob[ind[1]]
endfor
print, 'END   - FIXED CTF, ROUND 1'
print, '=============================================='

sav_vars = ['AD_FOBV','FOBV']
sav_inds = []

;; now rerun this script knowing the CTF for the free modeling to find f23/fct
ad_fctv = dblarr(niter)
f23v = dblarr(niter)
fctv = dblarr(niter)
;; fix obscured fraction
fob = median(fobv)
print, '=============================================='
print, 'BEGIN - FIXED CTF, ROUND 2'
for n = 0,niter-1 do begin
    if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
    ;; resample observed NH distribution to increase data density
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp lt 23.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    ;; fix CT fraction
    nob = round((nthin/(1.-fob))*fob)
    ;; NH_RESAMP: structure of increased CT sources
    ;; NH_MOD: draw N sources from NH_RESAMP
    ;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where RX model > RX limits from Carroll+21
    ;; AD: Anderson-Darling test statistic    
    nh_mod = dblarr(nsrc,nfixed)
    rx_mod = dblarr(nsrc,nfixed)
    iimod = bytarr(nsrc,nfixed)
    ad = dblarr(2,nfixed)
    for j = 0,nfixed-1 do begin
        nct = round(nob*fct[j])
        n23 = nob-nct
        nh_resamp = [nh_samp[ithin],23.+randomu(seed,n23),24.+2.*randomu(seed,nct)]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,j] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,j] = rx2nh(nh_mod[*,j],/rx_out,scat=rx_scat)
        iimod[*,j] = rx_mod[*,j] gt rxl
        idet = where(iimod[*,j] eq 1,detct)
        if (detct ge 5) then begin
            adtwo,rxd,rx_mod[idet,j],ad_stat,ad_crit
            ad[*,j] = [ad_stat,ad_crit]
        endif else if (detct gt 0) then begin
            ad[*,j] = 99.
        endif else message, 'NO MODELED DETECTIONS.'
    endfor
    ad_fctv[n] = min(ad[0,*],imin)
    f23v[n] = total(nh_mod[*,imin] ge 23. and nh_mod[*,imin] lt 24.)/nsrc
    fctv[n] = total(nh_mod[*,imin] ge 24.)/nsrc
endfor
print, 'END   - FIXED CTF, ROUND 2'
print, '=============================================='

sav_vars = [sav_vars,'AD_FCTV','F23V','FCTV']
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_fixed.sav"')


END







