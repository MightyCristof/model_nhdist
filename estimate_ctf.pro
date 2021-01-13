PRO estimate_ctf


common _data
common _nhobs


;; STDDEV observed in LX-LMIR relation of Chen+17
rx_scat = 0.2

;; separate WISE AGN, detections and non-detections
iixd = xdet ne '' and iiwac
iixn = xnon ne '' and iiwac
rxd = rldet[where(iixd,ndet)]
e_rxd = e_rldet[where(iixd)]
iwagn = where(iiwac,nsrc)
rxl = rxlim[iwagn]


;; run this script NITER times and look at the distribution in CTF
niter = 1000
ctf_ksv = dblarr(niter)
ctf_adv = dblarr(niter)

;; empirical distribution function for observed sources
;edf,rxd,x_data,edf_data

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
    div = 50.
    ctf = dindgen(div-1,start=1)/div
    nfrac = n_elements(ctf)
    ;; NH_RESAMP: structure of increased CT sources. scale1_0 == nh_samp
    ;; NH_RESAMP has different number of sources per tag, cannot use SOA2AOS
    ;; NH_MOD: draw (from NH_RESAMP) a number of sources == X-ray non-detections
    ;; RX_MOD: convert the NH model to RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where rx_mod > RX non-detection from Carroll+20
    ;; KS: KS test between RX model detections and RX detections from Carroll+20
    ;; AD: Anderson-Darling test statistic
    nct = (nthin/(1.-ctf[0]))*ctf[0]
    tag = 'frac'+string(rnd(ctf[0]*100,0),format='(i02)')
    re = execute('nh_resamp = {'+tag+':[nh_samp[ithin],24.+randomu(seed,nct)*nh_diff]}')
    re = execute('nh_mod = {'+tag+':(nh_resamp.(0))[randomi(nsrc,n_elements(nh_resamp.(0)))]}')
    re = execute('rx_mod = {'+tag+':rx2nh_z(nh_mod.(0),z[iwagn],/rx_out,scat=rx_scat)}')
    re = execute('iimod = {'+tag+':rx_mod.(0) gt rxl}')
    ks = dblarr(2,nfrac)
    kstwo,rxd,(rx_mod.(0))[where(iimod.(0) eq 1)],ks_stat,ks_prob
    ks[*,0] = [ks_stat,ks_prob]
    ad = dblarr(2,nfrac)
    adtwo,rxd,(rx_mod.(0))[where(iimod.(0) eq 1)],ad_stat,ad_crit
    ad[*,0] = [ad_stat,ad_crit]
    ;edf,(rx_mod.(0))[where(iimod.(0) eq 1)],x_model,edf_model
    ;edf_model = interpol(edf_model,x_model,x_data)
    ;adv = (edf_data-edf_model)^2./(edf_model*(1-edf_model))
    ;ad[0] = total(adv,/nan)/total(finite(adv))
    ;; repeat for all scalings
    for i = 1,nfrac-1 do begin
        tag = 'frac'+string(rnd(ctf[i]*100,0),format='(i02)')
        nct = (nthin/(1.-ctf[i]))*ctf[i]
        nh_resamp = create_struct(nh_resamp,tag,[nh_samp[ithin],24.+randomu(seed,nct)*nh_diff])
        nh_mod = create_struct(nh_mod,tag,(nh_resamp.(i))[randomi(nsrc,n_elements(nh_resamp.(i)))])
        rx_mod = create_struct(rx_mod,tag,rx2nh_z(nh_mod.(i),z[iwagn],/rx_out,scat=rx_scat))
        iimod = create_struct(iimod,tag,rx_mod.(i) gt rxl)
        kstwo,rxd,(rx_mod.(i))[where(iimod.(i) eq 1)],ks_stat,ks_prob
        ks[*,i] = [ks_stat,ks_prob]
        adtwo,rxd,(rx_mod.(i))[where(iimod.(i) eq 1)],ad_stat,ad_crit
        ad[*,i] = [ad_stat,ad_crit]
        ;edf,(rx_mod.(i))[where(iimod.(i) eq 1)],x_model,edf_model
        ;edf_model = interpol(edf_model,x_model,x_data)
        ;adv = (edf_data-edf_model)^2./(edf_model*(1-edf_model))
        ;ad[i] = total(adv,/nan)/total(finite(adv))
    endfor
    ks_min = min(ks[0,*],iks)
    ctf_ksv[n] = ctf[iks]
    ad_min = min(ad[0,*],iad)
    ctf_adv[n] = ctf[iad]
endfor
print, 'ESTIMATE_CTF - COMPLETE'

sav_vars = ['CTF_KSV','CTF_ADV']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_estimate.sav"')


END









;    ;; now do the same as above, but allow the range of NH to differ in bins
;    ;; CTF must sum to 1 across bins
;    div2 = 20.
;    ctf25 = dindgen(div2-1,start=1)/div2
;    ctf24 = 1.-ctf25
;    nsplit = n_elements(ctf24)
;    for i = 0,nfrac-1 do begin
;        nct = (nthin/(1.-ctf[i]))*ctf[i]
;        tag2 = 'split'+string(rnd(ctf24[0]*100,0),format='(i02)')+'_'+string(rnd(ctf25[0]*100,0),format='(i02)')
;        re = execute('nh_resamp2 = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct*ctf24[0]),25.+randomu(seed,nct*ctf25[0])]}')
;        re = execute('nh_mod2 = {'+tag2+':(nh_resamp2.(0))[randomi(nsrc,n_elements(nh_resamp2.(0)))]}')
;        re = execute('rx_mod2 = {'+tag2+':rl2nh(nh_mod2.(0),model="BORUS",/rl_out,scat=rx_scat)}')
;        re = execute('iimod2 = {'+tag2+':rx_mod2.(0) gt rxl}')
;        ks2 = dblarr(2,nsplit)
;        kstwo,rxd,(rx_mod2.(0))[where(iimod2.(0) eq 1)],d,prob
;        ks2[*,0] = [d,prob]
;        for j = 1,nsplit-1 do begin
;            tag2 = 'split'+string(rnd(ctf24[j]*100,0),format='(i02)')+'_'+string(rnd(ctf25[j]*100,0),format='(i02)')
;            ;nct = (nthin/(1.-ctf[i]))*ctf[i]
;            nh_resamp2 = create_struct(nh_resamp2,tag2,[nh_samp[ithin],24.+randomu(seed,nct*ctf24[j]),25.+randomu(seed,nct*ctf25[j])])
;            nh_mod2 = create_struct(nh_mod2,tag2,(nh_resamp2.(j))[randomi(nsrc,n_elements(nh_resamp2.(j)))])
;            rx_mod2 = create_struct(rx_mod2,tag2,rl2nh(nh_mod2.(j),model='BORUS',/rl_out,scat=rx_scat))
;            iimod2 = create_struct(iimod2,tag2,rx_mod2.(j) gt rxl)
;            kstwo,rxd,(rx_mod2.(j))[where(iimod2.(j) eq 1)],d,prob
;            ks2[*,j] = [d,prob]
;        endfor
;        re = execute("frac"+string(rnd(ctf[i]*100,0),format="(i02)")+' = {nh_resamp:nh_resamp2,nh_mod:nh_mod2,rx_mod:rx_mod2,iimod:iimod2,ks:ks2}')
;        ;if (ctf[i]*100 mod 10 eq 0) then print, 'MODEL_RXDIST FREE - COMPLETING FRAC: '+string(ctf[i],format='(f0.1)')
;    endfor
;    ;print, 'MODEL_RXDIST FREE - COMPLETE'
;
;    frac_str = 'frac'+string(rnd(ctf*100,0),format='(i02)')
;    frac_str = frac_str+':'+frac_str+'.'
;    re = execute('nh_resamp2 = {'+strjoin(frac_str+"nh_resamp",",")+'}')
;    re = execute('nh_mod2 = {'+strjoin(frac_str+"nh_mod",",")+'}')
;    re = execute('rx_mod2 = {'+strjoin(frac_str+"rx_mod",",")+'}')
;    re = execute('iimod2 = {'+strjoin(frac_str+"iimod",",")+'}')
;    re = execute('ks2 = {'+strjoin(frac_str+"ks",",")+'}')
;
;
;    ks2m = dblarr(nfrac)
;    iks2m = intarr(nfrac)
;    for i = 0,nfrac-1 do begin
;        ks2m[i] = min(ks2.(i)[0,*],imin)
;        iks2m[i] = imin
;    endfor
;
;    !null = min(ks2m,imin2)
;    iks2 = [imin2,iks2m[imin2]]
;
;    ctf_ksv[1:2,n] = iks2
;
;endfor
;;sav_vars = [sav_vars,'CTF24','CTF25','NSPLIT', $
;;                     'NH_RESAMP2','NH_MOD2','RX_MOD2','KS2']
;;sav_inds = [sav_inds,'IIMOD2','IKS2']
;;
;;
;;sav_str = strjoin([sav_vars,sav_inds],',')
;;re = execute('save,'+sav_str+',/compress,file="rx_model.sav"')
;save,ctf_ksv,file='ctf_ksv.sav'
;