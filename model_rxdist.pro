PRO model_rxdist, PLT = plt


common _data
common _nhobs
common _ctfest
common _split

;; STDDEV observed in LX-LMIR relation of Chen+17
rx_scat = 0.2

;; separate WISE AGN, detections and non-detections
iixd = xdet ne '' and iiwac
iixn = xnon ne '' and iiwac
rxd = rldet[where(iixd,ndet)]
e_rxd = e_rldet[where(iixd)]
rxl = rxlim[where(iiwac,nsrc)]

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
rx_mod_ks = rl2nh(nh_mod_ks,model="BORUS",/rl_out,scat=rx_scat)
iimod_ks = rx_mod_ks gt rxl
kstwo,rxd,rx_mod_ks[where(iimod_ks eq 1)],d,prob
ks = [d,prob]
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
;nct_ks = (nthin/(1.-ctf))*ctf
;tag2 = 'split'+string(rnd(ctf24[0]*100,0),format='(i02)')+'_'+string(rnd(ctf25[0]*100,0),format='(i02)')
re = execute('tag2 = "iter"+string(0,format="('+fmt2+')")')
re = execute('nh_resamp2_ks = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct_ks*ctf24_ks[0]),25.+randomu(seed,nct_ks*ctf25_ks[0])]}')
re = execute('nh_mod2_ks = {'+tag2+':(nh_resamp2_ks.(0))[randomi(nsrc,n_elements(nh_resamp2_ks.(0)))]}')
re = execute('rx_mod2_ks = {'+tag2+':rl2nh(nh_mod2_ks.(0),model="BORUS",/rl_out,scat=rx_scat)}')
re = execute('iimod2_ks = {'+tag2+':rx_mod2_ks.(0) gt rxl}')
ks2 = dblarr(2,niter)
kstwo,rxd,(rx_mod2_ks.(0))[where(iimod2_ks.(0) eq 1)],d,prob
ks2[*,0] = [d,prob]
for j = 1,niter-1 do begin
    re = execute('tag2 = "iter"+string(j,format="('+fmt2+')")')
    nh_resamp2_ks = create_struct(nh_resamp2_ks,tag2,[nh_samp[ithin],24.+randomu(seed,nct_ks*ctf24_ks),25.+randomu(seed,nct_ks*ctf25_ks)])
    nh_mod2_ks = create_struct(nh_mod2_ks,tag2,(nh_resamp2_ks.(j))[randomi(nsrc,n_elements(nh_resamp2_ks.(j)))])
    rx_mod2_ks = create_struct(rx_mod2_ks,tag2,rl2nh(nh_mod2_ks.(j),model='BORUS',/rl_out,scat=rx_scat))
    iimod2_ks = create_struct(iimod2_ks,tag2,rx_mod2_ks.(j) gt rxl)
    kstwo,rxd,(rx_mod2_ks.(j))[where(iimod2_ks.(j) eq 1)],d,prob
    ks2[*,j] = [d,prob]
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

;; empirical distribution function for observed sources
edf,rxd,x_data,edf_data

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
rx_mod_ad = rl2nh(nh_mod_ad,model="BORUS",/rl_out,scat=rx_scat)
iimod_ad = rx_mod_ad gt rxl
edf,rx_mod_ad[where(iimod_ad eq 1)],x_model,edf_model
edf_model = interpol(edf_model,x_model,x_data)
adv = (edf_data-edf_model)^2./(edf_model*(1-edf_model))
ad = total(adv,/nan)/total(finite(adv))

sav_vars = [sav_vars,'CTF_AD','NH_RESAMP_AD','NH_MOD_AD','RX_MOD_AD','KS']
sav_inds = [sav_inds,'IIMOD_AD','IKS']


;; now do the same as above, but allow the range of NH to differ in bins
;; CTF must sum to 1 across bins
ctf24_ad = median(ctf24_adv[ifrac_ad,*])
ctf25_ad = 1.-ctf24_ad;rnd(median(ctf25_adv[ifrac_ad,*]),2)

niter = 1000
fmt2 = 'i0'+strtrim(strlen(strtrim(niter,2)),2)

print, 'AD: MODEL_RXDIST FREE - 0% COMPLETE'
;nct_ad = (nthin/(1.-ctf))*ctf
;tag2 = 'split'+string(rnd(ctf24[0]*100,0),format='(i02)')+'_'+string(rnd(ctf25[0]*100,0),format='(i02)')
re = execute('tag2 = "iter"+string(0,format="('+fmt2+')")')
re = execute('nh_resamp2_ad = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct_ad*ctf24_ad[0]),25.+randomu(seed,nct_ad*ctf25_ad[0])]}')
re = execute('nh_mod2_ad = {'+tag2+':(nh_resamp2_ad.(0))[randomi(nsrc,n_elements(nh_resamp2_ad.(0)))]}')
re = execute('rx_mod2_ad = {'+tag2+':rl2nh(nh_mod2_ad.(0),model="BORUS",/rl_out,scat=rx_scat)}')
re = execute('iimod2_ad = {'+tag2+':rx_mod2_ad.(0) gt rxl}')
ad2 = dblarr(niter)
edf,(rx_mod2_ad.(0))[where(iimod2_ad.(0) eq 1)],x_model,edf_model
edf_model = interpol(edf_model,x_model,x_data)
adv2 = (edf_data-edf_model)^2./(edf_model*(1-edf_model))
ad2[0] = total(adv2,/nan)/total(finite(adv2))
for j = 1,niter-1 do begin
    re = execute('tag2 = "iter"+string(j,format="('+fmt2+')")')
    nh_resamp2_ad = create_struct(nh_resamp2_ad,tag2,[nh_samp[ithin],24.+randomu(seed,nct_ad*ctf24_ad),25.+randomu(seed,nct_ad*ctf25_ad)])
    nh_mod2_ad = create_struct(nh_mod2_ad,tag2,(nh_resamp2_ad.(j))[randomi(nsrc,n_elements(nh_resamp2_ad.(j)))])
    rx_mod2_ad = create_struct(rx_mod2_ad,tag2,rl2nh(nh_mod2_ad.(j),model='BORUS',/rl_out,scat=rx_scat))
    iimod2_ad = create_struct(iimod2_ad,tag2,rx_mod2_ad.(j) gt rxl)
    edf,(rx_mod2_ad.(j))[where(iimod2_ad.(j) eq 1)],x_model,edf_model
    edf_model = interpol(edf_model,x_model,x_data)
    adv2 = (edf_data-edf_model)^2./(edf_model*(1-edf_model))
    ad2[j] = total(adv2,/nan)/total(finite(adv2))
    if (j mod (niter/2.) eq 0) then print, 'AD: MODEL_RXDIST FREE - 50% COMPLETE'
endfor
print, 'AD: MODEL_RXDIST FREE - 100% COMPLETE'

ad2med = median(ad2)
isort = sort(ad2)
iloc = value_locate(ad2[isort],ad2med)
iad2 = isort[iloc]


sav_vars = [sav_vars,'CTF24_AD','CTF25_AD', $
                     'NH_RESAMP2_AD','NH_MOD2_AD','RX_MOD2_AD','AD2']
sav_inds = [sav_inds,'IIMOD2_AD','IAD2']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="rx_model2.sav"')


END












;if keyword_set(add) then begin
;    ;; combine RX for detections and non-detections
;    rl = rldet>rlnon
;    e_rl = e_rldet>e_rlnon
;    ;; split WISE AGN + secondary
;    rlw = rl[where(iiwac)]
;    e_rlw = e_rl[where(iiwac)]
;    rls = rl[where(iiwac eq 0)]
;    e_rls = e_rl[where(iiwac eq 0)]
;
;    ;; here we draw up a distribution in NH
;    nw = n_elements(rlw)
;    ns = n_elements(rls)
;    npull = nw
;    rl_obs = rlw
;    nh_obs = nh_mc(nh_lan_nst,npull)
;
;    ;; here we create random distributions to add to the 
;    peak = 22.+dindgen(8)/2.
;    npeak = n_elements(peak)
;    samp = round([reform(findgen(9,start=1)#10.^findgen(4,start=-2),9*4),1e2] * npull)
;    nsamp = n_elements(samp)
;    nh_add = dblarr(max(samp),npeak,nsamp)
;    for i = 0,npeak-1 do $
;        for j = 0,nsamp-1 do $
;            nh_add[0:samp[j]-1,i,j] = peak[i] + randomn(seed,samp[j])*0.5
;
;    ;; plot for results
;    ;plothist, nh_add[*,-1,-1],xra=[20,28]
;    ;for i = 0,npeak-1 do begin
;    ;    for j = 0,nsamp-1 do begin
;    ;        ig = where(nh_add[*,i,j] gt 0)
;    ;        plothist, nh_add[ig,i,j],/ov
;    ;    endfor
;    ;endfor
;
;    ;; here we convert to RX and test the CDF
;    ksd = dblarr(npeak,nsamp)
;    ksp = dblarr(npeak,nsamp)
;    ;; STDDEV = 0.3
;    ;; MAD = 0.2
;    for i = 0,npeak-1 do begin
;        for j = 0,nsamp-1 do begin
;            nh_mod = [nh_obs,nh_add[0:samp[j]-1,i,j]]
;            rx_mod = rl2nh(nh_mod,model='BORUS',/rl_out) + randomn(seed,n_elements(nh_mod))*rx_scat
;            kstwo,rl_obs,rx_mod,d,prob
;            ksd[i,j] = d
;            ksp[i,j] = prob
;        endfor
;    endfor
;
;    dmin = min(ksd,imin)
;    ind = array_indices(ksd,imin)
;
;    nh_mod = [nh_obs,nh_add[0:samp[ind[1]]-1,ind[0],ind[1]]]
;    rx_mod = rl2nh(nh_mod,model='BORUS',/rl_out)
;
;    ;; plot comparison
;    ;col = [[80,193,128],[83,5,124]]
;    col = ['teal','purple']
;    bin = 0.1
;    plothist,rl_obs,/peak,bin=bin,col=col[0],xra=[-3.,1.],yra=[0.,1.2]
;    plothist,rx_mod,/peak,bin=bin,col=col[1],/ov
;
;    yobs = cdf(rl_obs,xloc=xobs)
;    ymod = cdf(rx_mod,xloc=xmod)
;    dx = [xobs,xmod]
;    dx = dx[uniq(dx,sort(dx))]
;    dyobs = interpol(yobs,xobs,dx)
;    dymod = interpol(ymod,xmod,dx)
;
;    pobs = plot(dx,dyobs,'-',thick=2,col=col[0],yra=[-0.05,1.05],name='Carroll+20')
;    pmod = plot(dx,dymod,'--',thick=2,col=col[1],/ov,name='Model')
;    l = legend(target=[pobs,pmod],position=[0.32,0.92],/relative)
;endif
;