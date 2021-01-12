PRO estimate_split


common _data
common _nhobs
common _ctfest

;; STDDEV observed in LX-LMIR relation of Chen+17
rx_scat = 0.2

;; separate WISE AGN, detections and non-detections
iixd = xdet ne '' and iiwac
iixn = xnon ne '' and iiwac
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



niter = 1000
ctf24_ksv = dblarr(nfrac,niter)
ctf25_ksv = dblarr(nfrac,niter)
ctf24_adv = dblarr(nfrac,niter)
ctf25_adv = dblarr(nfrac,niter)

;; empirical distribution function for observed sources
edf,rxd,x_data,edf_data

for n = 0,niter-1 do begin
    if n mod (niter/10) eq 0 then print, strtrim(n/(niter/100),2)+'% complete'

    ;print, 'MODEL_RXDIST FREE - 0% COMPLETE'
    for i = 0,nfrac-1 do begin
        nct = (nthin/(1.-ctf2sig[i]))*ctf2sig[i]
        tag2 = 'split'+string(rnd(ctf24[0]*100,0),format='(i02)')+'_'+string(rnd(ctf25[0]*100,0),format='(i02)')
        re = execute('nh_resamp2 = {'+tag2+':[nh_samp[ithin],24.+randomu(seed,nct*ctf24[0]),25.+randomu(seed,nct*ctf25[0])]}')
        re = execute('nh_mod2 = {'+tag2+':(nh_resamp2.(0))[randomi(nsrc,n_elements(nh_resamp2.(0)))]}')
        re = execute('rx_mod2 = {'+tag2+':rx2nh_z(nh_mod2.(0),z[iwagn],/rx_out,scat=rx_scat)}')
        re = execute('iimod2 = {'+tag2+':rx_mod2.(0) gt rxl}')
        ks2 = dblarr(2,nsplit)
        kstwo,rxd,(rx_mod2.(0))[where(iimod2.(0) eq 1)],d,prob
        ks2[*,0] = [d,prob]
        ad = dblarr(nsplit)
        edf,(rx_mod2.(0))[where(iimod2.(0) eq 1)],x_model,edf_model
        cdf_model = interpol(edf_model,x_model,x_data)
        ad_stat = (edf_data-cdf_model)^2./(cdf_model*(1.-cdf_model))
        ad[0] = total(ad_stat,/nan)/total(finite(ad_stat))
        for j = 1,nsplit-1 do begin
            tag2 = 'split'+string(rnd(ctf24[j]*100,0),format='(i02)')+'_'+string(rnd(ctf25[j]*100,0),format='(i02)')
            ;nct = (nthin/(1.-ctf[i]))*ctf[i]
            nh_resamp2 = create_struct(nh_resamp2,tag2,[nh_samp[ithin],24.+randomu(seed,nct*ctf24[j]),25.+randomu(seed,nct*ctf25[j])])
            nh_mod2 = create_struct(nh_mod2,tag2,(nh_resamp2.(j))[randomi(nsrc,n_elements(nh_resamp2.(j)))])
            rx_mod2 = create_struct(rx_mod2,tag2,rx2nh_z(nh_mod2.(j),z[iwagn],/rx_out,scat=rx_scat))
            iimod2 = create_struct(iimod2,tag2,rx_mod2.(j) gt rxl)
            kstwo,rxd,(rx_mod2.(j))[where(iimod2.(j) eq 1)],d,prob
            ks2[*,j] = [d,prob]
            edf,(rx_mod2.(j))[where(iimod2.(j) eq 1)],x_model,edf_model
            edf_model = interpol(edf_model,x_model,x_data)
            ad_stat = (edf_data-edf_model)^2./(edf_model*(1-edf_model))
            ad[j] = total(ad_stat,/nan)/total(finite(ad_stat))
        endfor
        ;re = execute("frac"+string(rnd(ctf[i]*100,0),format="(i02)")+' = {nh_resamp:nh_resamp2,nh_mod:nh_mod2,rx_mod:rx_mod2,iimod:iimod2,ks:ks2}')
        ;if (i gt 0 and i mod (nfrac/2.) eq 0) then print, 'MODEL_RXDIST FREE - 50% COMPLETE'
        ksm = min(ks2[0,*],imin)
        ctf24_ksv[i,n] = ctf24[imin]
        ctf25_ksv[i,n] = ctf25[imin]
        adm = min(ad,imin)
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