PRO estimate_ctf


common _data
common _nhobs
common _rxnh

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
    re = execute('rx_mod = {'+tag+':rx2nh(nh_mod.(0),/rx_out,scat=rx_scat)}')
    re = execute('iimod = {'+tag+':rx_mod.(0) gt rxl}')
    ks = dblarr(2,nfrac)
    kstwo,rxd,(rx_mod.(0))[where(iimod.(0) eq 1)],ks_stat,ks_prob
    ks[*,0] = [ks_stat,ks_prob]
    ad = dblarr(2,nfrac)
    adtwo,rxd,(rx_mod.(0))[where(iimod.(0) eq 1)],ad_stat,ad_crit
    ad[*,0] = [ad_stat,ad_crit]
    for i = 1,nfrac-1 do begin
        tag = 'frac'+string(rnd(ctf[i]*100,0),format='(i02)')
        nct = (nthin/(1.-ctf[i]))*ctf[i]
        nh_resamp = create_struct(nh_resamp,tag,[nh_samp[ithin],24.+randomu(seed,nct)*nh_diff])
        nh_mod = create_struct(nh_mod,tag,(nh_resamp.(i))[randomi(nsrc,n_elements(nh_resamp.(i)))])
        rx_mod = create_struct(rx_mod,tag,rx2nh(nh_mod.(i),/rx_out,scat=rx_scat))
        iimod = create_struct(iimod,tag,rx_mod.(i) gt rxl)
        kstwo,rxd,(rx_mod.(i))[where(iimod.(i) eq 1)],ks_stat,ks_prob
        ks[*,i] = [ks_stat,ks_prob]
        adtwo,rxd,(rx_mod.(i))[where(iimod.(i) eq 1)],ad_stat,ad_crit
        ad[*,i] = [ad_stat,ad_crit]
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







