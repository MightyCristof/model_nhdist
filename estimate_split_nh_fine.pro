PRO estimate_split_nh


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
;common _fixed
common _free


cvm = 0

;; run this script NITER times 
niter = 1000;0
f24v2 = dblarr(niter)
f25v2 = dblarr(niter)
statv2 = dblarr(6,niter)
nhmv2 = dblarr(nsrc,niter)
rxmv2 = dblarr(nsrc,niter)
iimv2 = bytarr(nsrc,niter)

;; free CT fraction split between NH=24-25 and 25-26
step = 0.05d
f25 = [0.05d:1.-0.05d:0.05d]
f24 = 1.-f25
nfree = n_elements(f24)

a2_fine2 = dblarr(nfree,niter)

;; set CT fraction from CTF FREE modeling
fct = mode(fctv1_fine,kde=kde_bandwidth(fctv1_fine))
;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
nrej2 = 0l
for n = 0,niter-1 do begin
    ;; resample observed NH distribution to increase data density
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp le 24.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    ;; fix CT fraction
    nct = round((nthin/(1.-fct))*fct)
    ;; NH_RESAMP: structure of increased CT sources
    ;; NH_MOD: draw N sources from NH_RESAMP
    ;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where RX model > RX limits from Carroll+21
    ;; AD: Anderson-Darling test statistic    
    nh_mod = dblarr(nsrc,nfree)
    rx_mod = dblarr(nsrc,nfree)
    iimod = bytarr(nsrc,nfree)
    a2 = dblarr(nfree)
    p_a2 = dblarr(nfree)
    rxd_scat = randomn(seed,ndet)*e_rxd
    ;rxm_scat = randomn(seed,nsrc)*rx_scat
    for j = 0,nfree-1 do begin
        n25 = round(nct*f25[j])>1
        n24 = nct-n25
        nh_resamp = [nh_samp[ithin],24.+randomu(seed,n24),25.+randomu(seed,n25)]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,j] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,j] = rx2nh(nh_mod[*,j],/rx_out,scat=rx_scat)
        iimod[*,j] = rx_mod[*,j] gt rxl
        idet = where(iimod[*,j] eq 1,detct)
        if (detct ge 5) then begin
            ;; if comparing test statistics, need p_a2
            a2[j] = ad_test(rxd+rxd_scat,rx_mod[idet,j],permute=(test eq 'JOINT'),prob=p,cvm=cvm)
            p_a2[j] = p
            ;kstwo,rxd,rx_mod[idet,j],d,prob
            ;a2[j] = d
            ;p_a2[j] = prob            
        endif else if (detct gt 0) then begin
            a2[j] = -1.
            p_a2[j] = -1.
        endif else message, 'NO MODELED DETECTIONS.'
    endfor    
    iia2 = p_a2 gt min(p_a2)
    a2_fine2[*,n] = a2
    
    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case test of 
        'AD': begin
            ibest = where(a2 eq min(a2[where(a2 gt 0,/null)]),nbest)
            ;if (nbest gt 1) then stop
            if (nbest gt 1) then ibest = ibest[0]
            stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
            end
        'JOINT': begin
            ;; constraints by X-ray stacked fluxes
            fx_est = estimate_fx(rx_mod,iimod,/iterate)
            ;; compare to X-ray stacked fluxes
            ;x2_soft = ((fxstak[1,0]-fx_est.soft)/e_fxstak[1,0])^2.
            ;x2_hard = ((fxstak[1,1]-fx_est.hard)/e_fxstak[1,1])^2.
            ;; uncertainties in quadrature
            unc_soft = abs(fxstak[1,0] * sqrt((e_fxstak[1,0]/fxstak[1,0])^2. + (fx_est.e_soft/fx_est.soft)^2.))
            unc_hard = abs(fxstak[1,1] * sqrt((e_fxstak[1,1]/fxstak[1,1])^2. + (fx_est.e_hard/fx_est.hard)^2.))
            x2_soft = ((fxstak[1,0]-fx_est.soft)/unc_soft)^2.
            x2_hard = ((fxstak[1,1]-fx_est.hard)/unc_hard)^2.
            x2 = x2_soft+x2_hard
            p_x2 = 1.-chisqr_pdf(x2,1.) ;; dof = 1 (2 X-ray data points - 1)
            ;; correct for p-value == 1
            p_x2 = p_x2 > min(p_x2[where(p_x2,/null)])/2.
            iix2 = p_x2 gt min(p_x2)
            ;; combined test statistic
            x2_joint = -2.*(alog(p_a2)+alog(p_x2))<99.
            p_joint = 1.-chisqr_pdf(x2_joint,4.) ;; dof = 2k, where k is the number of tests being combined            
            ibest = where(x2_joint eq min(x2_joint[where(iia2 and iix2,/null)]),nbest)
            if (nbest gt 1) then stop
            stat = [a2[ibest],p_a2[ibest],x2[ibest],p_x2[ibest],x2_joint[ibest],p_joint[ibest]]
            end
        else: message, 'IMPROPER INPUT OF TEST METHOD.'
    endcase
    ;; skip if no good fits during iteration
    if (nbest eq 0) then begin
        nrej2++
        continue
    ;; record best fit statistics
    endif else begin
        f24v2[n] = f24[ibest]
        f25v2[n] = f25[ibest]
        statv2[*,n] = stat
        nhmv2[*,n] = nh_mod[*,ibest]
        rxmv2[*,n] = rx_mod[*,ibest]
        iimv2[*,n] = iimod[*,ibest]
    endelse
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - NH SPLIT'
    endif else if (n mod (ncount/10) eq 0) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - NH SPLIT'
print, '=============================================='

sav_vars = ['F24V2','F25V2','STATV2','NREJ2','NHMV2','RXMV2','IIMV2', $
            'A2_FINE2']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="nh_split.sav"')


END






