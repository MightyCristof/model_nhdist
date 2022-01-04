PRO ctf_variable, CHISQ = chisq


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
common _uniform


;; AD vs CvM test
cvm = 0

;; run this script NITER times 
niter = 1000

;; modeling variables
f24v_ = dblarr(niter)
f25v_ = dblarr(niter)
statv_ = dblarr(6,niter)
nhmv_ = dblarr(nsrc,niter)
rxmv_ = dblarr(nsrc,niter)
iimv_ = bytarr(nsrc,niter)
;; data variables
rxdv_ = dblarr(ndet,niter)

;; variable CT fraction split between NH=24-25 and 25-26
step_ = 0.05d
f25_steps = [step_:1.-step_:step_]
f24_steps = 1.-f25_steps
nfree = n_elements(f24_steps)

a2v_ = dblarr(nfree,niter)

;; set CT fraction from CTF FREE modeling
fct = mode(fctv,kde=kde_bandwidth(fctv))
;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
for n = 0,niter-1 do begin
    ;; resample observed NH distribution to increase data density
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp le 24.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    ;; fix CT fraction
    nct = round((nthin/(1.-fct))*fct)
    ;; iterative modeling arrays
    nh_mod = dblarr(nsrc,nfree)
    rx_mod = dblarr(nsrc,nfree)
    iimod = bytarr(nsrc,nfree)
    mdetf = dblarr(nfree)
    a2 = dblarr(nfree)
    p_a2 = dblarr(nfree)
    rxdv_[*,n] = rxd+randomn(seed,ndet)*0.23;rx_scat
    ;rxdv_[*,n] = rxd+randomn(seed,ndet)*e_rxd
    for j = 0,nfree-1 do begin
        n25 = round(nct*f25_steps[j])>1
        n24 = nct-n25
        nh_resamp = [nh_samp[ithin],24.+randomu(seed,n24),25.+randomu(seed,n25)]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,j] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,j] = rx2nh(nh_mod[*,j],/rx_out,scat=rx_scat)
        iimod[*,j] = rx_mod[*,j] gt rxl
        idet = where(iimod[*,j] eq 1,moddet,ncomplement=modnon)
        rxmn = dblarr(modnon)-9999.
        if keyword_set(chisq) then begin
            a2[j] = 1.
            p_a2[j] = 1.
        endif else begin
            if (moddet ge 5) then begin
                a2[j] = ad_test(rxdv_[*,n],rx_mod[idet,j],permute=1,prob=p,cvm=cvm)
                a2[j] = ad_test([rxdn,rxdv_[*,n]],[rxmn,rx_mod[idet,j]],permute=0,prob=p,weight=1)
                p_a2[j] = p
            endif else if (moddet gt 0) then begin
                a2[j] = -1.
                p_a2[j] = -1.
            endif else message, 'NO MODELED DETECTIONS.'
        endelse
    endfor
    ;; optional penalty
    if keyword_set(penalty) then begin
        ;; penalize large disparity in NH split
        spenalty = (f24_steps-f25_steps)^2.
        sp = spenalty/nfree
        pp = a2/(a2+sp)
        a2 += sp
        p_a2 *= pp
    endif
    ;; finite values only
    iia2 = finite(p_a2) and p_a2 gt 0.
    a2v_[*,n] = a2
    ;; determine "best-fit" model per iteration
    ;; constraints by X-ray stacked fluxes
    fx_est = estimate_fx(rx_mod,iimod,/iterate)
    ;; X-ray stack data minus model
    del_soft = abs(fxstak[1,0]-fx_est.csoft)
    del_hard = abs(fxstak[1,1]-fx_est.chard)
    ;; uncertainties
    sig_soft = e_fxstak[1,0]
    sig_hard = e_fxstak[1,1]
    ;; chi-square
    x2_soft = (del_soft/sig_soft)^2.
    x2_hard = (del_hard/sig_hard)^2.    
    x2 = x2_soft+x2_hard
    p_x2 = 1.-chisqr_pdf(x2,1.) ;; dof = 1 (2 X-ray data points - 1)
    ;; flag chisqr_pdf value == 1
    iix2 = finite(p_x2) and p_x2 gt 0.
    if keyword_set(chisq) then begin
        ibest = where(x2 eq min(x2[where(iix2,/null)]),nbest)
        if (nbest gt 1) then stop;ibest = ibest[0]
        stat = [a2[ibest],p_a2[ibest],x2[ibest],p_x2[ibest],0.,0.]
    endif else begin
        ;; combined test statistic
        x2_joint = -2.*(alog(p_a2)+alog(p_x2))
        p_joint = dblarr(nfree)
        p_joint[where(iia2 and iix2,/null)] = 1.-chisqr_pdf(x2_joint[where(iia2 and iix2,/null)],4.)  ;; dof = 2k, where k == 2, the number of tests being combined
        ;; choose best model
        ibest = where(x2_joint eq min(x2_joint[where(iia2 and iix2,/null)]),nbest)
        if (nbest gt 1) then ibest = ibest[0]
        stat = [a2[ibest],p_a2[ibest],x2[ibest],p_x2[ibest],x2_joint[ibest],p_joint[ibest]]
    endelse
    ;; halt if no best fit found
    if (nbest eq 0) then stop
    ;; record best fit statistics
    f24v_[n] = f24_steps[ibest]
    f25v_[n] = f25_steps[ibest]
    statv_[*,n] = stat
    nhmv_[*,n] = nh_mod[*,ibest]
    rxmv_[*,n] = rx_mod[*,ibest]
    iimv_[*,n] = iimod[*,ibest]
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - VARIABLE CTF'
    endif else if (n mod (ncount/10) eq 0) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - VARIABLE CTF'
print, '=============================================='

sav_vars = ['F24_STEPS','F25_STEPS', $
            'RXDV_', $
            'F24V_','F25V_','STATV_','NHMV_','RXMV_','A2V_']
sav_inds = ['IIMV_']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="fct_variable.sav"')


END






