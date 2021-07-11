FUNCTION ad_test, data, $
                  model, $
                  PERMUTE = permute, $
                  PROB = p, $
                  CVM = cvm, $
                  WEIGHT = weight

                 
FORWARD_FUNCTION ad_test

sz = n_elements(data)
if sz lt 2 then message, 'data must contain at least two elements'

edf = edf(data,xloc=xdat)
edf_mod = edf(model,xloc=xmod)

;; if weighting arrays with -9999, remove to continue or summing issues arise
if keyword_set(weight) then begin
    idat = where(xdat ne -9999.,dlen)
    if (dlen gt 0) then begin
        xdat = xdat[idat]
        edf = edf[idat]
    endif
    imod = where(xmod ne -9999.,mlen)
    if (mlen gt 0) then begin
        xmod = xmod[imod]
        edf_mod = edf_mod[imod]
    endif
endif

cdf = interpol(edf_mod,xmod,xdat)
;; correct interpolation outside of data bounds
ilt = where(cdf lt 0.,nlt)
if (nlt gt 0) then cdf[ilt] = 0.;interpol([0.,cdf[ilt[-1]+1]],[xdat[ilt[0]],xdat[ilt[-1]+1]],xdat[ilt])
igt = where(cdf gt 1.,ngt)
if (ngt gt 0) then cdf[igt] = 1.;interpol([cdf[igt[0]-1],1.],[xdat[igt[0]-1],xdat[igt[-1]]],xdat[igt])

;; squared distance between Fn and Fx
delta = edf - cdf
;; weighting (w = 1 for CvM, w = Fx(1-Fx) for AD)
w = cdf * (1. - cdf)
if keyword_set(cvm) then w = 1.

;ks = max(abs(delta), /nan)
;ky = max(delta, /nan) + max(-delta, /nan)
a2 = delta^2. / w
a2 = total(a2, /nan) / total(finite(a2))
;mad = total(abs(delta), /nan) / total(finite(delta))

;; probability placeholder
nreps = 1000
p = 1./nreps

if keyword_set(permute) then begin
    ndata = n_elements(data)
    ad_perm = dblarr(nreps)
    joint = [data,model]
    njoint = n_elements(joint)
    ;; permutation test (well, my quick approximation of one anyway)
    for i = 0,nreps-1 do begin
        p_joint = joint[randomi(njoint,njoint,/nodup)]
        p_dat = p_joint[0:ndata-1]
        p_mod = p_joint[ndata:-1]                
        ad_perm[i] = ad_test(p_dat,p_mod,weight=weight)
    endfor

    p = total(ad_perm gt a2)/nreps > 1./nreps
endif

return, a2


END

                  

