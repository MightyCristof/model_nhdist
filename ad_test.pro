FUNCTION ad_test, data, $
                  model, $
                  PERMUTE = permute, $
                  PROB = p, $
                  CVM = cvm, $
                  CDFM = cdf

                 
FORWARD_FUNCTION ad_test

sz = n_elements(data)
if sz lt 2 then message, 'data must contain at least two elements'

ndata = n_elements(data)
nmodel = n_elements(model)

edf = edf(data,xloc=xdata)
edf_model = edf(model,xloc=xmodel)
cdf = interpol(edf_model,xmodel,xdata)
;; correct interpolation outside of data bounds
ilt = where(cdf lt 0.,nlt)
if (nlt gt 0) then cdf[ilt] = 0.;interpol([0.,cdf[ilt[-1]+1]],[xdata[ilt[0]],xdata[ilt[-1]+1]],xdata[ilt])
igt = where(cdf gt 1.,ngt)
if (ngt gt 1) then cdf[igt] = 1.;interpol([cdf[igt[0]-1],1.],[xdata[igt[0]-1],xdata[igt[-1]]],xdata[igt])

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
    ad_perm = dblarr(nreps)
    joint = [data,model]
    njoint = n_elements(joint)
    ;for i = 0,nreps-1 do ad_perm[i] = ad_test(data,joint[randomi(nmodel,njoint,/nodup)])

    ;; permutation test (well, my quick approximation of one anyway)
    for i = 0,nreps-1 do begin
        pjoint = joint[randomi(njoint,njoint,/nodup)]
        pdata = pjoint[0:ndata-1]
        pmodel = pjoint[ndata:-1]                
        ad_perm[i] = ad_test(pdata,pmodel)
    endfor

    p = total(ad_perm gt a2)/nreps > 1./nreps
endif

return, a2


END

                  

