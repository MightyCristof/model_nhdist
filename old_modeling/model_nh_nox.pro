PRO model_nh_nox, lof, hif

;; load NuSTAR Analysis 
common _xray_soft
common _xray_hard
common _xray_nox



;; create a random two-peaked NH distribution
nsamp = n_elements(iagn_nox)
nlo = nsamp/lof
nhi = nsamp/hif
modlo = nh_mc(nlo)
modhi = nh_mc(nhi)+2.5
;nhmod = [modlo,modhi[where(modhi gt 25.)]]
nhmod = [modlo,modhi]

;remove everything N_H > 1e26
ibh = where(nhmod gt 27.)
nhmod[ibh] = randomu(seed,n_elements(ibh))*3.+24.

nmod = n_elements(nhmod)


;nhbin = freedman(nhmod)
nhbins = scott(nhmod)
hnhm = histogram(nhmod,locations=xnhm,bin=nhbin)
hnhmlo = histogram(modlo,locations=xnhmlo,bin=nhbin,min=min(nhmod),max=max(nhmod))
hnhmhi = histogram(modhi,locations=xnhmhi,bin=nhbin,min=min(nhmod),max=max(nhmod))
e = {xtitle:'$Model log( N_H )$',buffer:1}
hnh = histplot(xnhm,hnhm,_extra=e)
hnh = histplot(xnhmlo,hnhmlo,col='blue',ov=hnh)
hnh = histplot(xnhmhi,hnhmhi,col='red',ov=hnh)
hnh.save,'modeled_nh.png'

;; convert model NH distribution to model luminosity ratio
lscat = 0.3
llmod = ll2nh(nhmod,/get_lum) + randomn(seed,nmod)*lscat



;; histogram limits and bin size
llmm = [min(llmod)<min(lllim_hb[iagn_nox]),max(llmod)>max(lllim_hb[iagn_nox])]
;llbin = freedman(lllim[iagn])
llbin = scott(lllim_hb[iagn_nox])

;; model LL histogram
hllm = histogram(llmod,locations=xllm,min=llmm[0],max=llmm[1],bin=llbin)
hllmlo = histogram(llmod[0:nlo-1],locations=xllmlo,bin=llbin,min=llmm[0],max=llmm[1])
hllmhi = histogram(llmod[nlo:-1],locations=xllmhi,bin=llbin,min=llmm[0],max=llmm[1])
e = {xtitle:'$Model L_{obs}/L_{int}$',buffer:1}
hll = histplot(xllm,hllm,_extra=e)
hll = histplot(xllmlo,hllmlo,col='blue',ov=hll)
hll = histplot(xllmhi,hllmhi,col='red',ov=hll)
hll.save,'modeled_ll1.png'

;; randomly match and compare llmod to ll-limits
;; such that:    llmod[i] ==> ll[iagn[imll[i]]]
;; create a detected flag such that:    detect[i] == 1 if llmod > ll-lim
llen = n_elements(iagn_nox)
imll = fix(randomu(seed,nmod)*llen)
detect = llmod gt lllim_hb[iagn_nox[imll]]+.5
id = where(detect)
iu = where(~detect)
;; match limit
im = iagn_nox[imll[iu]]
;; histograms!
;; limits and detections VS. model limits and detections
hllmd = histogram(llmod[id],locations=xllmd,bin=llbin,min=llmm[0],max=llmm[1])
;hllmodu = histogram(llmod[where(~detect)],locations=xllmodu,min=llmm[0],max=llmm[1],bin=h)
hllmu = histogram(lllim_hb[im],locations=xllmu,bin=llbin,min=llmm[0],max=llmm[1])

;hp = histplot(xllmod,hllmod,col='teal',_extra=e)
e = {xtitle:'$Model L_{obs}/L_{int}$'}
hm = histplot(xllm,hllm,col='black',_extra=e,buffer=1)
hm = histplot(xllmd,hllmd,col='teal',ov=hm)
hm = histplot(xllmu,hllmu,col='orange',ov=hm)
hm.save,'modeled_ll2.png'


print, 'Detected:   '+strtrim(n_elements(where(detect)),2)
print, 'Undetected: '+strtrim(n_elements(where(~detect)),2)

e = {xra:[-3,1],xtitle:'$Model L_{obs}/L_{int}$'}
hd = histplot(xlldet_hb,hlldet_hb,col='dodger blue',_extra=e,layout=[2,1,1])
hd = histplot(xllmd,hllmd,col='teal',ov=hd)
hd.title = '$Detections$'
hu = histplot(xlllim_hb,hlllim_hb,col='orange',_extra=e,/current,layout=[2,1,2])
hu = histplot(xllmu,hllmu,col='teal',ov=hu)
hu.title = '$Non-detections$'
hd.save,'modeled_ll_det.png'




END





