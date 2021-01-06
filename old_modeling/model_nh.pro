

;; load NuSTAR Analysis 
common _nustar1
common _nustar2

;; create a random two-peaked NH distribution
nsamp = n_elements(iagn)
nlo = nsamp/100.
nhi = nsamp/(100./99.)
modlo = nh_mc(nlo)
modhi = nh_mc(nhi)+2.5
;nhmod = [modlo,modhi[where(modhi gt 25.)]]
nhmod = [modlo,modhi]

;remove everything N_H > 1e26
ibh = where(nhmod gt 27.)
nhmod[ibh] = randomu(seed,n_elements(ibh))*3.+24.

nmod = n_elements(nhmod)


nhbin = freedman(nhmod)
hnhm = histogram(nhmod,locations=xnhm,bin=nhbin)
hnhmlo = histogram(modlo,locations=xnhmlo,bin=nhbin,min=min(nhmod),max=max(nhmod))
hnhmhi = histogram(modhi,locations=xnhmhi,bin=nhbin,min=min(nhmod),max=max(nhmod))
e = {xtitle:'$N_H$',dimension:[1200,400],margin:0.2}
hnh = histplot(xnhm,hnhm,_extra=e,layout=[3,1,1])
hnh = histplot(xnhmlo,hnhmlo,col='blue',ov=hnh)
hnh = histplot(xnhmhi,hnhmhi,col='red',ov=hnh)

;; convert model NH distribution to model luminosity ratio
lscat = 0.3
llmod = ll2nh(nhmod,/get_lum) + randomn(seed,nmod)*lscat -0.3



;; histogram limits and bin size
llmm = [min(llmod)<min(lllim[iagn]),max(llmod)>max(lllim[iagn])]
llbin = freedman(lllim[iagn])

;; model LL histogram
hllm = histogram(llmod,locations=xllm,min=llmm[0],max=llmm[1],bin=llbin)
hllmlo = histogram(llmod[0:nlo-1],locations=xllmlo,bin=llbin,min=llmm[0],max=llmm[1])
hllmhi = histogram(llmod[nlo:-1],locations=xllmhi,bin=llbin,min=llmm[0],max=llmm[1])
e = {xtitle:'$L(obs)/L(int)$',margin:0.2}
hll = histplot(xllm,hllm,_extra=e,layout=[3,1,2],/current)
hll = histplot(xllmlo,hllmlo,col='blue',ov=hll)
hll = histplot(xllmhi,hllmhi,col='red',ov=hll)

;; randomly match and compare llmod to ll-limits
;; such that:    llmod[i] ==> ll[iagn[imll[i]]]
;; create a detected flag such that:    detect[i] == 1 if llmod > ll-lim
llen = n_elements(iagn)
imll = fix(randomu(seed,nmod)*llen)
detect = llmod gt lllim[iagn[imll]]
id = where(detect)
iu = where(~detect)
;; match limit
im = iagn[imll[iu]]
;; histograms!
;; limits and detections VS. model limits and detections
hllmd = histogram(llmod[id],locations=xllmd,bin=llbin,min=llmm[0],max=llmm[1])
;hllmodu = histogram(llmod[where(~detect)],locations=xllmodu,min=llmm[0],max=llmm[1],bin=h)
hllmu = histogram(lllim[im],locations=xllmu,bin=llbin,min=llmm[0],max=llmm[1])

;hp = histplot(xllmod,hllmod,col='teal',_extra=e)
e = {xtitle:'$Model L(obs)/L(int)$',margin:0.2}
hm = histplot(xllm,hllm,col='black',_extra=e,layout=[3,1,3],/current)
hm = histplot(xllmd,hllmd,col='teal',ov=hm)
hm = histplot(xllmu,hllmu,col='orange',ov=hm)

print, 'Detected:   '+strtrim(n_elements(where(detect)),2)
print, 'Undetected: '+strtrim(n_elements(where(~detect)),2)


hd = histplot(xlldet,hlldet,col='dodger blue',_extra=e,layout=[2,1,1])
hd = histplot(xllmd,hllmd,col='teal',ov=hd)
hd.title = '$Detections$'
hu = histplot(xlllim,hlllim,col='orange',_extra=e,layout=[2,1,2],/current)
hu = histplot(xllmu,hllmu,col='teal',ov=hu)
hu.title = '$Non-detections$'















END





