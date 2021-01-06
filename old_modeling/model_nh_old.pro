;PRO model_nh

;;IM READY UPDATE ME
;; create a random two-peaked NH distribution
nsamp = 20000.
modlo = nh_mc(nsamp/100.)-1.
modhi = nh_mc(nsamp/2.)+1.
modnh = [modlo,modhi]
nmod = n_elements(modnh)
;; convert model NH distribution to luminosity ratio
lscat = 0.3
modrl = rl2nh(modnh,model='BORUS',/rl_out) + randomn(seed,n_elements(modnh))*lscat
ymodrl = histogram(modrl,locations=xmodrl,binsize=scott(modrl))
;; randomly match and compare MODRL to RL such that:    modrl[i] ==> RL[irand[i]]
;; create a detected flag such that:    detect[i] == 1 if modrl > LL
;rl = rldet>rlnon
;nsrc = n_elements(rl)
;irand = randomi()

inon = where(xnon,nsrc)
irand = randomi(nmod,nsrc-1)
dmod = modrl gt rl[inon[irand]]
iii = irand[where(dmod)]
;; histogram luminosity ratio of limits
yrl = histogram(rl[inon],locations=xrl,bin=freedman(rlnon[inon]))
yrld = histogram(modrl[where(dmod)],locations=xrld,min=min(xmodrl),bin=width(xmodrl),nbins=n_elements(xmodrl))
yrln = histogram(rl[inon[iii]],locations=xrlu,min=min(xmodrl),bin=width(xmodrl),nbins=n_elements(xmodrl))


e = {xtitle:'$!8R!7_L$',fill_transparency:75}
hp = histplot(xmodrl,ymodrl,fill_color='light grey',_extra=e)
hp = histplot(xrld,yrld,fill_color='dodger blue',ov=hp)
hp = histplot(xrlu,yrlu,fill_color='orange',ov=hp)




END





