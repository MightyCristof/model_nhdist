PRO test_mc_sampling


yval = 3.
yerr = 1.
niter = 10000

noise_list = randomn(seed, niter)
noise_list = yerr*noise_list
y_list = yval + noise_list
x_list = 2*y_list

hist = histogram(y_list,binsize=0.1,min=-1.,locations=xx)
;xx = findgen(n_elements(hist))*0.1 -1.
fit = gaussfit(xx,hist,Acoeffs,nterms=3)
plot,xx,hist,xrange=[-2,10]
oplot,xx,fit,color=FSC_COLOR('blue'),thick=3
print,'y: ',Acoeffs[1],Acoeffs[2]

hist = histogram(x_list,binsize=0.1,min=-1.,locations=xx)
;xx = findgen(n_elements(hist))*0.1 -1.
fit = gaussfit(xx,hist,Acoeffs,nterms=3)
oplot,xx,hist,color=FSC_COLOR('red')
oplot,xx,fit,color=FSC_COLOR('green'),thick=3.
print,'x: ',Acoeffs[1],Acoeffs[2]



;; HOW TO RUN THIS ON MY NH DISTRIBUTION

readcol,'lansbury+15_nhdist.csv',ylan,delim=',',format='x,d'
xlan = findgen(5,start=20)
ylan = [0,ylan/max(ylan),0]
xlan = [19,xlan,25]
fit = gaussfit(xlan,ylan,coeff,nterms=nterms)

niter = 10000
noise = randomu(seed, niter)
xfit = noise + coeff[1]
yfit = coeff[2]*noise + mean(fit)

y = fit[round(x*n_elements(fit))]
hist = histogram(y_list,binsize=0.1,min=-1.,locations=xx)

p = plot(xlan,ylan,/stairstep,yra=[0,1.2])
p = plot(xlan,fit,'red',/ov)
p = plot(xx,hist/max(hist),'green',/stairstep,/ov)
END



yfit = randomu(seed,niter)*coeff[2] + mean(fit)




;; SPLINE METHOD
probabilities = [0, 1, 10, 2, 0, 0, 4, 5, 3, 1, 0]
x = [1:nel(probabilities)]

xq = [1:nel(probabilities):0.05]
pdf = spline(x,probabilities,xq)
pdf = pdf>0.
pdf /= total(pdf)
cdf = total(pdf,/cumulative)
iu = uniq(cdf)
xx = xq[iu]
cdf = cdf[iu]
randomvalues = randomu(seed,2500)
projection = interpol(xx,cdf,randomvalues)
proj = histogram(projection,bin=.1,locations=xproj)
p = plot(x,probabilities,xtitle='Original Distribution',yst=2,/stairstep,layout=[1,3,1])
p = plot(xq,pdf,'blue',xtitle = 'Probability Distribution',/current,layout=[1,3,2])
p = plot(xproj,proj,'green',xtitle = 'Randomly Sampled Distribution',/stairstep,/current,layout=[1,3,3])

;; TRY IT WITH A SPLINE INSTEAD
readcol,'lansbury+15_nhdist.csv',ylan,delim=',',format='x,d'
ylan = [0,ylan,0]
xlan = findgen(7,start=19)

xq = [xlan[0]:xlan[-1]:0.05]
pdf = spline(xlan,ylan,xq)
pdf = pdf>0.
pdf /= total(pdf)
cdf = total(pdf,/cumulative)
iu = uniq(cdf)
xx = xq[iu] 
cdf = cdf[iu]
niter = 1e4
randval = randomu(seed,niter)
projection = interpol(xx,cdf,randval)
proj = histogram(projection,bin=scott(projection),locations=xproj)
p = plot(xlan,ylan,xtitle='Original Distribution',yr=[0,1.2],/stairstep,layout=[1,3,1])
p = plot(xq,pdf,'blue',xtitle = 'Probability Distribution',/current,layout=[1,3,2])
p = plot(xproj,proj,'green',xtitle = 'Randomly Sampled Distribution',/stairstep,/current,layout=[1,3,3])



;; Lansbury+15 luminosity ratio
readcol,'lansbury+15_lx_lir.csv',lir_lan,lx_lan,delim=',',format='d,d'
ll_lan = lx_lan/lir_lan
plothist,ll_lan,/autobin,peak=1,xhist,yhist,/noplot
xq = [xhist[0]:xhist[-1]:1e-4]
pdf = spline(xhist,yhist,xq)
pdf = pdf>0.
pdf /= total(pdf)
cdf = total(pdf,/cumulative)
iu = uniq(cdf)
xx = xq[iu]
cdf = cdf[iu]
niter = 1e5
randval = randomu(seed,niter)
projection = interpol(xx,cdf,randval)
projhist = histogram(projection,bin=scott(projection),locations=xproj)
p = plot(xhist,yhist,xtitle='$Original L_X/L_{IR} Distribution$',yra=[0,1.2],/stairstep,layout=[1,3,1])
p = plot(xq,pdf,'blue',xtitle = 'Probability Distribution',/current,layout=[1,3,2])
p = plot(xproj,projhist,'green',xtitle = 'Randomly Sampled Distribution',/stairstep,/current,layout=[1,3,3])












