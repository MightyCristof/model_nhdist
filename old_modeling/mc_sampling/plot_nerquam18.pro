




;; iagn_NOX
p1 = plot(xnox,ynox,'orange',/stairstep, $
         /fill_background,fill_color='orange',fill_transparency=75,name='Candidate AGN w.o. X-ray')
;; Lansbury+15 sampled LL
p2 = plot(xsampll,ysampll,/stairstep,/ov,/fill_background,fill_transparency=25,name='Lansbury+15')
;; iagn_NUX
p3 = errorplot([-1.1,xhist,0.5],[0,yhist_norm,0],[0,yhist_err,0],'o',color='purple', $
              linestyle='-',/stairstep,/ov,/fill_background,fill_color='purple',fill_transparency=60,name='Candidate AGN w. X-ray')
;; sample corrected LL
p4 = plot([-1.1400000,outx,0.65999991],[0.,outy,0.],/ov,/stairstep,linestyle='-','teal',/fill_background,fill_color='teal',fill_transparency=75,name='$Corrected N_H model$')
p4.xr=[-2,1]
p4.yr=[0,1.2]
leg = legend(target=[p2,p1,p3,p4],position=[-1.,1.0],/data,/auto_text_color,sample_width=.1)

p4.xtitle = '$Log(L_{10-40 keV}[observed] / L_{10-40 keV}[expected])$'
p4.ytitle = 'Frequency'
p4.save,'LL_model.png'



readcol,'lansbury+15_nhdist.csv',ylan,delim=',',format='x,d'
ylan = [0,ylan,0]
xlan = [19.5,findgen(5,start=20),24.5]
p1 = plot(xlan,ylan/max(ylan*1.),/stairstep,yra=[0.,1.2],xra=[19.26],/fill_background,fill_transparency=25,name='Lansbury+15')
p2 = plot(xoutnh,youtnh/max(youtnh*1.),/stairstep,color='teal',/ov,/fill_background,fill_color='teal',fill_transparency=75,name='$Corrected N_H model$')
leg = legend(target=[p1,p2],position=[21.5,1.],/data,/auto_text_color,sample_width=.1)
p2.xtitle = '$Log(N_H/cm^{-2})$'
p2.ytitle = 'Frequency'
p2.save,'NH_model.png'