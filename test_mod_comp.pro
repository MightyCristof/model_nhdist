PRO test_mod_comp, TEST = test


common _group
common _model


dim = [700,300]
ind = randomi(100,n_elements(rx_modv[0,*]),/nodup)

for i = 0,99 do begin
    yy = edf(rx_modv[where(iimodv[*,ind[i]]),ind[i]],xloc=xx)
    if (i eq 0) then p = plot(xx,yy,transparency=90,dim=dim,layout=[2,1,1]) else $
                     p = plot(xx,yy,transparency=90,/ov)
endfor
yd = edf(rxd,xloc=xd)
p = plot(xd,yd,col='dodger blue',thick=2,/ov)

for i = 0,99 do begin
    yy_ = edf(rx_modv_[where(iimodv_[*,ind[i]]),ind[i]],xloc=xx_)
    if (i eq 0) then p = plot(xx_,yy_,transparency=90,current=1,layout=[2,1,2]) else $
                     p = plot(xx_,yy_,transparency=90,/ov)
endfor
p = plot(xd,yd,col='dodger blue',thick=2,/ov)




;; so anyway...
if keyword_set(oldcomparisonthatsnotvalidanymorebutyoumightneedatsomepoint) then begin
    if (n_elements(test) eq 0) then test = 1

    case test of 
        3: begin
            gfile = 'variant_models/'+['src_wac/select_group.sav','src_wac_lownh/select_group.sav','select_group.sav']
            ffile = 'variant_models/'+['src_wac/ctf_estimate2.sav','src_wac_lownh/ctf_estimate2.sav','rx_model.sav']
            mfile = 'variant_models/'+['src_wac/rx_model.sav','src_wac_lownh/rx_model.sav','rx_model.sav']
            var = ['rx_mod2_adv','rx_mod2_adv','rx_mod2v']
            ivar = ['iimod2_adv','iimod2_adv','iimod2v']
            title = ['Lansbury+2015 N$_H$','Added low N$_H$','Post-modeled low N$_H$']
            dim = [1000,300]
            end
        2: begin
            gfile = 'variant_models/'+['src_wac/select_group.sav','src_wac_lownh/select_group.sav']
            ffile = 'variant_models/'+['src_wac/ctf_estimate2.sav','src_wac_lownh/ctf_estimate2.sav']
            mfile = 'variant_models/'+['src_wac/rx_model.sav','src_wac_lownh/rx_model.sav']
            var = ['rx_mod2_adv','rx_mod2_adv']
            ivar = ['iimod2_adv','iimod2_adv']
            title = ['Lansbury+2015 N$_H$','Added low N$_H$']
            dim = [700,300]
            end
        1: begin
            gfile = 'variant_models/'+['select_group.sav']
            ffile = 'variant_models/'+['ctf_fixed.sav']
            mfile = 'variant_models/'+['rx_model.sav']
            var = ['rx_mod2_adv','rx_mod2_adv']
            ivar = ['iimod2_adv','iimod2_adv']
            title = ['Lansbury+2015 N$_H$','Added low N$_H$']
            dim = [400,300]
            end
        else: message, "Incorrect number of test cases."
    endcase

    ntest = n_elements(gfile)

    marg = [0.18,0.15,0.12,0.1]
    curr = [0,1,1]
    for f = 0,n_elements(gfile)-1 do begin
        restore,gfile[f]
        restore,ffile[f]
        restore,mfile[f]
        ;niter = n_elements(rx_modv[0,*])
        re = execute('niter = n_elements('+var[f]+'[0,*])')
        ndraw = 100
        irand = randomi(ndraw,niter,/nodup)
        for i = 0,n_elements(irand)-1 do begin
            ;edf_model = edf(rx_modv[where(iimodv[*,irand[i]] eq 1),irand[i]],xloc=x_model)
            re = execute('index = where('+ivar[f]+'[*,irand[i]] eq 1)')
            re = execute('edf_model = edf('+var[f]+'[index,irand[i]],xloc=x_model)')
            if (i eq 0) then pmodel = plot(x_model,edf_model,transparency=90,xra=[-3.5,1], $
                                           xtitle='$R_X$',ytitle='Empirical Dist. Func.',name='Models', $
                                           layout=[ntest,1,f+1],margin=marg,dim=dim,current=curr[f]) else $
                             pmodel = plot(x_model,edf_model,/ov,transparency=90,name='Models')
        endfor
        edf_data = edf(rxd,xloc=x_data)
        pdata = plot(x_data,edf_data,col='red',thick=2,name='Data',/ov,title=title[f])
        l = legend(target=[pdata,pmodel],position=[0.45,0.94],/relative)
    endfor

    for f = 0,n_elements(gfile)-1 do begin
        restore,gfile[f]
        restore,ffile[f]
        restore,mfile[f]
        ;niter = n_elements(rx_modv[0,*])
        re = execute('niter = n_elements('+var[f]+'[0,*])')
        irand = randomi(ndraw,niter,/nodup)
        for i = 0,n_elements(irand)-1 do begin
            ;yhm = histogram(rx_modv[where(iimodv[*,irand[i]] eq 1),irand[i]],bin=freedman(rx_modv[where(iimodv[*,irand[i]] eq 1),irand[i]]),location=xhm)
            re = execute('index = where('+ivar[f]+'[*,irand[i]] eq 1)')
            re = execute('yhm = histogram('+var[f]+'[index,irand[i]],bin=freedman('+var[f]+'[index,irand[i]]),location=xhm)')
            if (i eq 0) then pmodel = plot(xhm+width(xhm)/2.,yhm,/stairstep,transparency=90,yra=[0,80],xtitle='$R_X$',ytitle='Frequency',name='Model', $
                                           layout=[ntest,1,f+1],margin=marg,current=curr[f],dim=dim) else $
                             pmodel = plot(xhm+width(xhm)/2.,yhm,/stairstep,transparency=90,/ov,name='Model')   
        endfor
        yhd = histogram(rxd,bin=scott(rxd),location=xhd)
        pdata = plot(xhd+width(xhd)/2.,yhd,/stairstep,col='red',thick=2,/ov,name='Data',title=title[f]) 
        l = legend(target=[pdata,pmodel],position=[0.45,0.94],/relative)
    endfor
endif


END
