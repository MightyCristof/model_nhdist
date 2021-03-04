PRO plot_model_nhdist, APPS = apps, COMP_NH = comp_nh


common _data
common _nhobs
common _rxnh
common _group
common _ctfest
common _split
common _rxmod
common _fxest


if keyword_set(apps) then begin
    ymod_ks = histogram(nh_mod_ks,locations=xmod_ks,bin=1.,min=19.,max=27.)
    ymod_ad = histogram(nh_mod_ad,locations=xmod_ad,bin=1.,min=19.,max=27.)
    ymod2_ks = histogram(nh_mod2_ks[*,iks2],locations=xmod2_ks,bin=1.,min=19.,max=27.)
    ymod2_ad = histogram(nh_mod2_ad[*,iad2],locations=xmod2_ad,bin=1.,min=19.,max=27.)

    ;; ANANNA LX LOW
    c_lo = total(nh_ana_lo.yh)    
    c_ks = c_lo/total(ymod_ks)
    c_ad = c_lo/total(ymod_ad)
    c_ks2 = c_lo/total(ymod2_ks)
    c_ad2 = c_lo/total(ymod2_ad)

    elo = {xra:[20,26],yra:[0,1.5], $
           font_name:'Times',font_size:14, $
           xtitle:'$log !8N!7_{H} [cm^{-2}]$', $
           stairstep:1,thick:4, $
           buffer:0}
    im1 = image('../data_prep/Ananna+2019_Fig10 LOW.jpg',transparency=60,dimension=[550,440],position=[50,65,530,415],/device)
    ;p = plot(nh_ana_lo.xh,nh_ana_lo.yh,_extra=elo,/current,position=im1.position)
    pks = plot(xmod_ks+width(xmod_ks)/2.,ymod_ks*c_ks,'__',col='teal',_extra=elo,/current,position=im1.position,name='KS')
    pad = plot(xmod_ad+width(xmod_ad)/2.,ymod_ad*c_ad,'__',col='purple',_extra=elo,/current,position=im1.position,name='AD')

    pks2 = plot(xmod2_ks+width(xmod2_ks)/2.,ymod2_ks*c_ks2,'__',col='teal',_extra=elo,/current,position=im1.position,name='KS')
    pad2 = plot(xmod2_ad+width(xmod2_ad)/2.,ymod2_ad*c_ad2,'__',col='purple',_extra=elo,/current,position=im1.position,name='AD')
    l = legend(target=pks,/auto_text_color,position=[22.9,0.88],/data,HORIZONTAL_SPACING=0.06) 
    t = text(0.18,0.45,"PRELIMINARY",col='red',font_name='Times',font_style='Bold',font_size=16)
    ;p.save,'nh_dist.png'

    ;; ANANNA LX HIGH
    c_hi = total(nh_ana_hi.yh)
    c_ks = c_hi/total(ymod_ks)
    c_ad = c_hi/total(ymod_ad)
    c_ks2 = c_hi/total(ymod2_ks)
    c_ad2 = c_hi/total(ymod2_ad)


    ehi = {xra:[20,26],yra:[0,0.7], $
           font_name:'Times',font_size:14, $
           xtitle:'$log !8N!7_{H} [cm^{-2}]$', $
           stairstep:1,thick:4, $
           buffer:0}
    im1 = image('data_prep/Ananna+2019_Fig10 HI.jpg',transparency=60,dimension=[550,440],position=[50,65,530,415],/device)
    ;p = plot(nh_ana_hi.xh,nh_ana_hi.yh,_extra=ehi,/current,position=im1.position)
    pks = plot(xmod_ks+width(xmod_ks)/2.,ymod_ks*c_ks,'__',col='teal',_extra=ehi,/current,position=im1.position,name='KS')
    pad = plot(xmod_ad+width(xmod_ad)/2.,ymod_ad*c_ad,'__',col='purple',_extra=ehi,/current,position=im1.position,name='AD')

    pks2 = plot(xmod2_ks+width(xmod2_ks)/2.,ymod2_ks*c_ks2,'__',col='teal',_extra=ehi,/current,position=im1.position,name='KS')
    pad2 = plot(xmod2_ad+width(xmod2_ad)/2.,ymod2_ad*c_ad2,'__',col='purple',_extra=ehi,/current,position=im1.position,name='AD')
    l = legend(target=pks,/auto_text_color,position=[22.9,0.88],/data,HORIZONTAL_SPACING=0.06) 
    t = text(0.18,0.45,"PRELIMINARY",col='red',font_name='Times',font_style='Bold',font_size=16)
    ;p.save,'nh_dist.png'



endif
stop

if keyword_set(nh_mod) then begin
    e = {xra:[20.,26.],yra:[0.,1.5],$
         xtitle='$!8N!7_H$',$
         font_size:16,font_name:'Times'}

    p = plot(nh_ana_lo.xh,nh_ana_lo.yh,/stairstep,xra=[20.,26.],yra=[0.,1.5],thick=4,xtitle='$!8N!7_H$')
    p = plot(nh_ana_hi.xh,nh_ana_hi.yh,/stairstep,'--',thick=4,/ov)
    p = plot(xks+0.5,yks/total(yks[where(xks lt 24.)]),/stairstep,col='teal','__',/ov,thick=4)
    p = plot(xad+0.5,yad/total(yad[where(xad lt 24.)]),/stairstep,col='purple','-.',/ov,thick=4)
    t = text(21.,1.3,'$log F_X [KS]: '+string(mean(median(logfx_full_ksv,dim=2)),format='(f6.2)')+'$',/data,col='teal')
    t = text(21.,1.15,'$log F_X [AD]: '+string(mean(median(logfx_full_adv,dim=2)),format='(f6.2)')+'$',/data,col='purple')
endif



        
if keyword_set(paper) then begin


    yks = histogram(nh_mod2_ksv[*,iks2],location=xks,min=20.,max=26.,bin=1.)
    yad = histogram(nh_mod2_adv[*,iad2],location=xad,min=20.,max=26.,bin=1.)
    
    elo = {xra:[20,26],yra:[0,1.5], $
           font_name:'Times',font_size:14, $
           xtitle:'$log !8N!7_{H} [cm^{-2}]$', $
           stairstep:1,thick:4, $
           buffer:0}
    imlo = image('../data_prep/Ananna+2019_Fig10 LOW.jpg',transparency=60,dimension=[550,440],position=[50,65,530,415],/device)
    p = plot(nh_ana_lo.xh,nh_ana_lo.yh,'-',col='purple',_extra=elo,/current,position=imlo.position,name='XLO');/stairstep,xra=[20.,26.],yra=[0.,1.6],thick=4,xtitle='$N_H$')
    p = plot(nh_ana_hi.xh,nh_ana_hi.yh,'--',col='purple',_extra=elo,/ov,name='XHI')    
    p = plot(xks+0.5,yks/total(yks[where(xks lt 24.)]),/stairstep,col=col[*,0],'__',_extra=elo,/ov,name='KS')
    p = plot(xad+0.5,yad/total(yad[where(xad lt 24.)]),/stairstep,col=col[*,1],'-.',_extra=elo,/ov,name='AD')
    
    p = plot()    
    
endif   







END


