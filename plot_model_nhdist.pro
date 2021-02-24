plot_model_nhdist, COMP_NH = comp_nh


common _data
common _nhobs
common _ctfest
common _split
common _rxmod
common _fxest


if keyword_set(postdoc_apps)
    ymod = histogram(nh_mod2_ks[*,iks2],locations=xmod,bin=1.,min=19.,max=27.)
    c_ana = total(nh_ana_lo.yh)
    c_me = c_ana/total(ymod)

    e = {xra:[20,26],yra:[0,1.45], $
         font_name:'Times',font_size:14, $
         xtitle:'$log !8N!7_{H} [cm^{-2}]$', $
         stairstep:1,thick:4, $
         buffer:0}
    im1 = image('../data_prep/Ananna+2019_Fig10 LOW.jpg',transparency=60,dimension=[550,440],position=[50,65,530,415],/device)
    ;p = plot(nh_ana_lo.xh,nh_ana_lo.yh,/stairstep,xra=,yra=[,col='orange',/current,position=im1.position)
    p = plot(xmod+width(xmod)/2.,ymod*c_me,'__',col='red',_extra=e,/current,position=im1.position,name='Carroll et al. (in prep.)')
    l = legend(target=p,/auto_text_color,position=[22.9,0.88],/data,HORIZONTAL_SPACING=0.06) 
    t = text(0.18,0.45,"PRELIMINARY",col='red',font_name='Times',font_style='Bold',font_size=16)
    p.save,'nh_dist.png'
endif


if keyword_set(comp_nh) then begin
    yks = histogram(nh_mod2_ks.(iks2),location=xks,min=20.,max=26.,bin=1.)
    yad = histogram(nh_mod2_ad.(iad2),location=xad,min=20.,max=26.,bin=1.)
    p = plot(nh_ana_lo.xh,nh_ana_lo.yh,/stairstep,xra=[20.,26.],yra=[0.,1.6],thick=4,xtitle='$N_H$')
    p = plot(nh_ana_hi.xh,nh_ana_hi.yh,/stairstep,'--',thick=4,/ov)
    p = plot(xks+0.5,yks/total(yks[where(xks lt 24.)]),/stairstep,col='teal','__',/ov,thick=4)
    p = plot(xad+0.5,yad/total(yad[where(xad lt 24.)]),/stairstep,col='purple','-.',/ov,thick=4)
    t = text(21.,1.45,'Model: 1% scatter ',/data)
    t = text(21.,1.3,'$log F_X [KS]: '+string(median(logfx_med_ksv),format='(f6.2)')+'$',/data,col='teal')
    t = text(21.,1.15,'$log F_X [AD]: '+string(median(logfx_med_adv),format='(f6.2)')+'$',/data,col='purple')
endif










END


