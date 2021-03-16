PRO plot_model_nhdist, NHMOD = nhmod, $
                       NHAVG = nhavg, $
                       HIDE = hide, $
                       LOW_RES = low_res, $
                       SAV = sav


common _data
common _nhobs
common _rxnh
common _group
common _ctfest
common _split
common _rxmod
common _fxest


if keyword_set(low_res) then res = 100 else res = 600

file_mkdir,'figures'
;;----------------------------------------------------------------------------------------
;; MODELED NH DISTRUBITUION w/ STACKING RESULTS
;;----------------------------------------------------------------------------------------
if keyword_set(nhmod) then begin
    ;; X-ray stack image path
    soft = file_search('../data_prep/nodet_wagn_052_stack.png')
    hard = file_Search('../data_prep/nodet_wagn_27_stack.png')
    
    ;; histogram normalization to match Ananna+2019
    ksnorm = total(nh_mod2_ks[where(nh_mod2_ks.xh lt 24.)].yh)
    adnorm = total(nh_mod2_ad[where(nh_mod2_ad.xh lt 24.)].yh)
    
    ;; X-ray fluxes    
    energy_center = [1.25,4.5]
    energy_range = [0.75,2.5]
    ksoff = [0.25,1.1]     
    adoff = [0.25,1.5]
    stack_flux = fluxv[1,0:1]
    stack_err = flux_errv[1,0:1]
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    dim = [840,740]
    sq = 180
    gap = 40
    pos = [[80,540,80+sq,540+sq],[80+sq+gap,540,80+2.*sq+gap,540+sq],[560,540,800,540+sq],$
                            [80,70,800,480]]
    e = {xra:[20.,26.],yra:[0.,1.4],$
         stairstep:1, $
         xtitle:'$log !8N!7_{H} [cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times', $
         dim:dim,device:1,buffer:0}
    if keyword_set(hide) then e.buffer = 1
    
    kscol = 'teal';[65,182,196];[99,172,190];
    adcol = 'purple';[37,52,148];[96,26,74];
    
    plo = plot(nh_ana_lo.xh,nh_ana_lo.yh,'-.',thick=2,_extra=e,pos=pos[*,3],name='$Ananna+2019 (log !8L!7_X < 43.6)$');,fill_background=0,fill_transparency=50)
    phi = plot(nh_ana_hi.xh,nh_ana_hi.yh,':',thick=2,_extra=e,/ov,name='$Ananna+2019 (log !8L!7_X > 43.6)$');fill_background=0,fill_transparency=50,/ov)
    pks = errorplot(nh_mod2_ks.xh+nh_mod2_ks.xoff,nh_mod2_ks.yh/ksnorm,nh_mod2_ks.sig/ksnorm, $
                    '__',thick=4,_extra=e,fill_background=1,fill_color=kscol,fill_transparency=50,/ov,name='This work (KS test)')
    pad = errorplot(nh_mod2_ad.xh+nh_mod2_ad.xoff,nh_mod2_ad.yh/adnorm,nh_mod2_ad.sig/adnorm, $
                   '-',thick=4,_extra=e,fill_background=1,fill_color=adcol,fill_transparency=50,/ov,name='This work (AD test)')
    ;; add legend
    l = legend(target=[plo,phi,pks,pad],position=[0.125,0.61],/normal,horizontal_alignment=0.,font_size=12,font_name='Times')
    ;; add X-ray stacked images
    psoft = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,0],/current,/device,ytitle='offset in Decl. [arcsec.]',font_name='Times')
    imsoft = image(soft,pos=pos[*,0],/current,/device)    
    lsoft = text(target=psoft,0.5,0.82,'  0.5$-$2 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    phard = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,1],/current,/device);,xtitle='offset in R.A. [arcsec.]                     ')
    imhard = image(hard,pos=pos[*,1],/current,/device)
    lhard = text(target=phard,0.5,0.82,'  2$-$7 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    xt = text(target=plo,(pos[2,0]+pos[0,1])/2.,pos[1,0]-42,/device,'offset in R.A. [arcsec.]',alignment=0.5,font_name='Times')
    ;; add X-ray flux estimates
    p = errorplot(energy_center,stack_flux,energy_range,stack_err, $
                   linestyle='',/xlog,/ylog,pos=pos[*,2],/current,/device, $
                   xtitle='$energy [keV]$',ytitle='$!8F!7_{2-10 keV}  [erg s^{-1} cm^{-2}]$', $
                   font_name='Times')
                   ;xtickunit='exponent'
    px = plot(energy_center,stack_flux,'s',col='black', $
              sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',/ov,name='Stack')
    p.xra=[1e-1,1e1]
    p.yra=[1e-16,1e-14]
    p = errorplot(energy_center-ksoff,[fx_soft_ksx,fx_hard_ksx],[e_fx_soft_ksx,e_fx_hard_ksx], $
                  linestyle='',errorbar_capsize=0.1,/ov)
    pks = plot(energy_center-ksoff,[fx_soft_ksx,fx_hard_ksx],'S',col='black', $
               sym_size=1.8,sym_thick=1.5,sym_filled=1,sym_fill_color=kscol,sym_transparency=0,/ov,name='KS test')
    p = errorplot(energy_center+adoff,[fx_soft_adx,fx_hard_adx],[e_fx_soft_adx,e_fx_hard_adx], $
                  linestyle='',errorbar_capsize=0.1,/ov)
    pad = plot(energy_center+adoff,[fx_soft_adx,fx_hard_adx],'o',col='black', $
               sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=adcol,sym_transparency=0,/ov,name='AD test')
    ;; add "legend"    
    po = plot([1.5e-1],[5.75e-15],'S',col='black',sym_size=1.8,sym_thick=1.5,sym_filled=1,sym_fill_color=kscol,sym_transparency=0,/ov)
    to = text(0.15,0.85,target=pks,/relative,'KS test',font_name='Times')
    po = plot([1.5e-1],[3.65e-15],'o',col='black',sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=adcol,sym_transparency=0,/ov)
    to = text(0.15,0.75,target=pad,/relative,'AD test',font_name='Times')
    po = plot([1.5e-1],[2.3e-15],'s',col='black',sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    to = text(0.15,0.65,target=px,/relative,'X-ray stack',font_name='Times')

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/nhmod.eps',/BITMAP else $
                                                     p.save,'figures/nhmod.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; MODELED NH DISTRUBITUION w/ STACKING RESULTS --- ANANNA+2019 AVERAGED NH
;;----------------------------------------------------------------------------------------
if keyword_set(nhavg) then begin
    ;; fraction of hi- and low- LX sources
    lxfrac = 1
    if keyword_set(lxfrac) then suffix = '_frac' else $
                                suffix = ''

    ;; X-ray stack image path
    soft = file_search('../data_prep/nodet_wagn_052_stack.png')
    hard = file_Search('../data_prep/nodet_wagn_27_stack.png')
    
    ;; histogram normalization to match Ananna+2019
    ksnorm = total(nh_mod2_ks[where(nh_mod2_ks.xh lt 24.)].yh)
    adnorm = total(nh_mod2_ad[where(nh_mod2_ad.xh lt 24.)].yh)
    ;; anything NH<20, add to NH=20-21 bin
    inhlt20 = where(nh_mod2_ks.xh lt 20.,nnhlt20)
    if (nnhlt20 eq 1) then begin
        nh_mod2_ks[inhlt20+1].yh += total(total(nh_mod2_ks[inhlt20].yh))
        nh_mod2_ks[inhlt20].yh = 0.
        nh_mod2_ad[inhlt20+1].yh += total(total(nh_mod2_ad[inhlt20].yh))
        nh_mod2_ad[inhlt20].yh = 0.
    endif
    
    ;; X-ray fluxes    
    energy_center = [1.25,4.5]
    energy_range = [0.75,2.5]
    modoff = [-0.25,-1.1]     
    stack_flux = fluxv[1,0:1]
    stack_err = flux_errv[1,0:1]
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    dim = [840,740]
    sq = 180
    gap = 40
    pos = [[80,540,80+sq,540+sq],[80+sq+gap,540,80+2.*sq+gap,540+sq],[560,540,800,540+sq],$
                            [80,70,800,480]]
    e = {xra:[20.,26.],yra:[0.,1.1],$
         stairstep:1, $
         xtitle:'$log !8N!7_{H} [cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times', $
         dim:dim,device:1,buffer:0}
    if keyword_set(hide) then e.buffer = 1
    
    kscol = 'teal';[65,182,196];[99,172,190];
    adcol = 'purple';[37,52,148];[96,26,74];
    
    ;; KS TEST
    pks = errorplot(nh_mod2_ks.xh+nh_mod2_ks.xoff,nh_mod2_ks.yh/ksnorm,nh_mod2_ks.sig/ksnorm, $
                    '-',thick=4,_extra=e,pos=pos[*,3],fill_background=1,fill_color=kscol,fill_transparency=50,name='This work')
    if keyword_set(lxfrac) then pavg = plot(nh_ana_lo.xh,nh_ana_lo.yh*frac_lox+nh_ana_hi.yh*frac_hix,'--',thick=2,_extra=e,/ov,name='$Ananna+2019 (averaged)$') else $
                                pavg = plot(nh_ana_lo.xh,(nh_ana_lo.yh+nh_ana_hi.yh)/2.,'--',thick=2,_extra=e,/ov,name='$Ananna+2019 (averaged)$')
    ;; CT fraction
    ctks = text(25.,0.35,'$!8f!7_{CT} = '+string(ctf2_ks,format='(d4.2)')+'$',/data,font_size=16,font_name='Times',alignment=0.5)
    ;; add legend
    l = legend(target=[pavg,pks],position=[0.125,0.61],/normal,horizontal_alignment=0.,font_size=12,font_name='Times')
    ;; add X-ray stacked images
    psoft = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,0],/current,/device,ytitle='offset in Decl. [arcsec.]',font_name='Times')
    imsoft = image(soft,pos=pos[*,0],/current,/device)    
    lsoft = text(target=psoft,0.5,0.82,'  0.5$-$2 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    phard = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,1],/current,/device);,xtitle='offset in R.A. [arcsec.]                     ')
    imhard = image(hard,pos=pos[*,1],/current,/device)
    lhard = text(target=phard,0.5,0.82,'  2$-$7 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    xt = text(target=plo,(pos[2,0]+pos[0,1])/2.,pos[1,0]-42,/device,'offset in R.A. [arcsec.]',alignment=0.5,font_name='Times')
    ;; add X-ray flux estimates
    p = errorplot(energy_center,stack_flux,energy_range,stack_err, $
                   linestyle='',/xlog,/ylog,pos=pos[*,2],/current,/device, $
                   xtitle='$energy [keV]$',ytitle='$!8F!7_{2-10 keV}  [erg s^{-1} cm^{-2}]$', $
                   font_name='Times')
                   ;xtickunit='exponent'
    px = plot(energy_center,stack_flux,'s',col='black', $
              sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',/ov,name='Stack')
    p.xra=[1e-1,1e1]
    p.yra=[1e-16,1e-14]
    p = errorplot(energy_center+modoff,[fx_soft_ksx,fx_hard_ksx],[e_fx_soft_ksx,e_fx_hard_ksx], $
                  linestyle='',errorbar_capsize=0.1,/ov)
    pks = plot(energy_center+modoff,[fx_soft_ksx,fx_hard_ksx],'S',col='black', $
               sym_size=1.8,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    pks = plot(energy_center+modoff,[fx_soft_ksx,fx_hard_ksx],'S',col='black', $
               sym_size=1.8,sym_thick=1.5,sym_filled=1,sym_fill_color=kscol,sym_transparency=20,/ov)
    ;; manual "legend"    
    po = plot([1.5e-1],[5.75e-15],'S',col='black',sym_size=1.8,sym_thick=1.5,sym_filled=1,sym_fill_color=kscol,sym_transparency=20,/ov)
    to = text(0.15,0.85,target=pks,/relative,'Modeled flux',font_name='Times')
    po = plot([1.5e-1],[3.65e-15],'s',col='black',sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    to = text(0.15,0.75,target=px,/relative,'X-ray stack',font_name='Times')

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/nhavg_ks'+suffix+'.eps',/BITMAP else $
                                                     p.save,'figures/nhavg_ks'+suffix+'.png',resolution=res
    endif
    
    ;; AD TEST
    e.yra=[0.,1.]
    pad = errorplot(nh_mod2_ad.xh+nh_mod2_ad.xoff,nh_mod2_ad.yh/adnorm,nh_mod2_ad.sig/adnorm, $
                   '-',thick=4,_extra=e,pos=pos[*,3],fill_background=1,fill_color=adcol,fill_transparency=40,name='This work')
    if keyword_set(lxfrac) then pavg = plot(nh_ana_lo.xh,nh_ana_lo.yh*frac_lox+nh_ana_hi.yh*frac_hix,'--',thick=2,_extra=e,/ov,name='$Ananna+2019 (averaged)$') else $
                               pavg = plot(nh_ana_lo.xh,(nh_ana_lo.yh+nh_ana_hi.yh)/2.,'--',thick=2,_extra=e,/ov,name='$Ananna+2019 (averaged)$')
    ;; CT fraction
    ctad = text(25.,0.35,'$!8f!7_{CT} = '+string(ctf2_ad,format='(d4.2)')+'$',/data,font_size=16,font_name='Times',alignment=0.5)
    ;; add legend
    l = legend(target=[pavg,pad],position=[0.125,0.61],/normal,horizontal_alignment=0.,font_size=12,font_name='Times')
    ;; add X-ray stacked images
    psoft = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,0],/current,/device,ytitle='offset in Decl. [arcsec.]',font_name='Times')
    imsoft = image(soft,pos=pos[*,0],/current,/device)    
    lsoft = text(target=psoft,0.5,0.82,'  0.5$-$2 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    phard = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,1],/current,/device);,xtitle='offset in R.A. [arcsec.]                     ')
    imhard = image(hard,pos=pos[*,1],/current,/device)
    lhard = text(target=phard,0.5,0.82,'  2$-$7 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    xt = text(target=plo,(pos[2,0]+pos[0,1])/2.,pos[1,0]-42,/device,'offset in R.A. [arcsec.]',alignment=0.5,font_name='Times')
    ;; add X-ray flux estimates
    p = errorplot(energy_center,stack_flux,energy_range,stack_err, $
                   linestyle='',/xlog,/ylog,pos=pos[*,2],/current,/device, $
                   xtitle='$energy [keV]$',ytitle='$!8F!7_{2-10 keV}  [erg s^{-1} cm^{-2}]$', $
                   font_name='Times')
                   ;xtickunit='exponent'
    px = plot(energy_center,stack_flux,'s',col='black', $
              sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',/ov,name='Stack')
    p.xra=[1e-1,1e1]
    p.yra=[1e-16,1e-14]
    p = errorplot(energy_center+modoff,[fx_soft_adx,fx_hard_adx],[e_fx_soft_adx,e_fx_hard_adx], $
                  linestyle='',errorbar_capsize=0.1,/ov)
    pad = plot(energy_center+modoff,[fx_soft_adx,fx_hard_adx],'o',col='black', $
               sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    pad = plot(energy_center+modoff,[fx_soft_adx,fx_hard_adx],'o',col='black', $
               sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=adcol,sym_transparency=20,/ov)
    ;; add "legend"    
    ;po = plot([1.5e-1],[5.75e-15],'S',col='black',sym_size=1.8,sym_thick=1.5,sym_filled=1,sym_fill_color=kscol,sym_transparency=20,/ov)
    ;to = text(0.15,0.85,target=pks,/relative,'Modeled flux',font_name='Times')
    po = plot([1.5e-1],[5.75e-15],'o',col='black',sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=adcol,sym_transparency=20,/ov)
    to = text(0.15,0.85,target=pad,/relative,'Modeled flux',font_name='Times')
    ;po = plot([1.5e-1],[2.3e-15],'s',col='black',sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    ;to = text(0.15,0.65,target=px,/relative,'X-ray stack',font_name='Times')
    po = plot([1.5e-1],[3.65e-15],'s',col='black',sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    to = text(0.15,0.75,target=px,/relative,'X-ray stack',font_name='Times')
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/nhavg_ad'+suffix+'.eps',/BITMAP else $
                                                     p.save,'figures/nhavg_ad'+suffix+'.png',resolution=res
    endif
endif










END


