PRO plot_model_nhdist, PROPERTIES = properties, $
                       NHDIST = nhdist, $
                       RXDIST = rxdist, $
                       XSPEC = xspec, $
                       NHMODEL = nhmodel, $
                       RXMODEL = rxmodel, $
                       RXNH = rxnh, $
                       STACK = stack, $
                       NEWNH = newnh, $
                       HIDE = hide, $
                       LOW_RES = low_res, $
                       SAV = sav


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
common _uniform
;common _variable
common _model


if keyword_set(low_res) then res = 100 else res = 600
file_mkdir,'figures'


;;----------------------------------------------------------------------------------------
;; DISTRIBUTIONS OF DATA PROPERTIES
;;----------------------------------------------------------------------------------------
if keyword_set(properties) then begin

    ;; shorthand indices
    ix = where(iiwac)
    id = where(iiwd)
    in = where(iiwn)
    
    ;; REDSHIFT
    bn = 0.05
    yh = histogram(z[ix],bin=bn,location=xh,min=0.,max=0.8)
    yhd = histogram(z[id],bin=bn,location=xhd,min=0.,max=0.8)
    yhn = histogram(z[in],bin=bn,location=xhn,min=0.,max=0.8)

    e = {xra:[-0.05,0.85],yra:[0.,120.],$
         stairstep:1,thick:2,fill_background:1, $
         xtitle:'$!8z!7$',ytitle:'Frequency', $
         font_name:'Times',font_size:14, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
    
    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources',fill_transparency=30)
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.',fill_transparency=20)
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected',fill_transparency=15)
    leg = legend(target=[p,pd,pn],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=14,font_name='Times')
    ;t = text(0.03,0.90,"a",/normal,font_name='Times',font_style='Bold',font_size=16)
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/z_dist.eps',/BITMAP else $
                                                     p.save,'figures/z_dist.png',resolution=res
    endif    

    ;; LX
    bn = 0.25
    yh = histogram(loglx[ix],bin=bn,location=xh,min=41.,max=46.)
    yhd = histogram(loglx[id],bin=bn,location=xhd,min=41.,max=46.)
    yhn = histogram(loglx[in],bin=bn,location=xhn,min=41.,max=46.)
    
    e = {xra:[41.5,45.5],yra:[0.,80.],$
         stairstep:1,thick:2,fill_background:1, $
         xtitle:'$log !8L!7_{X}  [erg s^{-1} cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times',font_size:14, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources',fill_transparency=30)
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.',fill_transparency=20)
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected',fill_transparency=15)
    leg = legend(target=[p,pd,pn],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=14,font_name='Times')
    ;leg = legend(target=[pd],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=12,font_name='Times')
    ;t = text(0.03,0.90,"b",/normal,font_name='Times',font_style='Bold',font_size=16)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/lx_dist.eps',/BITMAP else $
                                                     p.save,'figures/lx_dist.png',resolution=res
    endif    

    ;; LIR
    bn = 0.25
    yh = histogram(loglir[ix],bin=bn,location=xh,min=41.,max=46.)
    yhd = histogram(loglir[id],bin=bn,location=xhd,min=41.,max=46.)
    yhn = histogram(loglir[in],bin=bn,location=xhn,min=41.,max=46.)
    
    e = {xra:[41.5,45.5],yra:[0.,200.],$
         stairstep:1,thick:2,fill_background:1, $
         xtitle:'$log !8L!7_{MIR}  [erg s^{-1} cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times',font_size:14, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources',fill_transparency=30)
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.',fill_transparency=20)
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected',fill_transparency=15)
    leg = legend(target=[p,pd,pn],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=14,font_name='Times')
    ;t = text(0.03,0.90,"c",/normal,font_name='Times',font_style='Bold',font_size=16)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/lir_dist.eps',/BITMAP else $
                                                     p.save,'figures/lir_dist.png',resolution=res
    endif    

    ;; LX(LIR)
    bn = 0.25
    yh = histogram(loglxir[ix],bin=bn,location=xh,min=41.,max=46.)
    yhd = histogram(loglxir[id],bin=bn,location=xhd,min=41.,max=46.)
    yhn = histogram(loglxir[in],bin=bn,location=xhn,min=41.,max=46.)
    
    e = {xra:[41.5,45.5],yra:[0.,200.],$
         stairstep:1,thick:2,fill_background:1, $
         xtitle:'$log !8L!7_{X}(!8L!7_{MIR})  [erg s^{-1} cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times',font_size:14, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources',fill_transparency=30)
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.',fill_transparency=20)
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected',fill_transparency=15)
    leg = legend(target=[p,pd,pn],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=14,font_name='Times')
    ;t = text(0.03,0.90,"d",/normal,font_name='Times',font_style='Bold',font_size=16)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/lxir_dist.eps',/BITMAP else $
                                                     p.save,'figures/lxir_dist.png',resolution=res
    endif    

endif


;;----------------------------------------------------------------------------------------
;; MODELED NH DISTRUBITUION w/ STACKING RESULTS --- ANANNA+2019 AVERAGED NH
;;----------------------------------------------------------------------------------------
if keyword_set(nhdist) then begin

    ;; X-ray stack image path
    prep_dir = file_search('data_prep')
    up = ''
    while (prep_dir eq '') do begin
        up = up+'../'
        prep_dir = file_search(up+'data_prep')
    endwhile
    soft = file_search(prep_dir+'/nodet_wagn_052_stack.png')
    hard = file_Search(prep_dir+'/nodet_wagn_27_stack.png')

    ;; X-ray stacked fluxes    
    energy_center = [1.25,4.5]
    energy_range = [0.75,2.5]
    modoff = [[0.80,0.80],[1.25,1.25]]
    stack_flux = fxstak[1,0:1]
    stack_err = e_fxstak[1,0:1]
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    ;; variables for FIXED vs FREE
    fx_models = ['FX_NON','FX_NON_']
    nh_models = ['NH_MOD','NH_MOD_']
    ctf_val = ['FCT','FCT_']
    
    ;; plot set up
    dim = [880,780]
    sq = 180
    gap = 40
    pos = [[100,560,100+sq,560+sq],[100+sq+gap,560,100+2.*sq+gap,560+sq],[600,560,840,560+sq],$
                                   [100,70,840,480]]
    e = {font_name:'Times', $
         dim:dim,device:1,buffer:0}
    enh = {xra:[20.,26.],yra:[0.,1.0], $
           stairstep:1, $
           xtitle:'$log !8N!7_{H}  [cm^{-2}]$',ytitle:'Frequency (normalized; $!8N!7_H < 10^{24} cm^{-2}$)', $;ytitle:'Frequency (normalized)', $
           font_name:'Times',font_size:14, $
           current:1,dim:dim,device:1}
    if keyword_set(hide) then e.buffer = 1

    ;; color
    col = [106,81,163];[84,39,143]
            ;; purple     [84,39,143]
            ;; teal       [28,144,153]
            ;; mint       [127,205,187]
            ;; aquamarine [28,144,153]

    ;; BEGIN FIGURE
    ;; plot X-ray stacked images
    psoft = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,0],_extra=e,/device,ytitle='offset in Decl. [arcsec.]',font_size=14)
    isoft = image(soft,pos=pos[*,0],/current,/device)    
    tsoft = text(target=psoft,0.5,0.82,'  0.5$-$2 keV  ',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    phard = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,1],/current,/device);,xtitle='offset in R.A. [arcsec.]                     ')
    ihard = image(hard,pos=pos[*,1],/current,/device)
    thard = text(target=phard,0.5,0.82,'  2$-$7 keV  ',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    xt = text((pos[2,0]+pos[0,1])/2.,pos[1,0]-42,/device,'offset in R.A. [arcsec.]',alignment=0.5,font_name='Times',font_size=14)

    ;; plot X-ray flux estimates
    pxd = errorplot(energy_center,stack_flux,energy_range,stack_err, $
                   linestyle='',/xlog,/ylog,pos=pos[*,2],/current,/device, $
                   xtitle='$energy  [keV]$',ytitle='$!8F!7_{X}  [erg s^{-1} cm^{-2}]$', $
                   font_name='Times',font_size=14)
    pxd = plot(energy_center,stack_flux,'s',col='black', $
               sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',/ov,name='Stack')
    pxd.xra=[1e-1,1e1]
    pxd.yra=[1e-16,1e-14]
;    for i = 0,1 do begin
;        re = execute('fxmod = '+fx_models[i])
;        ;; FULL FLUX
;        fxm = [mode(fxmod.soft,kde=kde_bandwidth(fxmod.soft)),mode(fxmod.hard,kde=kde_bandwidth(fxmod.hard))]
;        e_fxm = [mode(fxmod.e_soft,kde=kde_bandwidth(fxmod.e_soft)),mode(fxmod.e_hard,kde=kde_bandwidth(fxmod.e_hard))]
;        ;fxm = [mean(fxmod.soft),mean(fxmod.hard)]
;        ;e_fxm = [mean(fxmod.e_soft),mean(fxmod.e_hard)]
;        if (i eq 0) then begin  ;; uniform model flux
;            pxm = errorplot(energy_center*modoff[*,i],fxm,e_fxm, $
;                            linestyle='',errorbar_capsize=0.1,/ov)
;            pxm = plot(energy_center*modoff[*,i],fxm,'o',col='black', $
;                       sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
;            pxm = plot(energy_center*modoff[*,i],fxm,'o',col='black', $
;                       sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=col,sym_transparency=20,fill_transparency=20,/ov)
;        endif else begin    ;; variable model flux
;            pxm = errorplot(energy_center*modoff[*,i],fxm,e_fxm, $
;                            linestyle='',errorbar_capsize=0.1,/ov)
;            pxm = plot(energy_center*modoff[*,i],fxm,'o',col='black', $
;                       sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
;            pxm = plot(energy_center*modoff[*,i],fxm,'o',col=col, $
;                       sym_size=1.3,sym_thick=1.5,sym_filled=0,sym_fill_color=col,sym_transparency=30,/ov)
;            pxm = plot(energy_center*modoff[*,i],fxm,'o',col=col, $
;                       sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,fill_transparency=0,/ov)
;        endelse
;    endfor    
    to = text(0.06,0.85,target=pxm,/relative,'X-ray non-det. sources',font_name='Times')
    po = plot([1.5e-1],[3.65e-15],'o',col='black',sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=col,sym_transparency=20,/ov)
    to = text(0.15,0.75,target=pxm,/relative,'Model (uniform)',font_name='Times')
    po = plot([1.5e-1],[2.32e-15],'o',col='black',sym_size=1.4,sym_thick=1.5,sym_filled=0,sym_fill_color=col,sym_transparency=0,/ov)
    po = plot([1.5e-1],[2.32e-15],'o',col=col,sym_size=1.3,sym_thick=1.5,sym_filled=0,sym_fill_color=col,sym_transparency=30,/ov)
    po = plot([1.5e-1],[2.32e-15],'o',col=col,sym_size=1.2,sym_thick=1.5,sym_filled=0,sym_fill_color='white',sym_transparency=0,/ov)
    to = text(0.15,0.65,target=pxm,/relative,'Model (variable)',font_name='Times')
    po = plot([1.5e-1],[1.47e-15],'s',col='black',sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    to = text(0.15,0.55,target=pxd,/relative,'X-ray stack',font_name='Times')
    
    ;; plot NH distribution
    pct = plot([24.,26.],[1.,1.]*enh.yra[1],_extra=enh,pos=pos[*,3], $
               fill_background=1,fill_color='light grey',fill_level=0.,fill_transparency=0)
    pct.stairstep = 0   ;; CT shading
    stop
    for i = 0,1 do begin
        re = execute('nhm = '+nh_models[i])
        ilo = where(nhm.xh lt 20.)
        nhm[ilo[-1]+1].yh += total(nhm[ilo].yh)
        nhm[ilo].yh = 0.
        nhm[ilo[-1]+1].sig += sqrt(total(nhm[ilo[0]:ilo[-1]+1].sig^2.))
        nhm[ilo].sig = 0.
        if (i eq 0) then nhm[where(nhm.xh eq 24.)].yh = nhm[where(nhm.xh eq 25.)].yh
        if (i eq 0) then begin  ;; uniform nh dist
            re = execute('ctf = string('+ctf_val[i]+',format="(d0.2)")+"\pm"+string(E_'+ctf_val[i]+',format="(d0.2)")')
            pnhu = errorplot(nhm.xhoff,nhm.yh,nhm.sig,'-',thick=4,errorbar_thick=4,_extra=enh, $
                             pos=pos[*,3],fill_background=1,fill_color=col,fill_transparency=30,name='Model (uniform)')
        endif else begin    ;; variable nh dist
            pnhv = plot(nhm.xhoff,nhm.yh,':',thick=2,_extra=enh, $
                        pos=pos[*,3],fill_background=0,fill_color=col,fill_transparency=30,name='Model (variable)')
        endelse
    
    endfor
    ;; Ananna+2019
    pan = plot(nh_ana_lox.xh,nh_ana_lox.yh*frac_lox+nh_ana_hix.yh*frac_hix,'__',thick=4,_extra=enh,/ov,name='$Ananna+2019 (weighted)$')
    tct = text(0.683,0.91,'CT',target=pct,/relative,font_name='Times',font_size=14)
    ctad = text(25.,enh.yra[1]/4.,'$!8f!7_{CT} = '+ctf+'$',target=pnhu,/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5)
    leg = legend(target=[pan,pnhu,pnhv],position=[0.125,0.60],/normal,horizontal_alignment=0.,sample_width=0.14,font_size=14,font_name='Times')

    ;; figure reference for Nature caption
    ;t = text(0.03,0.97,"a",/normal,font_name='Times',font_style='Bold',font_size=16)
    ;t = text(0.03,0.62,"b",/normal,font_name='Times',font_style='Bold',font_size=16)
    
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then pct.save,'figures/nh_dist.eps',/BITMAP else $
                                                     pct.save,'figures/nh_dist.png',resolution=res
    endif    
endif


;;----------------------------------------------------------------------------------------
;; DISTRIBUTION OF BEST-FIT RX FROM MODELING
;;----------------------------------------------------------------------------------------
if keyword_set(rxdist) then begin

    data = ['RXDV','RXDV_']
    models = ['RXMV','RXMV_']
    indices = ['IIMV','IIMV_']
    ;models = ['RX_MODV','RX_MODV_']
    ;indices = ['IIMODV','IIMODV_']
    letter = ['a','b']
    file = ['rx_uniform','rx_variable']
    
    rxbin = 0.2d
    
    ;for i = 0,1 do begin
    for i = 0,0 do begin
        re = execute('dat = '+data[i])
        re = execute('model = '+models[i])
        re = execute('index = '+indices[i])
        
        ;; data
        dat = hist2d_avg(dat,rxbin,normalize=1)
        ;; models
        model = hist2d_avg(model,rxbin,iidet=index,normalize=1,confidence=1)

        e = {xra:[-2.6,1.],yra:[0.,0.3],$
             stairstep:1,fill_background:1, $
             xtitle:'$!8R!7_{!8L!7_X}$',ytitle:'Frequency (normalized)', $
             font_name:'Times',font_size:14, $
             buffer:0}
        if keyword_set(hide) then e.buffer = 1

        ;; blue
        ;col1 = [4,90,141]
        ;col2 = [54,144,192]
        ;col3 = [116,169,207]
        ;; purple
        col1 = [106,81,163]
        col2 = [128,125,186]
        col3 = [188,189,220]
        
        p3u = plot(model.xhoff,model.sig3u_det,_extra=e,color='black',transparency=75,fill_color=col3,fill_transparency=0)
        p2u = plot(model.xhoff,model.sig2u_det,_extra=e,color='black',transparency=75,fill_color=col2,fill_transparency=0,/ov)
        p1u = plot(model.xhoff,model.sig1u_det,_extra=e,color='black',transparency=75,fill_color=col1,fill_transparency=0,/ov)
        p1l = plot(model.xhoff,model.sig1l_det,_extra=e,color='black',transparency=75,fill_color=col2,fill_transparency=0,/ov)
        p2l = plot(model.xhoff,model.sig2l_det,_extra=e,color='black',transparency=75,fill_color=col3,fill_transparency=0,/ov)
        p3l = plot(model.xhoff,model.sig3l_det,_extra=e,color='black',transparency=75,fill_color='white',/ov)
        p = plot(model.xhoff,model.yh_det,/stairstep,thick=2,/ov,name='$!8R!7_X^{model}$')
        ;pd = errorplot(xhd[indd]+0.1,yhd[indd]*normd,ehd[indd]*normd,'o',sym_size=1.,sym_filled=1,/ov,name='$!8R!7_{!8L!7_X} (data)$')
        pd = errorplot(dat.xhoff,dat.yh,dat.sig,'o',sym_size=1.,sym_filled=1,/ov,name='$!8R!7_{!8L!7_X} (data)$')
        plim = plot([1.,1.]#rx2nh(23.,/rx_out),e.yra,'--',thick=2,/ov)    
        psym3 = plot([0.,0.],[-1.,-1.],'s',sym_size=1.5,sym_filled=1,col='black',sym_fill_color=col3,name='3$\sigma$ conf.',/ov)
        psym2 = plot([0.,0.],[-1.,-1.],'s',sym_size=1.5,sym_filled=1,col='black',sym_fill_color=col2,name='2$\sigma$ conf.',/ov)
        psym1 = plot([0.,0.],[-1.,-1.],'s',sym_size=1.5,sym_filled=1,col='black',sym_fill_color=col1,name='1$\sigma$ conf.',/ov)
        leg = legend(target=[psym1,psym2,psym3,pd],position=[0.35,0.85],/normal,sample_width=0.,horizontal_spacing=0.06,font_name='Times',font_size=14)
        ;; NH limit
        t = text(0.62,0.8,/normal,'$!8N!7_H \geq$ 10$^{23}$ cm$^{-2}$',alignment=1.,font_name='Times',font_size=14)
        p = text(0.62,0.75,/normal,'$\leftarrow$',alignment=1.,font_size=26)
        ;; figure reference for Nature caption
        ;t = text(0.03,0.90,letter[i],/normal,font_name='Times',font_style='Bold',font_size=16)

        ;if keyword_set(sav) then begin
        ;    print, '    SAVING PLOT'
        ;    if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/'+file[i]+'.eps',/BITMAP else $
        ;                                                 p.save,'figures/'+file[i]+'.png',resolution=res
        ;endif        
        if keyword_set(sav) then begin
            print, '    SAVING PLOT'
            if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/rx_dist.eps',/BITMAP else $
                                                         p.save,'figures/rx_dist.png',resolution=res
        endif        
    endfor    
endif


;;----------------------------------------------------------------------------------------
;; XSPEC MODEL WITH INCREASING NH
;;----------------------------------------------------------------------------------------
if keyword_set(xspec) then begin

    ;; XSPEC model component output
    file = '/Users/ccarroll/Research/projects/model_nhdist/workspace/data_prep/borus_modelcomps.txt'

    ;; read data
    readcol,file,energy,tot21,bor21,pl21,scpl21,tot22,bor22,pl22,scpl22,tot23,bor23,pl23,scpl23,tot24,bor24,pl24,scpl24,tot25,bor25,pl25,scpl25,/silent
    
    ;; concat variables
    ;nhs = ['21','22','23','24','25']
    nhs = ['21','23','25']
    nnh = n_elements(nhs)
    tot = 'tot'+nhs
    bor = 'bor'+nhs
    pl = 'pl'+nhs
    scpl = 'scpl'+nhs

    ;; plot variables
    letter = ['a','b','c']
    line = ['-','__','-.',':','__']
    ;col = [[213,94,0],[204,121,167],[0,114,178],[240,228,66],[0,158,115]]

    ;;      BORUS         SCPL        PL         TOTAL
    ;;      MAGENTA       GREEN       ORANGE     BLUE
    col = [[204,121,167],[0,158,115],[213,94,0],[0,114,178]]

    e = {xra:[0.1,300],yra:[1e-4,100], $
         xlog:1,ylog:1, $
         xtickname:['0.1','1','10','100'],ytickname:'10!U'+strtrim([-4:2],2), $
         ytickinterval:1, $
         xtitle:'Energy  [keV]',ytitle:'keV$^2$  [Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$]', $
         font_name:'Times',font_size:14}
    thick = 2
    
    for i = 0,n_elements(nhs)-1 do begin
        p = plot(energy,tot21,_extra=e,/nodata)
        re = execute('pb = plot(energy,'+bor[i]+',linestyle=line[3],col=col[*,0],thick=thick,/ov,name="Torus")')
        re = execute('ps = plot(energy,'+scpl[i]+',linestyle=line[2],col=col[*,1],thick=thick,/ov,name="Scattered")')
        re = execute('pp = plot(energy,'+pl[i]+',linestyle=line[1],col=col[*,2],thick=thick,/ov,name="Primary")')
        re = execute('pt = plot(energy,'+tot[i]+',linestyle=line[0],col=col[*,3],thick=thick,/ov,name="Total emission")')
        leg = legend(target=[pt,pp,ps,pb],position=[0.410,0.866],sample_width=0.16,horizontal_spacing=0.06,font_name='Times',font_size=14)
        tx = text(0.8,0.76,'!8N!7$_H = 10^{'+nhs[i]+'} cm^{-2}$',/relative,alignment=1.,font_name='Times',font_style='Bold',font_size=16)
        ;t = text(0.03,0.90,letter[i],/normal,font_name='Times',font_style='Bold',font_size=16)

        if keyword_set(sav) then begin
            print, '    SAVING PLOT'
            if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/xspec_nh'+nhs[i]+'.eps',/BITMAP else $
                                                         p.save,'figures/xspec_nh'+nhs[i]+'.png',resolution=res
        endif                
    endfor
endif
;    ;; component SEDs
;    re = execute('pb = plot(energy,'+bor[0]+',linestyle=line[0],col=col[*,0],thick=thick,/ov,name="Torus")')
;    for i = 1,nnh-1 do re = execute('p = plot(energy,'+bor[i]+',linestyle=line[i],col=col[*,0],thick=thick,/ov)')
;    re = execute('ps = plot(energy,'+scpl[0]+',linestyle=line[0],col=col[*,1],thick=thick,/ov,name="Scattered")')
;    for i = 1,nnh-1 do re = execute('p = plot(energy,'+scpl[i]+',linestyle=line[i],col=col[*,1],thick=thick,/ov)')
;    re = execute('pp = plot(energy,'+pl[0]+',linestyle=line[0],col=col[*,2],thick=thick,/ov,name="Primary")')
;    for i = 1,nnh-1 do re = execute('p = plot(energy,'+pl[i]+',linestyle=line[i],col=col[*,2],thick=thick,/ov)')
;    ;; total SEDs
;    re = execute('pt21 = plot(energy,'+tot[0]+',linestyle=line[0],col=col[*,3],thick=thick,/ov,name="Total (!8N!7$_H = 10^{21} cm^{-2})$")')
;    re = execute('pt23 = plot(energy,'+tot[1]+',linestyle=line[1],col=col[*,3],thick=thick,/ov,name="Total (!8N!7$_H = 10^{23} cm^{-2})$")')
;    re = execute('pt25 = plot(energy,'+tot[2]+',linestyle=line[2],col=col[*,3],thick=thick,/ov,name="Total (!8N!7$_H = 10^{25} cm^{-2})$")')
;    ;for i = 0,nnh-1 do re = execute('pt = plot(energy,'+tot[i]+',linestyle=line[i],col=col[*,0],/ov)')
;    ;; legend
;    leg = legend(target=[pt21,pt23,pt25],position=[0.48,0.87],horizontal_spacing=0.06,font_name='Times',font_size=14)
;    leg = legend(target=[pp,ps,pb],position=[0.85,0.30],horizontal_alignment=1.,horizontal_spacing=0.06,font_name='Times',font_size=14)
;    if keyword_set(sav) then begin
;        print, '    SAVING PLOT'
;        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/xspec_model.eps',/BITMAP else $
;                                                     p.save,'figures/xspec_model.png',resolution=res
;    endif        

    
;;----------------------------------------------------------------------------------------
;; EXAMPLE OF NH MODELS WITH INCREASING FCT
;;----------------------------------------------------------------------------------------
if keyword_set(nhmodel) then begin

    nsamp = 10000.
    nh_samp = nh_mc(nh_obs,nsamp)
    ithin = where(nh_samp lt 24.,nthin)
    nct = round((nthin/(1.-dindgen(9,start=1)/10.))*(dindgen(9,start=1)/10.))     ;; determine the number of CT sources
    nh_fct10 = [nh_samp[ithin],24.+2.*randomu(seed,nct[0])]
    nh_fct20 = [nh_samp[ithin],24.+2.*randomu(seed,nct[1])]
    nh_fct30 = [nh_samp[ithin],24.+2.*randomu(seed,nct[2])]
    nh_fct40 = [nh_samp[ithin],24.+2.*randomu(seed,nct[3])]
    nh_fct50 = [nh_samp[ithin],24.+2.*randomu(seed,nct[4])]
    nh_fct60 = [nh_samp[ithin],24.+2.*randomu(seed,nct[5])]
    nh_fct70 = [nh_samp[ithin],24.+2.*randomu(seed,nct[6])]
    nh_fct80 = [nh_samp[ithin],24.+2.*randomu(seed,nct[7])]
    nh_fct90 = [nh_samp[ithin],24.+2.*randomu(seed,nct[8])]

    y10 = histogram(nh_fct10,bin=0.2,locations=x10,min=19.,max=27.)
    y20 = histogram(nh_fct20,bin=0.2,locations=x20,min=19.,max=27.)
    y30 = histogram(nh_fct30,bin=0.2,locations=x30,min=19.,max=27.)
    y40 = histogram(nh_fct40,bin=0.2,locations=x40,min=19.,max=27.)
    y50 = histogram(nh_fct50,bin=0.2,locations=x50,min=19.,max=27.)
    y60 = histogram(nh_fct60,bin=0.2,locations=x60,min=19.,max=27.)
    y70 = histogram(nh_fct70,bin=0.2,locations=x70,min=19.,max=27.)
    y80 = histogram(nh_fct80,bin=0.2,locations=x80,min=19.,max=27.)
    y90 = histogram(nh_fct90,bin=0.2,locations=x90,min=19.,max=27.)

    nrm = nthin*1.
    e = {xra:[18.9,27.1],yra:[0.,0.18], $
         stairstep:1,fill_background:1,color:'black', $
         xtitle:'$log !8N!7_H^{ model} / cm^2$',ytitle:'Frequency (arbitrary)', $
         font_name:'Times',font_size:14}
         
    ;; purples
;    col = [[63,0,125], $
;           [84,39,143], $
;           [106,81,163], $
;           [128,125,186], $
;           [158,154,200], $
;           [188,189,220], $
;           [218,218,235]]
    ;; blues
;    col = [[8,48,107], $
;           [8,81,156], $
;           [33,113,181], $
;           [66,146,198], $
;           [107,174,214], $
;           [158,202,225], $
;           [198,219,239], $
;           [222,235,247]]
    ;; reds
    col = [[127,0,0], $
           [179,0,0], $
           [215,48,31], $
           [239,101,72], $
           [252,141,89], $
           [253,187,132], $
           [253,212,158]]

    p = plot(x70,y70/nrm,_extra=e,fill_color=col[*,0])
    p = plot(x60,y60/nrm,_extra=e,fill_color=col[*,1],/ov)
    p = plot(x50,y50/nrm,_extra=e,fill_color=col[*,2],/ov)
    p = plot(x40,y40/nrm,_extra=e,fill_color=col[*,3],/ov)
    p = plot(x30,y30/nrm,_extra=e,fill_color=col[*,4],/ov)
    p = plot(x20,y20/nrm,_extra=e,fill_color=col[*,5],/ov)
    p = plot(x10,y10/nrm,_extra=e,fill_color=col[*,6],/ov)
    ;p = plot(x80,y80/nrm,/stairstep,/ov)
    ;p = plot(x90,y90/nrm,/stairstep,/ov)

    ar = ARROW([26.18,26.18], [0.05,0.13], COLOR=col[*,2], /DATA, /CURRENT,thick=2)
    ct = text(26.48,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ct = text(26.48,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ct = text(26.48,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ct10 = text(25.5,0.017,'0.10',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ct20 = text(25.5,0.033,'0.20',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ct30 = text(25.5,0.051,'0.30',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ct40 = text(25.5,0.075,'0.40',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ct50 = text(25.5,0.108,'0.50',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ct60 = text(25.5,0.158,'0.60',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/nh_model.eps',/BITMAP else $
                                                     p.save,'figures/nh_model.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; EXAMPLE OF RX MODELS WITH INCREASING FCT
;;----------------------------------------------------------------------------------------
if keyword_set(rxmodel) then begin

    nsamp = 10000.
    nh_samp = nh_mc(nh_obs,nsamp)
    ithin = where(nh_samp lt 24.,nthin)
    nct = round((nthin/(1.-dindgen(9,start=1)/10.))*(dindgen(9,start=1)/10.))     ;; determine the number of CT sources
    nh_fct10 = [nh_samp[ithin],24.+2.*randomu(seed,nct[0])]
    nh_fct20 = [nh_samp[ithin],24.+2.*randomu(seed,nct[1])]
    nh_fct30 = [nh_samp[ithin],24.+2.*randomu(seed,nct[2])]
    nh_fct40 = [nh_samp[ithin],24.+2.*randomu(seed,nct[3])]
    nh_fct50 = [nh_samp[ithin],24.+2.*randomu(seed,nct[4])]
    nh_fct60 = [nh_samp[ithin],24.+2.*randomu(seed,nct[5])]
    nh_fct70 = [nh_samp[ithin],24.+2.*randomu(seed,nct[6])]
    nh_fct80 = [nh_samp[ithin],24.+2.*randomu(seed,nct[7])]
    nh_fct90 = [nh_samp[ithin],24.+2.*randomu(seed,nct[8])]

    rx_fct10 = rx2nh(nh_fct10,/rx_out,scat=0.3)
    rx_fct20 = rx2nh(nh_fct20,/rx_out,scat=0.3)
    rx_fct30 = rx2nh(nh_fct30,/rx_out,scat=0.3)
    rx_fct40 = rx2nh(nh_fct40,/rx_out,scat=0.3)
    rx_fct50 = rx2nh(nh_fct50,/rx_out,scat=0.3)
    rx_fct60 = rx2nh(nh_fct60,/rx_out,scat=0.3)
    rx_fct70 = rx2nh(nh_fct70,/rx_out,scat=0.3)
    rx_fct80 = rx2nh(nh_fct80,/rx_out,scat=0.3)
    rx_fct90 = rx2nh(nh_fct90,/rx_out,scat=0.3)

    y10 = histogram(rx_fct10,bin=0.2,locations=x10,min=-4.4,max=1.4)
    y20 = histogram(rx_fct20,bin=0.2,locations=x20,min=-4.4,max=1.4)
    y30 = histogram(rx_fct30,bin=0.2,locations=x30,min=-4.4,max=1.4)
    y40 = histogram(rx_fct40,bin=0.2,locations=x40,min=-4.4,max=1.4)
    y50 = histogram(rx_fct50,bin=0.2,locations=x50,min=-4.4,max=1.4)
    y60 = histogram(rx_fct60,bin=0.2,locations=x60,min=-4.4,max=1.4)
    y70 = histogram(rx_fct70,bin=0.2,locations=x70,min=-4.4,max=1.4)
    y80 = histogram(rx_fct80,bin=0.2,locations=x80,min=-4.4,max=1.4)
    y90 = histogram(rx_fct90,bin=0.2,locations=x90,min=-4.4,max=1.4)

    e = {xra:[-4.2,1.2],yra:[0.,0.3], $
         stairstep:1,fill_background:1,color:'black', $
         xtitle:'$!8R_{L!7_X}^{ model}$',ytitle:'Frequency (arbitrary)', $
         font_name:'Times',font_size:14}
         
    ;; purples
;    col = [[63,0,125], $
;           [84,39,143], $
;           [106,81,163], $
;           [128,125,186], $
;           [158,154,200], $
;           [188,189,220], $
;           [218,218,235]]
    ;; blues
    col = [[8,48,107], $
           [8,81,156], $
           [33,113,181], $
           [66,146,198], $
           [107,174,214], $
           [158,202,225], $
           [198,219,239], $
           [222,235,247], $
           [247,251,255]]
    ;; reds
;    col = [[127,0,0], $
;           [179,0,0], $
;           [215,48,31], $
;           [239,101,72], $
;           [252,141,89], $
;           [253,187,132], $
;           [253,212,158]]
;    
    p = plot(x90,y90/10000.,_extra=e,fill_color=col[*,0])
    p = plot(x80,y80/10000.,_extra=e,fill_color=col[*,1],/ov)
    p = plot(x70,y70/10000.,_extra=e,fill_color=col[*,2],/ov)
    p = plot(x60,y60/10000.,_extra=e,fill_color=col[*,3],/ov)
    p = plot(x50,y50/10000.,_extra=e,fill_color=col[*,4],/ov)
    p = plot(x40,y40/10000.,_extra=e,fill_color=col[*,5],/ov)
    p = plot(x30,y30/10000.,_extra=e,fill_color=col[*,6],/ov)
    p = plot(x20,y20/10000.,_extra=e,fill_color=col[*,7],/ov)
    p = plot(x10,y10/10000.,_extra=e,fill_color=col[*,8],/ov)

    ar = ARROW([0.6,0.6], [0.08,0.22], COLOR=col[*,2], /DATA, /CURRENT,thick=2)
    ct = text(0.8,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ct = text(0.8,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ct = text(0.8,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ;ct10 = text(25.5,0.017,'0.10',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct20 = text(25.5,0.033,'0.20',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct30 = text(25.5,0.051,'0.30',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct40 = text(25.5,0.075,'0.40',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct50 = text(25.5,0.108,'0.50',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct60 = text(25.5,0.158,'0.60',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
stop
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/rx_model.eps',/BITMAP else $
                                                     p.save,'figures/rx_model.png',resolution=res
    endif
endif




;;----------------------------------------------------------------------------------------
;; RX-NH conversion
;;----------------------------------------------------------------------------------------
if keyword_set(rxnh) then begin

    restore,'/Users/ccarroll/Research/projects/model_nhdist/workspace/data_prep/rx_mc_borus.sav'    
    nhp = rebin(nh_fine,n_elements(nh_fine)*100d)
    rxgp = interpol(reform(rx_fine[*,0,-2]),nh_fine,nhp)
    rxhi = interpol(reform(rx_fine[*,0,-1]),nh_fine,nhp)
    rxlo = interpol(reform(rx_fine[*,0,0]),nh_fine,nhp)

    ;; index where upper and lower sigma intersect
    !null = min(abs(rxhi[100:-1]-rxlo[100:-1]),ind)
    ;; index where upper sigma intersects CT shading
    ict = value_locate(nhp,24.) 

    e = {xra:[21.0,25.0],yra:[-3.5,0.5], $
         thick:2,fill_background:1, $
         xtitle:'$log !8N!7_H / cm^2$',ytitle:'$!8R_{L_{!7X}}$', $
         font_name:'Times',font_size:14}

    ;; additive RGB purple
    col = [123,73,99]
    ; pretty purple
    col = [128,125,186]
    
    ;; CT shading
    pct = plot([24.,25.],e.yra[1]*[1.,1.],_extra=e,fill_color='light grey',fill_level=e.yra[0])
    ;; RX–NH shading
    plo = plot(nhp[0:ind],rxlo[0:ind],_extra=e,linestyle='',/ov,fill_color=col,fill_transparency=0,name='$\sigma_{scatt} = -2.0$')
    phi = plot(nhp[0:ind],rxhi[0:ind],_extra=e,linestyle='',/ov,fill_color='white',name='$\sigma_{scatt} = 0.5$')
    p = plot(nhp[ind:-1],rxlo[ind:-1],_extra=e,linestyle='',/ov,fill_color=col,fill_level=rxlo[ind],fill_transparency=0)
    p = plot(nhp[ind:-1],rxhi[ind:-1],_extra=e,linestyle='',/ov,fill_level=rxlo[ind],fill_color='white')
    ;; fix CT shading
    pct = plot(nhp[ict:-1],rxhi[ict:-1],_extra=e,linestyle='',fill_color='light grey',fill_level=rxlo[ind],/ov)
    ;; plot data
    plo = plot(nhp,rxlo,'-.',/ov,thick=2,name='$\sigma_{!8f!7_{scatt}} = -2.0$')
    phi = plot(nhp,rxhi,':',/ov,thick=2,name='$\sigma_{!8f!7_{scatt}} = 0.5$')
    pgp = plot(nhp,rxgp,'-',thick=2,/ov,name='$!8f!7_{scatt}$ (Gupta+2021)')
    ;; CT label
    t = text(24.10,0.15,'CT',col='black',font_size=14,font_name='Times',font_style='Bold',/data)
    ;; legend
    leg = legend(target=[pgp,phi,plo],sample_width=0.16,horizontal_spacing=0.06,font_name='Times',font_size=14)
    if (col[0] eq 123) then leg.position = [0.590,0.350] else $
                            leg.position = [0.540,0.350]

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/rx_nh.eps',/BITMAP else $
                                                     p.save,'figures/rx_nh.png',resolution=res
    endif


endif



;;----------------------------------------------------------------------------------------
;; X-RAY STACKED FLUXES
;;----------------------------------------------------------------------------------------

;; PER TONIMA
;; out_ is for  all sources
;; ind_ is for individual sources
;; cts means counts
;; ctss means counts in soft band
;; cths is counts in hard band
;; Ind_cts_w_bg_s means individual counts in soft band including background photons
;; ind_ctsb is background counts in soft band
;; ind_cthb is background in hard band
;; ind_net_count_s has background removed
;; ind_cnt_ps_s is individual counts per second in soft band
;; ind_softflux and is fluxes in soft band for each source
if keyword_set(stack) then begin
    ;; STACKFAST results for 91 non-detected Chandra sources
    dir = '/Users/ccarroll/Research/projects/model_nhdist/workspace/scratch/StackFast_7/'
    
    src = mrdfits(dir+'FINAL_stackdata_individual_columns.fits',1)
    stack = mrdfits(dir+'FINAL_stackdata.fits',1)
    
    bkg = [stack.out_ctsb,stack.out_cthb]
    net = [stack.out_net_ctss,stack.out_net_cths]
    cts = bkg+net
    cps = [stack.out_count_per_ss,stack.out_count_per_sh]
    ;; BACKGROUND CHECK
    ;; OUT_CTSB        DOUBLE           634.00000
    ;; OUT_CTHB        DOUBLE           1219.0000
    ;; IDL> print, total([[src.ind_ctsb],[src.ind_cthb]],1)
    ;;        634.00000       1219.0000    
    ;; NET CHECK
    ;; OUT_NET_CTSS    DOUBLE           106.00000
    ;; OUT_NET_CTHS    DOUBLE           109.00000
    ;; IDL> print, total([[src.ind_cts_w_bg_s],[src.ind_cts_w_bg_h]],1)
    ;;        106.00000       109.00000    
    texp = total(src.ind_exp_sky)
    ct2fx_cycle12_gamma14 = [3.612E-12,7.439E-12]
    ct2fx_cycle12_gamma18 = [4.356E-12,5.345E-12]
    ct2fx_current_gamma14 = [5.936E-12,1.223E-11]
    ct2fx_current_gamma18 = [8.141E-12,2.016E-11]

    fxs  = net/texp*ct2fx_cycle12_gamma14    
    fxs2 = net/texp*ct2fx_cycle12_gamma18
    fxs3 = net/texp*ct2fx_current_gamma14
    fxs4 = net/texp*ct2fx_current_gamma18
    fxs_all = [[fxs],[fxs2]];,[fxs3],[fxs4]]
    unc_fxs = stddev(fxs_all,dim=2)/mean(fxs_all,dim=2)
    
    
    err = [stack.out_errcntpss/stack.out_count_per_ss,stack.out_errcntpsh/stack.out_count_per_sh]
    unc_phot = 0.011/0.550  ;; model
    unc_vari = 0.016/0.584  ;; model
    unc_simu = 0.018/0.550  ;; model
    unc_c21x = [1.180/3.990,0.360/1.320]  ;; stack
    
    e_fxs = (sqrt(err^2. + unc_fxs^2. + unc_c21x^2.)) * fxs
    stop
endif



;;----------------------------------------------------------------------------------------
;; MODELED NH DISTRUBITUION w/ STACKING RESULTS --- ANANNA+2019 AVERAGED NH
;;----------------------------------------------------------------------------------------
if keyword_set(newnh) then begin
    ;; X-ray stack image path
    prep_dir = file_search('data_prep')
    up = ''
    while (prep_dir eq '') do begin
        up = up+'../'
        prep_dir = file_search(up+'data_prep')
    endwhile
    soft = file_search(prep_dir+'/nodet_wagn_052_stack.png')
    hard = file_Search(prep_dir+'/nodet_wagn_27_stack.png')

    ;; STACKFAST FLUXES
    energy_center = [1.25,4.5]
    energy_range = [0.75,2.5]
    modoff = [[0.76,0.77],[1.31,1.32]]
    ;fxstak[1,0:1]
    ;stack_flux = [4.2155345e-16,8.9311185e-16]  ;; Chandra current cycle gamma 1.4
    ;stack_flux = [5.7814466e-16,1.4722106e-15]  ;; Chandra current cycle gamma 1.8
    stack_flux = [2.5651130e-16,5.4324276e-16]  ;; Chandra cycle 12 gamma 1.4 (previous step FXS)
    ;stack_flux = [3.0934751e-16,3.9032567e-16]    ;; Chandra cycle 12 gamma 1.8
    stack_err = [8.9192778e-17,2.1103027e-16]   ;; additional uncertainties based on photon index
    ;stack_err = [8.2510708e-17,1.6940501e-16]   ;; poisson error from STACKFAST (previous step E_FXS)
    ;stack_err = [1.24e-16,4.00e-16];[7.3145295e-17,2.2261186e-16];e_fxstak[1,0:1]
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    ;; PARAMETER AND ASSUMPTION UNCERTAINTIES
    unc_phot = 0.011/0.550  ;; model
    unc_vari = 0.016/0.584  ;; model
    unc_simu = 0.018/0.550  ;; model
    unc_lxir = 0.3
    unc_c21s = 1.180/3.990  ;; stack
    unc_c21h = 0.360/1.320  ;; stack

    ;; MODELING FLUXES
    ;; uniform
    ;; soft: 1.825227808798384e-16 \pm 6.617095252140283e-18
    ;; hard: 1.4501219449132314e-15 \pm 2.4950863035050844e-16
    ;; variable 
    ;; soft: 1.4777313317360868e-16 \pm 6.408453646783788e-18
    ;; hard: 1.1265466583256214e-15 \pm 2.399513231681352e-16
    sfx = [1.825227808798384e-16,1.4777313317360868e-16]
    hfx = [1.4501219449132314e-15,1.1265466583256214e-15]
    e_sfx = [6.617095252140283e-18,6.408453646783788e-18]
    e_hfx = [2.4950863035050844e-16,2.399513231681352e-16]
    fx_mod = [mean(sfx),mean(hfx)]
    e_fx_mod = sqrt((e_sfx/sfx)^2. + (e_hfx/hfx)^2. + unc_phot^2. + unc_vari^2. + unc_simu^2. + unc_lxir^2.) * fx_mod
    
    ;; NH MODELING UNCERTAINTIES
    err = nh_mod.sig/nh_mod.yh > 0.
    e_nhmod = sqrt(err^2. + unc_phot^2. + unc_vari^2. + unc_simu^2.) * nh_mod.yh

    ;; variables for FIXED vs FREE
    fx_models = ['FX_NON','FX_NON_']
    nh_models = ['NH_MOD','NH_MOD_']
    ctf_val = ['FCT','FCT_']
    
    ;; plot set up
    dim = [880,780]
    sq = 180
    gap = 40
    pos = [[100,560,100+sq,560+sq],[100+sq+gap,560,100+2.*sq+gap,560+sq],[600,560,840,560+sq],$
                                   [100,70,840,480]]
    e = {font_name:'Times', $
         dim:dim,device:1,buffer:0}
    enh = {xra:[20.,26.],yra:[0.,1.0], $
           stairstep:1, $
           xtitle:'$log !8N!7_{H}  [cm^{-2}]$',ytitle:'Frequency (normalized; $!8N!7_H < 10^{24} cm^{-2}$)', $;ytitle:'Frequency (normalized)', $
           font_name:'Times',font_size:14, $
           current:1,dim:dim,device:1}
    if keyword_set(hide) then e.buffer = 1

    ;; color
    ;col = [106,81,163];[84,39,143]
            ;; purple     [84,39,143]
            ;; teal       [28,144,153]
            ;; mint       [127,205,187]
            ;; aquamarine [28,144,153]
    ;; reds
    col = [[127,0,0], $
           [179,0,0], $
           [215,48,31], $
           [239,101,72], $
           [252,141,89], $
           [253,187,132], $
           [253,212,158]]

    ;; BEGIN FIGURE
    ;; plot STACKFAST FLUXES
    psoft = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,0],_extra=e,/device,ytitle='offset in Decl. [arcsec.]',font_size=14)
    isoft = image(soft,pos=pos[*,0],/current,/device)    
    tsoft = text(target=psoft,0.5,0.82,'  0.5–2 keV  ',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    phard = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,1],/current,/device);,xtitle='offset in R.A. [arcsec.]                     ')
    ihard = image(hard,pos=pos[*,1],/current,/device)
    thard = text(target=phard,0.5,0.82,'  2–7 keV  ',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
    xt = text((pos[2,0]+pos[0,1])/2.,pos[1,0]-42,/device,'offset in R.A. [arcsec.]',alignment=0.5,font_name='Times',font_size=14)

    ;; plot X-ray flux estimates
    pxd = errorplot(energy_center,stack_flux,energy_range,stack_err, $
                   linestyle='',/xlog,/ylog,pos=pos[*,2],/current,/device, $
                   xtitle='$energy  [keV]$',ytitle='$!8F!7_{X}  [erg s^{-1} cm^{-2}]$', $
                   font_name='Times',font_size=14)
    pxd = plot(energy_center,stack_flux,'s',col='black', $
               sym_size=1.5,sym_thick=1.5,sym_filled=1,sym_fill_color='white',/ov,name='Stack')
    pxd.xra=[1e-1,1e1]
    pxd.yra=[5e-17,1e-14]
    ;; plot MODEL FLUXES
    fxm = fx_mod
    e_fxm = e_fx_mod
    ;;; UNIFORM
    ;pxm = errorplot(energy_center*modoff[*,0],fxm[*,0],e_fxm[*,0], $
    ;                    linestyle='',errorbar_capsize=0.1,/ov)
    ;pxm = plot(energy_center*modoff[*,0],fxm[*,0],'S',col='black', $
    ;               sym_size=2.0,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    ;pxm = plot(energy_center*modoff[*,0],fxm[*,0],'S',col='black', $
    ;               sym_size=2.0,sym_thick=1.5,sym_filled=1,sym_fill_color=col[*,2],sym_transparency=20,fill_transparency=20,/ov)
    ;;; VARIABLE
    ;pxm2 = errorplot(energy_center*modoff[*,1],fxm[*,1],e_fxm[*,1], $
    ;                     linestyle='',errorbar_capsize=0.1,/ov)
    ;pxm2 = plot(energy_center*modoff[*,1],fxm[*,1],'S',col='black', $
    ;                sym_size=2.0,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    ;pxm2 = plot(energy_center*modoff[*,1],fxm[*,1],'S',col='black', $
    ;                sym_size=2.0,sym_thick=1.5,sym_filled=1,sym_fill_color=col[*,2],sym_transparency=20,fill_transparency=20,/ov)
    ;; AVERAGED
    pxm = errorplot(energy_center*modoff[*,0],fxm[*,0],e_fxm[*,0], $
                        linestyle='',errorbar_capsize=0.1,/ov)
    pxm = plot(energy_center*modoff[*,0],fxm[*,0],'S',col='black', $
                   sym_size=2.0,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    pxm = plot(energy_center*modoff[*,0],fxm[*,0],'S',col='black', $
                   sym_size=2.0,sym_thick=1.5,sym_filled=1,sym_fill_color=col[*,2],sym_transparency=20,fill_transparency=20,/ov)
    
    ps = plot([1.5e-1],[5.30e-15],'s',col='black',sym_size=1.5,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov) 
    ts = text(0.15,0.84,target=pxd,/relative,'X-ray stacking',font_name='Times',font_size=14)
    po = plot([1.5e-1],[2.80e-15],'S',col='black',sym_size=2.0,sym_thick=1.5,sym_filled=1,sym_fill_color=col[*,2],sym_transparency=10,/ov)
    to = text(0.15,0.72,target=pxm,/relative,'Model fluxes',font_name='Times',font_size=14)
    
    ;; plot NH DISTRIBUTION
    pct = plot([24.,26.],[1.,1.]*enh.yra[1],_extra=enh,pos=pos[*,3], $
               fill_background=1,fill_color='light grey',fill_level=0.,fill_transparency=0)
    pct.stairstep = 0   ;; CT shading
    
    pnh = errorplot(nh_mod.xh+0.5,nh_mod.yh/total(nh_mod[where(nh_mod.xh lt 24.)].yh),e_nhmod,'-',thick=4,errorbar_thick=4,_extra=enh, $
                              pos=pos[*,3],fill_background=1,fill_color=col[*,2],fill_transparency=20,name='This work')
    
    ;pnh2 = errorplot(nh_mod_.xh+0.5,nh_mod_.yh/total(nh_mod_[where(nh_mod_.xh lt 24.)].yh),e_nhmod,':',thick=4,errorbar_thick=4,_extra=enh, $
    ;                           pos=pos[*,3],fill_background=0,fill_color=col[*,2],fill_transparency=20,name='This work again!')

    ;; Ananna+2019
    pan = plot(nh_ana_lox.xh,nh_ana_lox.yh*frac_lox+nh_ana_hix.yh*frac_hix,'__',thick=4,_extra=enh,/ov,name='$Ananna+2019 (weighted)$')
    tct = text(0.683,0.91,'CT',target=pct,/relative,font_name='Times',font_size=14)
    ;p = plot(nh_ric_int.xh,nh_ric_int.yh/total(nh_ric_int[where(nh_ric_int.xh lt 24.)].yh),':',/ov,/stairstep,thick=4)
    ;; CT Fraction text and legend
    ctad = text(25.,enh.yra[1]/3.,'$!8f!7_{CT} = 0.562^{+0.066}_{-0.090}$',target=pnh,/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5)
    leg = legend(target=[pnh,pan],position=[0.138,0.59],/normal,horizontal_alignment=0.,sample_width=0.14,font_size=14,font_name='Times')
    
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then pct.save,'figures/nh_dist.eps',/BITMAP else $
                                                     pct.save,'figures/nh_dist.png',resolution=res
    endif    
endif






if keyword_set(new_rxnh) then begin

    rx = mrdfits('rx.fits',0,hd)
    nh = mrdfits('rx.fits',3)
    phot = mrdfits('rx.fits',4)
    scat = mrdfits('rx.fits',5)
    refl = mrdfits('rx.fits',6)
    open = mrdfits('rx.fits',7)

    ;; LEFT SIDE, 6 CHANNELS    
    nrx = rx[7,*,*,*,*]
    ind_left = []
    ;; 1
    il1 = where(nrx lt -0.75 and nrx gt -0.85,nl1)
    indl1 = array_indices(rx,il1)
    mm =  minmax( nrx[il1] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_left = [[ind_left],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; uniform photon index
    ;; uniform scattering
    ;; uniform reflection 0>=R>=1
    ;; 10° opening angle only
    ;; 2
    il2 = where(nrx lt -0.6 and nrx gt -0.75,nl2)
    indl2 = array_indices(rx,il2)
    mm = minmax( nrx[il2] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_left = [[ind_left],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 3
    il3 = where(nrx lt -0.45 and nrx gt -0.6,nl3)
    indl3 = array_indices(rx,il3)
    mm = minmax( nrx[il3] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_left = [[ind_left],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 4
    il4 = where(nrx lt -0.35 and nrx gt -0.45,nl4)
    indl4 = array_indices(rx,il4)
    mm = minmax( nrx[il4] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_left = [[ind_left],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 5
    il5 = where(nrx lt -0.24 and nrx gt -0.35,nl5)
    indl5 = array_indices(rx,il5)
    mm = minmax( nrx[il5] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_left = [[ind_left],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 6
    il6 = where(nrx lt -0.15 and nrx gt -0.24,nl6)
    indl6 = array_indices(rx,il6)
    mm = minmax( nrx[il6] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_left = [[ind_left],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; COMBINE
    nl = nl1+nl2+nl3+nl4+nl5+nl6
    rx_left = []
    for i = 0,n_elements(ind_left[0,*])-1 do begin
        str = '*,'+strjoin(strtrim(ind_left[1:-1,i],2),',')
        re = execute('rx_left = [[rx_left],[rx['+str+']]]')
    endfor


    ;; RIGHT SIDE, 9 CHANNELS
    nrx = rx[16,*,*,*,*]
    ind_right = []
    ;; 1
    ir1 = where(nrx lt -2.8 and nrx gt -3.0,nr1)
    indr1 = array_indices(rx,ir1)
    mm =  minmax( nrx[ir1] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 
    ;; 
    ;; 
    ;; only 10° opening angle
    ;; 2
    ir2 = where(nrx lt -2.52 and nrx gt -2.7,nr2)
    indr2 = array_indices(rx,ir2)
    mm = minmax( nrx[ir2] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 3
    ir3 = where(nrx lt -2.4 and nrx gt -2.52,nr3)
    indr3 = array_indices(rx,ir3)
    mm = minmax( nrx[ir3] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 4
    ir4 = where(nrx lt -2.25 and nrx gt -2.4,nr4)
    indr4 = array_indices(rx,ir4)
    mm = minmax( nrx[ir4] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 5
    ir5 = where(nrx lt -2.15 and nrx gt -2.25,nr5)
    indr5 = array_indices(rx,ir5)
    mm = minmax( nrx[ir5] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 6
    ir6 = where(nrx lt -2.1 and nrx gt -2.15,nr6)
    indr6 = array_indices(rx,ir6)
    mm = minmax( nrx[ir6] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 7
    ir7 = where(nrx lt -1.95 and nrx gt -2.1,nr7)
    indr7 = array_indices(rx,ir7)
    mm = minmax( nrx[ir7] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 8
    ir8 = where(nrx lt -1.91 and nrx gt -1.95,nr8)
    indr8 = array_indices(rx,ir8)
    mm = minmax( nrx[ir8] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; 9
    ir9 = where(nrx lt -1.0 and nrx gt -1.91,nr9)
    indr9 = array_indices(rx,ir9)
    mm = minmax( nrx[ir9] )
    imin = where(rx eq mm[0])
    imax = where(rx eq mm[1])
    ind_right = [[ind_right],[array_indices(rx,imin)],[array_indices(rx,imax)]]
    ;; COMBINE
    nr = nr1+nr2+nr3+nr4+nr5+nr6+nr7+nr8+nr9
    rx_right = []
    for i = 0,n_elements(ind_right[0,*])-1 do begin
        str = '*,'+strjoin(strtrim(ind_right[1:-1,i],2),',')
        re = execute('rx_right = [[rx_right],[rx['+str+']]]')
    endfor






endif












END










