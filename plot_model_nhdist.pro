PRO plot_model_nhdist, PROPERTIES = properties, $
                       NHDIST = nhdist, $
                       RXDIST = rxdist, $
                       XSPEC = xspec, $
                       NHMODEL = nhmodel, $
                       RXMODEL = rxmodel, $
                       RXNH = rxnh, $
                       STACK = stack, $
                       NEW_NHDIST = new_nhdist, $
                       PARAM_COMP = param_comp, $
                       MRLX = mrlx, $
                       CORNER = corner, $
                       HIDE = hide, $
                       LOW_RES = low_res, $
                       SAV = sav


;common _data
;common _nhdist
;common _nhobs
;common _rxnh
;common _group
;common _uniform
;;common _variable
;common _model


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
    file = '/Users/ccarroll/Research/projects/model_nhdist/workspace/data_prep/mc_modelcomps.txt'

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

    e = {xra:[0.1,300],yra:[5e-5,5e3], $
         xlog:1,ylog:1, $
         thick:3, $
         xtickname:['0.1','1','10','100'],ytickname:'10!U'+strtrim([-4:3],2), $
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

    readcol,'/Users/ccarroll/Research/projects/model_nhdist/workspace/data_prep/Ricci+2017_Fig23_data.csv', $
            nh,freq,format='d,d'
    nh_obs = soa2aos({xh:nh,yh:freq})
    nh_obs[where(nh_obs.xh eq 20.5,/null)].yh = 0.0
    
    nsamp = 50000.
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

    ;ar = ARROW([26.18,26.18], [0.05,0.13], COLOR=col[*,2], /DATA, /CURRENT,thick=2)
    ar = ARROW([0.835,0.835], [0.35,0.65], COLOR=col[*,2], /NORMAL, /CURRENT,thick=2)
    ct = text(0.85,0.5,'$!8f!7_{CT}^{  model}$',/NORMAL,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=1.0,orientation=90)
    ;ct = text(26.48,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ;ct = text(26.48,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
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

    readcol,'/Users/ccarroll/Research/projects/model_nhdist/workspace/data_prep/Ricci+2017_Fig23_data.csv', $
            nhric,freq,format='d,d'
    nh_obs = soa2aos({xh:nhric,yh:freq})
    nh_obs[where(nh_obs.xh eq 20.5,/null)].yh = 0.0

    nsamp = 50000.
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

    e = {xra:[-4.2,1.2],yra:[0.,0.25], $
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
    p = plot(x90,y90/nsamp,_extra=e,fill_color=col[*,0])
    p = plot(x80,y80/nsamp,_extra=e,fill_color=col[*,1],/ov)
    p = plot(x70,y70/nsamp,_extra=e,fill_color=col[*,2],/ov)
    p = plot(x60,y60/nsamp,_extra=e,fill_color=col[*,3],/ov)
    p = plot(x50,y50/nsamp,_extra=e,fill_color=col[*,4],/ov)
    p = plot(x40,y40/nsamp,_extra=e,fill_color=col[*,5],/ov)
    p = plot(x30,y30/nsamp,_extra=e,fill_color=col[*,6],/ov)
    p = plot(x20,y20/nsamp,_extra=e,fill_color=col[*,7],/ov)
    p = plot(x10,y10/nsamp,_extra=e,fill_color=col[*,8],/ov)



    ar = ARROW([0.235,0.235], [0.35,0.65], COLOR=col[*,2], /NORMAL, /CURRENT,thick=2)
    ct = text(0.21,0.5,'$!8f!7_{CT}^{  model}$',/NORMAL,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.,orientation=90)
    ;ct = text(0.8,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ;ct = text(0.8,e.yra[1]/2.,'$!8f!7_{CT}^{  model}$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5,orientation=90)
    ;ct10 = text(25.5,0.017,'0.10',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct20 = text(25.5,0.033,'0.20',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct30 = text(25.5,0.051,'0.30',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct40 = text(25.5,0.075,'0.40',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct50 = text(25.5,0.108,'0.50',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)
    ;ct60 = text(25.5,0.158,'0.60',/data,font_size=14,font_name='Times',font_style='Bold',alignment=0.5,vertical_alignment=0.5)

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

    ;; read data
    rx = mrdfits('rx.fits',0,hd)
    nh = mrdfits('rx.fits',3)
    phot = mrdfits('rx.fits',4)
    scat = mrdfits('rx.fits',5)
    refl = mrdfits('rx.fits',6)
    open = mrdfits('rx.fits',7)
    
    ;; rebin for finer elements so the plot lines match up when shading
    nhp = rebin(nh,n_elements(nh)*10d)
    sz = size(rx,/dim)
    nsrc = product(sz[1:-1])
    vec_rx = reform(rx,sz[0],nsrc)
    rxp = dblarr(n_elements(nhp),nsrc)
    for i = 0,nsrc-1 do rxp[*,i] = interpol(vec_rx[*,i],nh,nhp)
    
    ;; find the upper and lower bounds
    ;for i = 0,n_elements(nhp)-1 do begin
    ;   mm = minmax(rxp[i,*],imm)
    ;   print, nhp[i], mm[0], imm[0], mm[1], imm[1]
    ;endfor
    ;% Compiled module: MINMAX.
    ;   21.000000       0.0000000           0       0.0000000           0
    ;   21.250000     -0.11406621        1205    -0.012380086        1566
    ;   21.500000     -0.22716184        1205    -0.025157160        1566
    ;   21.750000     -0.33922480        1205    -0.040565469        1566
    ;   22.000000     -0.45035115        1403    -0.058008511         180
    ;   22.250000     -0.56108265        1403    -0.087032824         180
    ;   22.500000     -0.67243711        1403     -0.13039931         180
    ;   22.750000     -0.78581145        1403     -0.19339732         180
    ;   23.000000     -0.90249016        1403     -0.27963530         180
    ;   23.250000      -1.0236432        1403     -0.39574794         180
    ;   23.500000      -1.1504725        1403     -0.55835479         180
    ;   23.750000      -1.2827198        1403     -0.80429493         180
    ;   24.000000      -1.4360859        1388      -1.1948127          18
    ;   24.250000      -2.0294886        1370      -1.5392724          15
    ;   24.500000      -2.5777558        1370      -1.6623632        1401
    ;   24.750000      -2.7985731        1370      -1.7803468        1401
    ;   25.000000      -2.9315420        1370      -1.8980953        1401    
    
    
    ;; bounds for left and right of intersection
    ileft = [180,1403]
    iright = [1401,1370]
    ;; index where upper and lower sigma intersect
    ioff = 25
    !null = min(abs(rxp[ioff:-1,ileft[0]]-rxp[ioff:-1,ileft[1]]),ind)
    inters_left = ind+ioff
    !null = min(abs(rxp[ioff:-1,iright[0]]-rxp[ioff:-1,iright[1]]),ind)
    inters_right = ind+ioff

    ;; index for CT shading
    ict = value_locate(nhp,24.) 

    e = {xra:[21.0,25.0],yra:[-3.5,0.5], $
         thick:1,fill_background:1,fill_level:-3.5, $
         xtitle:'$log !8N!7_H / cm^2$',ytitle:'$!8R_{L_{!7X}}$', $
         font_name:'Times',font_size:14}
    transp = 40
    ;; similar red as NH_DIST.PNG
    ;col = [223,89,75]
    col = [215,48,31]

    ;; CT shading
    pct = plot([24.,25.],e.yra[1]*[1.,1.],_extra=e,fill_color='light grey')
    ;; remove CT shading in plot area for transparency to look right
    pct = plot(nhp[ict:-1],rxp[inters_right:-1,iright[0]],_extra=e,color='white',fill_color='white',/ov)
    pct = plot(nhp[ict:inters_left],rxp[ict:inters_left,ileft[0]],_extra=e,col='white',fill_col='white',/ov)
    ;; RX–NH shading
    p = plot(nhp[0:inters_left],rxp[0:inters_left,ileft[0]],_extra=e,col=col,fill_color=col,transparency=transp,fill_transparency=transp,/ov)
    p = plot(nhp[0:inters_left],rxp[0:inters_left,ileft[1]],_extra=e,col=col,transparency=transp,fill_color='white',/ov)
    ;; fix CT shading
    ;pct = plot(nhp[ict:inters_left],rxp[ict:inters_left,ileft[1]],_extra=e,linestyle='',fill_color='light grey',/ov)
    ;; continue RX-NH shading
    p = plot(nhp[inters_right:-1],rxp[inters_right:-1,iright[0]],_extra=e,col=col,fill_color=col,transparency=transp,fill_transparency=transp,/ov)
    p = plot(nhp[inters_right:-1],rxp[inters_right:-1,iright[1]],_extra=e,col=col,transparency=transp,fill_color='white',/ov)
    ;; fix CT shading
    pct = plot(nhp[ict:-1],rxp[ict:-1,iright[1]],_extra=e,linestyle='',fill_color='light grey',/ov)
    ;; plot data
    for i = 0,nsrc-1 do p = plot(nh,vec_rx[*,i],col=col,transparency=98,/ov)

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
if keyword_set(new_nhdist) then begin
    restore,'../observed_nhdist.sav'
    restore,'select_nhobs.sav'
    restore,'select_group.sav'
    restore,'rx_model.sav'

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
    stack_flux = [5.7814466e-16,1.4722106e-15]  ;; Chandra current cycle gamma 1.8
    ;stack_flux = [2.5651130e-16,5.4324276e-16]  ;; Chandra cycle 12 gamma 1.4 (previous step FXS)
    stack_err = [3.0934751e-16,5.9032567e-16]    ;; Chandra cycle 12 gamma 1.8
    ;stack_err = [8.9192778e-17,2.1103027e-16]   ;; additional uncertainties based on photon index
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
    ;sfx = [1.825227808798384e-16,1.4777313317360868e-16]
    ;hfx = [1.4501219449132314e-15,1.1265466583256214e-15]
    ;e_sfx = [6.617095252140283e-18,6.408453646783788e-18]
    ;e_hfx = [2.4950863035050844e-16,2.399513231681352e-16]
    ;fx_mod = [mean(sfx),mean(hfx)]
    ;e_fx_mod = sqrt((e_sfx/sfx)^2. + (e_hfx/hfx)^2. + unc_phot^2. + unc_vari^2. + unc_simu^2. + unc_lxir^2.) * fx_mod
    
    ;; MCMC
    ;fx_mod = [mean(fx_non.csoft),mean(fx_non.chard)]
    ;e_fx_mod = [mean(fx_non.e_csoft),mean(fx_non.e_chard)]
    fx_mod = [mean(fx_non.csoft),mean(fx_non.chard)]
    e_fx_mod = [mean(fx_non.e_csoft),mean(fx_non.e_chard)]
    
    ;; NH MODELING UNCERTAINTIES
    ;err = nh_mod.sig/nh_mod.yh > 0.
    ;e_nhmod = sqrt(err^2. + unc_phot^2. + unc_vari^2. + unc_simu^2.) * nh_mod.yh
    
    print, "LOG (FX MODEL - FX STACK) [dex]"
    print, alog10(fx_mod/stack_flux)
    
    ;; MCMC
    nhm = nh_mod
    ;nhm = nh_mod_
    err = nhm.sig/nhm.yh > 0.
    ;; uncertainties on NH MODEL, FCT, FSCAT
    ;e_nhm = sqrt(err^2. + (0.037/0.555)^2. + (0.13/1.04)^2. )
    ;print, e_nhm
    e_nhm = sqrt(err*0. + (0.037/0.555)^2. + (0.13/1.04)^2. )
    ;print, e_nhm
    
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
                   linestyle='',errorbar_thick=2,/xlog,/ylog,pos=pos[*,2],/current,/device, $
                   xtitle='$energy  [keV]$',ytitle='$!8F!7_{X}  [erg s^{-1} cm^{-2}]$', $
                   font_name='Times',font_size=14)
    pxd = plot(energy_center,stack_flux,'s',col='black', $
               sym_size=1.5,sym_thick=2,sym_filled=1,sym_fill_color='white',/ov,name='Stack')
    pxd.xra=[1e-1,1e1]
    pxd.yra=[5e-17,3e-14]
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
                        linestyle='',errorbar_thick=2,errorbar_capsize=0.1,/ov)
    pxm = plot(energy_center*modoff[*,0],fxm[*,0],'S',col='black', $
                   sym_size=2.0,sym_thick=2,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
    pxm = plot(energy_center*modoff[*,0],fxm[*,0],'S',col='black', $
                   sym_size=2.0,sym_thick=2,sym_filled=1,sym_fill_color=col[*,2],sym_transparency=20,fill_transparency=20,/ov)
    
    ps = plot([1.5e-1],[1.40e-14],'s',col='black',sym_size=1.5,sym_thick=2,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov) 
    ts = text(0.15,0.84,target=pxd,/relative,'X-ray stacking',font_name='Times',font_size=14)
    po = plot([1.5e-1],[6.50e-15],'S',col='black',sym_size=2.0,sym_thick=2,sym_filled=1,sym_fill_color=col[*,2],sym_transparency=10,/ov)
    to = text(0.15,0.72,target=pxm,/relative,'Model fluxes',font_name='Times',font_size=14)
    
    ;; plot NH DISTRIBUTION
    pct = plot([24.,26.],[1.,1.]*enh.yra[1],_extra=enh,pos=pos[*,3], $
               fill_background=1,fill_color='light grey',fill_level=0.63,fill_transparency=0)
    pct.stairstep = 0   ;; CT shading
    
    pnh = errorplot(nhm.xh+0.5,nhm.yh/total(nhm[where(nhm.xh lt 24.)].yh),e_nhm,'-',thick=4,errorbar_thick=4,_extra=enh, $
                              pos=pos[*,3],fill_background=1,fill_color=col[*,2],fill_transparency=20,name='This work')
    
    ;pnh2 = errorplot(nh_mod_.xh+0.5,nh_mod_.yh/total(nh_mod_[where(nh_mod_.xh lt 24.)].yh),e_nhmod,':',thick=4,errorbar_thick=4,_extra=enh, $
    ;                           pos=pos[*,3],fill_background=0,fill_color=col[*,2],fill_transparency=20,name='This work again!')

    ;; Ananna+2019
    pan = plot(nh_ana_lox.xh,nh_ana_lox.yh*frac_lox+nh_ana_hix.yh*frac_hix,'__',thick=4,_extra=enh,/ov,name='$Ananna+2019 (weighted)$')
    tct = text(0.683,0.91,'CT',target=pct,/relative,font_name='Times',font_size=14)
    ;p = plot(nh_ric_int.xh,nh_ric_int.yh/total(nh_ric_int[where(nh_ric_int.xh lt 24.)].yh),':',/ov,/stairstep,thick=4)
    ;; CT Fraction text and legend
    ctad = text(25.,enh.yra[1]/3.5,'$!8f!7_{CT} = 0.555^{+0.037}_{-0.032}$',target=pnh,/data,font_size=18,font_name='Times',alignment=0.5,vertical_alignment=0.5)
    leg = legend(target=[pnh,pan],position=[0.138,0.59],/normal,horizontal_alignment=0.,sample_width=0.14,font_size=14,font_name='Times')
    
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then pct.save,'figures/nh_dist.eps',/BITMAP else $
                                                     pct.save,'figures/nh_dist.png',resolution=res
    endif    
endif



;;----------------------------------------------------------------------------------------
;; FCT vs F24 and FSCATT vs F24
;;----------------------------------------------------------------------------------------
if keyword_set(param_comp) then begin

    ;; PLOT FINAL PARAM SET
    f24 = [0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    f25 = 1.-f24

    fct = [0.552,0.550,0.565,0.555,0.559,0.559,0.551]
    fct_upp = [0.035,0.038,0.026,0.037,0.035,0.029,0.040]
    fct_low = [0.033,0.030,0.041,0.032,0.034,0.038,0.030]
    fct_err = transpose([[fct_low],[fct_upp]])

    fscat = [-0.91,-0.96,-0.99,-1.04,-1.08,-1.14,-1.14]
    fscat_upp = [0.09,0.11,0.10,0.13,0.14,0.14,0.12]
    fscat_low = [0.14,0.14,0.15,0.12,0.12,0.11,0.12]
    fscat_err = transpose([[fscat_low],[fscat_upp]])

    ;; FTC vs F24
;   efct = {xra:[0.1,0.9],yra:[0.41,0.69], $
;           ;xtitle:'$!8f!7_{24-25}$',ytitle:'$!8f!7_{CT}$', $
;           xtitle:'$!8f!7_{25-26}$',ytitle:'$!8f!7_{CT}$', $
;           font_name:'Times',font_size:14}
;
;   pfct = errorplot(f25,fct,fct_err,linestyle='',_extra=efct,errorbar_thick=2)
;   p = plot(f25,fct,'o',_extra=efct,col='black',sym_filled=1,sym_fill_color=[215,48,39],sym_size=1.5,sym_thick=2,/ov)
;   
;   ;; save image
;   if keyword_set(sav) then begin
;       print, '    SAVING PLOT'
;       if (strupcase(strtrim(sav,2)) eq 'EPS') then pfct.save,'figures/fct_vs_f25.eps',/BITMAP else $
;                                                    pfct.save,'figures/fct_vs_f25.png',resolution=res
;   endif    
;
;   ;; FSCAT vs F24
;   efscat = {xra:[0.1,0.9],yra:[-1.5,-0.5], $
;             sym_filled:1, $
;             ;xtitle:'$!8f!7_{24-25}$',ytitle:'$\sigma_{scatt}$', $
;             xtitle:'$!8f!7_{25-26}$',ytitle:'$\sigma_{scatt}$', $
;             font_name:'Times',font_size:14}
;
;   pfscat = errorplot(f25,fscat,fscat_err,'o',_extra=efscat,errorbar_thick=2)
;   p = plot(f25,fscat,'o',_extra=efscat,col='black',sym_filled=1,sym_fill_color=[69,117,180],sym_size=1.5,sym_thick=2,/ov)
;   
;   ;; save image
;   if keyword_set(sav) then begin
;       print, '    SAVING PLOT'
;       if (strupcase(strtrim(sav,2)) eq 'EPS') then pfscat.save,'figures/fscat_vs_f25.eps',/BITMAP else $
;                                                    pfscat.save,'figures/fscat_vs_f25.png',resolution=res
;   endif    

    ;; PLOT ADDITIONAL PARAMETER TESTS
    par = read_csv('all_params_in_table_format.csv',header=hd)
    for i = 0,n_elements(hd)-1 do re = execute(hd[i]+' = par.field'+string(i+1,'(i02)'))
    ;; f24 is read in here and overwrites the previous entry. need to create f25 again
    f25 = 1.-f24
    
    iio30 = oa eq 30.d
    iio60 = oa eq 60.d
    iir05 = refl eq 0.5d
    iir10 = refl eq 1.0d
    iig18 = gam eq 1.8d
    iig19 = gam eq 1.9d
    iig20 = gam eq 2.0d
    
    fct_err = abs(transpose([[fct_low],[fct_upp]]))
    fscat_err = abs(transpose([[fscat_low],[fscat_upp]]))
    
    ;; FTC vs F24 BY PARAM
    efct = {xra:[0.1,0.9],yra:[0.46,0.74], $
            sym_color:'black', $;sym_filled:1,sym_size:1.2, $
            ;xtitle:'$!8f!7_{24-25}$',ytitle:'$!8f!7_{CT}$', $
            xtitle:'$!8f!7_{25-26}$',ytitle:'$!8f!7_{CT}$', $
            font_name:'Times',font_size:14}

    ;; OPENING ANGLE 30°
    ii = where(iio30 and iir05 and iig18)
    poa30 = errorplot(f25[ii]-0.025,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',errorbar_thick=2)
    p1 = plot(f25[ii]-0.025,fct_mean[ii],'s',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.8')
   
    ii = where(iio30 and iir10 and iig18)
    p = errorplot(f25[ii]-0.015,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p4 = plot(f25[ii]-0.015,fct_mean[ii],'s',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[215,48,39],/ov,name='$R$ = 1.0, $\Gamma$ = 1.8')
   
    ii = where(iio30 and iir05 and iig19)
    p = errorplot(f25[ii]-0.005,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p2 = plot(f25[ii]-0.005,fct_mean[ii],'o',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.9')
    
    ii = where(iio30 and iir10 and iig19)
    p = errorplot(f25[ii]+0.005,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p5 = plot(f25[ii]+0.005,fct_mean[ii],'o',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[215,48,39],/ov,name='$R$ = 1.0, $\Gamma$ = 1.9')
    
    ii = where(iio30 and iir05 and iig20)
    p = errorplot(f25[ii]+0.015,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p3 = plot(f25[ii]+0.015,fct_mean[ii],'S',_extra=efct,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 2.0')
    
    ii = where(iio30 and iir10 and iig20)
    p = errorplot(f25[ii]+0.025,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p6 = plot(f25[ii]+0.025,fct_mean[ii],'S',_extra=efct,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color=[215,48,39],/ov,name='$R$ = 1.0, $\Gamma$ = 2.0')

    title = text(0.18,0.80,/normal,'Opening Angle 30$\deg$',font_name='Times',font_size=15)
    l = legend(target=[p1,p2,p3,p4,p5,p6],position=[0.88,0.85],/normal,horizontal_spacing=0.06,sample_width=0.0,font_size=14,font_name='Times')
     ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then poa30.save,'figures/param_fct_oa30.eps',/BITMAP else $
                                                     poa30.save,'figures/param_fct_oa30.png',resolution=res
    endif    

    ;; OPENING ANGLE 60°
    ii = where(iio60 and iir05 and iig18)
    poa60 = errorplot(f25[ii]-0.025,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',errorbar_thick=2)
    p1 = plot(f25[ii]-0.025,fct_mean[ii],'s',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.8')
   
    ii = where(iio60 and iir10 and iig18)
    p = errorplot(f25[ii]-0.015,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p4 = plot(f25[ii]-0.015,fct_mean[ii],'s',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[215,48,39],/ov,name='$R$ = 1.0, $\Gamma$ = 1.8')
   
    ii = where(iio60 and iir05 and iig19)
    p = errorplot(f25[ii]-0.005,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p2 = plot(f25[ii]-0.005,fct_mean[ii],'o',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.9')
    
    ii = where(iio60 and iir10 and iig19)
    p = errorplot(f25[ii]+0.005,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p5 = plot(f25[ii]+0.005,fct_mean[ii],'o',_extra=efct,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[215,48,39],/ov,name='$R$ = 1.0, $\Gamma$ = 1.9')
    
    ii = where(iio60 and iir05 and iig20)
    p = errorplot(f25[ii]+0.015,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p3 = plot(f25[ii]+0.015,fct_mean[ii],'S',_extra=efct,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 2.0')
    
    ii = where(iio60 and iir10 and iig20)
    p = errorplot(f25[ii]+0.025,fct_mean[ii],fct_err[*,ii],_extra=efct,linestyle='',/ov,errorbar_thick=2)
    p6 = plot(f25[ii]+0.025,fct_mean[ii],'S',_extra=efct,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color=[215,48,39],/ov,name='$R$ = 1.0, $\Gamma$ = 2.0')

    title = text(0.18,0.80,/normal,'Opening Angle 60$\deg$',font_name='Times',font_size=15)
    l = legend(target=[p1,p2,p3,p4,p5,p6],position=[0.88,0.85],/normal,horizontal_spacing=0.06,sample_width=0.0,font_size=14,font_name='Times')
     ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then poa60.save,'figures/param_fct_oa60.eps',/BITMAP else $
                                                     poa60.save,'figures/param_fct_oa60.png',resolution=res
    endif    


    ;; FSCAT vs F24 BY PARAM
    efscat = {xra:[0.1,0.9],yra:[-1.5,-0.5], $
              sym_color:'black', $;sym_filled:1,sym_size:1.2, $
              ;xtitle:'$!8f!7_{24-25}$',ytitle:'$\sigma_{scatt}$', $
              xtitle:'$!8f!7_{25-26}$',ytitle:'$\sigma_{scatt}$', $
              font_name:'Times',font_size:14}

    ;; OPENING ANGLE 30°
    ii = where(iio30 and iir05 and iig18)
    poa30 = errorplot(f25[ii]-0.025,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',errorbar_thick=2)
    p1 = plot(f25[ii]-0.025,fscat_mean[ii],'s',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.8')
   
    ii = where(iio30 and iir10 and iig18)
    p = errorplot(f25[ii]-0.015,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p4 = plot(f25[ii]-0.015,fscat_mean[ii],'s',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[69,117,180],/ov,name='$R$ = 1.0, $\Gamma$ = 1.8')
   
    ii = where(iio30 and iir05 and iig19)
    p = errorplot(f25[ii]-0.005,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p2 = plot(f25[ii]-0.005,fscat_mean[ii],'o',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.9')
    
    ii = where(iio30 and iir10 and iig19)
    p = errorplot(f25[ii]+0.005,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p5 = plot(f25[ii]+0.005,fscat_mean[ii],'o',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[69,117,180],/ov,name='$R$ = 1.0, $\Gamma$ = 1.9')
    
    ii = where(iio30 and iir05 and iig20)
    p = errorplot(f25[ii]+0.015,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p3 = plot(f25[ii]+0.015,fscat_mean[ii],'S',_extra=efscat,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 2.0')
    
    ii = where(iio30 and iir10 and iig20)
    p = errorplot(f25[ii]+0.025,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p6 = plot(f25[ii]+0.025,fscat_mean[ii],'S',_extra=efscat,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color=[69,117,180],/ov,name='$R$ = 1.0, $\Gamma$ = 2.0')

    title = text(0.62,0.80,/normal,'Opening Angle 30$\deg$',font_name='Times',font_size=15)
    l = legend(target=[p1,p2,p3,p4,p5,p6],position=[0.42,0.85],/normal,horizontal_spacing=0.06,sample_width=0.0,font_size=14,font_name='Times')
     ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then poa30.save,'figures/param_fscat_oa30.eps',/BITMAP else $
                                                     poa30.save,'figures/param_fscat_oa30.png',resolution=res
    endif    

    ;; OPENING ANGLE 60°
    ii = where(iio60 and iir05 and iig18)
    poa60 = errorplot(f25[ii]-0.025,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',errorbar_thick=2)
    p1 = plot(f25[ii]-0.025,fscat_mean[ii],'s',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.8')
   
    ii = where(iio60 and iir10 and iig18)
    p = errorplot(f25[ii]-0.015,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p4 = plot(f25[ii]-0.015,fscat_mean[ii],'s',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[69,117,180],/ov,name='$R$ = 1.0, $\Gamma$ = 1.8')
   
    ii = where(iio60 and iir05 and iig19)
    p = errorplot(f25[ii]-0.005,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p2 = plot(f25[ii]-0.005,fscat_mean[ii],'o',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 1.9')
    
    ii = where(iio60 and iir10 and iig19)
    p = errorplot(f25[ii]+0.005,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p5 = plot(f25[ii]+0.005,fscat_mean[ii],'o',_extra=efscat,sym_size=1.5,sym_filled=1,sym_thick=2,col='black',sym_fill_color=[69,117,180],/ov,name='$R$ = 1.0, $\Gamma$ = 1.9')
    
    ii = where(iio60 and iir05 and iig20)
    p = errorplot(f25[ii]+0.015,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p3 = plot(f25[ii]+0.015,fscat_mean[ii],'S',_extra=efscat,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color='white',/ov,name='$R$ = 0.5, $\Gamma$ = 2.0')
    
    ii = where(iio60 and iir10 and iig20)
    p = errorplot(f25[ii]+0.025,fscat_mean[ii],fscat_err[*,ii],_extra=efscat,linestyle='',/ov,errorbar_thick=2)
    p6 = plot(f25[ii]+0.025,fscat_mean[ii],'S',_extra=efscat,sym_size=2.0,sym_thick=2,col='black',sym_filled=1,sym_fill_color=[69,117,180],/ov,name='$R$ = 1.0, $\Gamma$ = 2.0')

    title = text(0.62,0.80,/normal,'Opening Angle 60$\deg$',font_name='Times',font_size=15)
    l = legend(target=[p1,p2,p3,p4,p5,p6],position=[0.42,0.85],/normal,horizontal_spacing=0.06,sample_width=0.0,font_size=14,font_name='Times')
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then poa60.save,'figures/param_fscat_oa60.eps',/BITMAP else $
                                                     poa60.save,'figures/param_fscat_oa60.png',resolution=res
    endif     
endif




;;----------------------------------------------------------------------------------------
;; MRLX plots
;;----------------------------------------------------------------------------------------
if keyword_set(mrlx) then begin
    
    
    
    ;; MRLX solo
    readcol,'mrlx_dist.csv',grid_of_rlx,mrlx_vals,format='d,d'
    
    e = {xra:[-3.2,0.5],yra:[0.,0.012], $
         thick:3, $
         xtitle:'$!8M!7(!8R_{L_{!7X}})$',ytitle:'Probability', $
         font_name:'Times',font_size:14}
         
    p = plot(grid_of_rlx,mrlx_vals/total(mrlx_vals),_extra=e,col=[69,117,180])
    t = text(0.18,0.80,'$!8M!7(!8R_{L_{!7X}})$',/normal,font_name='Times',font_style='Bold',font_size=14,fill_background=1)
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/mrlx_dist.eps',/BITMAP else $
                                                     p.save,'figures/mrlx_dist.png',resolution=res
    endif
    stop
    ;; convolution with MRLX
    
    ;; RX and RX-LIMIT data for sample
    rx_data = read_csv('rx_data.csv',header=hd)
    for i = 0,n_elements(hd)-1 do re = execute(hd[i]+' = rx_data.field'+string(i+1,'(i1)'))
    ;; detected sources
    iid = rx ne -99.

    ;; evenly spaced grid of rlx values
    grid_num = 600.
    grid_of_rlx = [-3.:1:4./grid_num]

    readcol,'mrlx_vals.csv',mrlx_vals,format='d'
    readcol,'first_term_rlx_grid.csv',f1,f2,f3,f4,f5,f6,format='d'
    first_term_rlx = [[f1],[f2],[f3],[f4],[f5],[f6]]
    readcol,'second_term_rlx_grid.csv',s1,s2,s3,s4,s5,s6,format='d'
    second_term_rlx = [[s1],[s2],[s3],[s4],[s5],[s6]]

    col = [ $
           ;[165,0,38],$
           [215,48,39],$
           [244,109,67],$
           [253,174,97],$
           ;[254,224,144],$
           ;[255,255,191],$
           ;[224,243,248],$
           [171,217,233],$
           [116,173,209],$
           [69,117,180] $
           ;[49,54,149] $
           ]
        
    e = {xra:[-3.,1.0],yra:[-1.,13.], $
         thick:3, $
         ytickname:['','','','','','','',''], $
         xtitle:'$!8R_{L_{!7X}}$',ytitle:'Likelihood (arbitrary)', $
         font_name:'Times',font_size:14}
    
    ;; FIRST TERM RLX
    p = plot(grid_of_rlx,first_term_rlx[*,0],_extra=e,col=col[*,0])  
    for i = 1,5 do p = plot(grid_of_rlx,first_term_rlx[*,i]+i*2,_extra=e,/ov,col=col[*,i],transparency=0)    
    t = text(0.18,0.80,'$!8N!7(!8R!7_{!8L!7_{X}, !8i!7} , 0.23)$',/normal,font_name='Times',font_style='Bold',font_size=14,fill_background=1)
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/first_term_rlx.eps',/BITMAP else $
                                                     p.save,'figures/first_term_rlx.png',resolution=res
    endif

    ;; MRLX CONVOLVED WITH FIRST TERM RLX
    p = plot(grid_of_rlx,mrlx_vals*first_term_rlx[*,0],_extra=e,col=col[*,0])  
    for i = 1,5 do p = plot(grid_of_rlx,mrlx_vals*first_term_rlx[*,i]+i*2,_extra=e,/ov,col=col[*,i],transparency=0)    
    t = text(0.18,0.80,'$!8M!7(!8R!7_{!8L!7_{X}}) !M* !8N!7(!8R!7_{!8L!7_{X}, !8i!7} , 0.23)$',/normal,font_name='Times',font_style='Bold',font_size=14,fill_background=1)
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/first_term_convolved.eps',/BITMAP else $
                                                     p.save,'figures/first_term_convolved.png',resolution=res
    endif

    ;; SECOND TERM RLX
    p = plot(grid_of_rlx,second_term_rlx[*,0],_extra=e,col=col[*,0])  
    for i = 1,5 do p = plot(grid_of_rlx,second_term_rlx[*,i]+i*2,_extra=e,/ov,col=col[*,i],transparency=0)    
    t = text(0.18,0.80,'$!8S_i!7$',/normal,font_name='Times',font_style='Bold',font_size=14,fill_background=1)
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/second_term_rlx.eps',/BITMAP else $
                                                     p.save,'figures/second_term_rlx.png',resolution=res
    endif

    ;; MRLX CONVOLVED WITH SECOND TERM RLX
    p = plot(grid_of_rlx,mrlx_vals*second_term_rlx[*,0],_extra=e,col=col[*,0])  
    for i = 1,5 do p = plot(grid_of_rlx,mrlx_vals*second_term_rlx[*,i]+i*2,_extra=e,/ov,col=col[*,i],transparency=0)    
    t = text(0.18,0.80,'$!8M!7(!8R!7_{!8L!7_{X}}) !M* !8S_i!7$',/normal,font_name='Times',font_style='Bold',font_size=14,fill_background=1)
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/second_term_convolved.eps',/BITMAP else $
                                                     p.save,'figures/second_term_convolved.png',resolution=res
    endif
endif



;;----------------------------------------------------------------------------------------
;; CORNER PLOT
;;----------------------------------------------------------------------------------------
if keyword_set(corner) then begin


    chain = read_csv('sampler_chain.csv')
    fscat = chain.field1
    fct = chain.field2

    e = {font_name:'Times',font_size:18}

    ;; load contour 
    im = image('figures/output.png',dimension=[1028,1036],position=[50,50,978,986],/device)

    ;; overplot contours
    p = plot(fscat,fct,/nodata,/current,position=[0.116,0.097,0.721,0.6975],_extra=e)
    ;; overplot fscat hist
    hfs = histogram(fscat,locations=xfs,binsize=0.2,min=-1.8,max=-0.6)
    pfs = plot(xfs,hfs/total(hfs),/stairstep,/nodata,/current,position=[0.116,0.718,0.721,0.918],_extra=e, $
               xmajor=7,ymajor=4,xshowtext=0,yshowtext=0)
    ;; overplot fct hist
    hfct = histogram(fct,locations=xfct,binsize=0.05,min=0.40,max=0.70)

    pfct = plot(xfct,hfct/total(hfct),/stairstep,/nodata,/current,position=[0.7415,0.097,0.9430,0.6975],_extra=e, $
                xmajor=4,ymajor=7,xshowtext=0,yshowtext=0)
    ;; remove ticks and text
    ;;axes
    t = text(0.400,0.050,'$\sigma_{scatt}$',fill_background=1,fill_color='white',font_size=22)
    t = text(0.042,0.385,'$!8f!7_{CT}$',fill_background=1,fill_color='white',font_size=22,color='black',orientation=90)  
    ;; labels
    t = text(0.345,0.936,'$\sigma_{scatt} = -1.04^{+0.13}_{-0.12}$',font_size=20,font_name='Times',fill_background=1,fill_color='white')
    t = text(0.772,0.716,'$!8f!7_{CT} = 0.555^{+0.037}_{-0.032}$',font_size=20,font_name='Times',fill_background=1,fill_color='white')
    
    ;; save image
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/corner_plot.eps',/BITMAP else $
                                                     p.save,'figures/corner_plot.png',resolution=res
    endif



endif





END










