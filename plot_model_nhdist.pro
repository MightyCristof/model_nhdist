PRO plot_model_nhdist, PROPERTIES = properties, $
                       NHDIST = nhdist, $
                       RXDIST = rxdist, $
                       XSPEC = xspec, $
                       HIDE = hide, $
                       LOW_RES = low_res, $
                       SAV = sav


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
common _fixed
common _free
common _split
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
    
    bn = 0.05
    yh = histogram(z[ix],bin=bn,location=xh,min=0.,max=0.8)
    yhd = histogram(z[id],bin=bn,location=xhd,min=0.,max=0.8)
    yhn = histogram(z[in],bin=bn,location=xhn,min=0.,max=0.8)

    e = {xra:[0.,0.85],yra:[0.,120.],$
         stairstep:1,thick:2,fill_background:1, $
         xtitle:'$!8z!7$',ytitle:'Frequency', $
         font_name:'Times',font_size:16, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
    
    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources',fill_transparency=30)
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected',fill_transparency=15)
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.',fill_transparency=20)
    leg = legend(target=[p,pd,pn],position=[0.61,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=12,font_name='Times')
    t = text(0.03,0.90,"a",/normal,font_name='Times',font_style='Bold',font_size=16)
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/z_dist.eps',/BITMAP else $
                                                     p.save,'figures/z_dist.png',resolution=res
    endif    

    bn = 0.25
    yh = histogram(loglir[ix],bin=bn,location=xh,min=41.,max=46.)
    yhd = histogram(loglir[id],bin=bn,location=xhd,min=41.,max=46.)
    yhn = histogram(loglir[in],bin=bn,location=xhn,min=41.,max=46.)
    
    e = {xra:[42.,46.],yra:[0.,200.],$
         stairstep:1,thick:2,fill_background:1, $
         xtitle:'$log !8L!7_{MIR}  [erg s^{-1} cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times', $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources',fill_transparency=30)
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected',fill_transparency=15)
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.',fill_transparency=20)
    leg = legend(target=[p,pd,pn],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=12,font_name='Times')
    t = text(0.03,0.90,"b",/normal,font_name='Times',font_style='Bold',font_size=16)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/lir_dist.eps',/BITMAP else $
                                                     p.save,'figures/lir_dist.png',resolution=res
    endif    

    bn = 0.25
    yh = histogram(loglxir[ix],bin=bn,location=xh,min=41.,max=46.)
    yhd = histogram(loglxir[id],bin=bn,location=xhd,min=41.,max=46.)
    yhn = histogram(loglxir[in],bin=bn,location=xhn,min=41.,max=46.)
    
    e = {xra:[42.,46.],yra:[0.,200.],$
         stairstep:1,thick:2,fill_background:1, $
         xtitle:'$log !8L!7_{X}(!8L!7_{MIR})  [erg s^{-1} cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times', $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources',fill_transparency=30)
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected',fill_transparency=15)
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.',fill_transparency=20)
    leg = legend(target=[p,pd,pn],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=12,font_name='Times')
    t = text(0.03,0.90,"c",/normal,font_name='Times',font_style='Bold',font_size=16)

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
    modoff = 0.85*[1.,1.]
    stack_flux = fxstak[1,0:1]
    stack_err = e_fxstak[1,0:1]
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    ;; variables for FIXED vs FREE
    fx_models = ['FX_NON','FX_NON_']
    nh_models = ['NH_MOD','NH_MOD_']
    ctf_values = ['FCT','FCT_']
    ctf_errors = ['E_FCT','E_FCT_']
    
    ;; plot set up
    dim = [880,780]
    sq = 180
    gap = 40
    pos = [[100,560,100+sq,560+sq],[100+sq+gap,560,100+2.*sq+gap,560+sq],[600,560,840,560+sq],$
                                   [100,70,840,480]]
    e = {font_name:'Times', $
         dim:dim,device:1,buffer:0}
    enh = {xra:[20.,26.],yra:[0.,0.8], $
           stairstep:1, $
           xtitle:'$log !8N!7_{H}  [cm^{-2}]$',ytitle:'Frequency (normalized)', $
           font_name:'Times',font_size:14, $
           current:1,dim:dim,device:1}
    if keyword_set(hide) then e.buffer = 1

    ;; color
    col = [106,81,163];[84,39,143]
            ;; purple     [84,39,143]
            ;; teal       [28,144,153]
            ;; mint       [127,205,187]
            ;; aquamarine [28,144,153]

    for i = 0,1 do begin
        ;; set FIXED vs FREE data for plotting
        re = execute('fxmod = '+fx_models[i])
        fxm = [mode(fxmod.soft,kde=kde_bandwidth(fxmod.soft)),mode(fxmod.hard,kde=kde_bandwidth(fxmod.hard))]
        ;e_fxm = [mode(fxmod.e_soft,kde=kde_bandwidth(fxmod.e_soft)),mode(fxmod.e_hard,kde=kde_bandwidth(fxmod.e_hard))]
        e_fxm = [stddev(fxmod.soft),stddev(fxmod.hard)]
        re = execute('nhm = '+nh_models[i])
        ilo = where(nhm.xh lt 20.)
        nhm[ilo[-1]+1].yh += total(nhm[ilo].yh)
        nhm[ilo].yh = 0.
        nhm[ilo[-1]+1].sig += sqrt(total(nhm[ilo[0]:ilo[-1]+1].sig^2.))
        nhm[ilo].sig = 0.
        re = execute('ctf = string('+ctf_values[i]+',format="(d0.2)")+"\pm"+string('+ctf_errors[i]+',format="(d0.2)")')
        ;; fix NH YRA
        enh.yra[1] = 0.8>ceil(max(nhm.yh))
    
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
        pxm = errorplot(energy_center*modoff,fxm,e_fxm, $
                        linestyle='',errorbar_capsize=0.1,/ov)
        pxm = plot(energy_center*modoff,fxm,'o',col='black', $
                   sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
        pxm = plot(energy_center*modoff,fxm,'o',col='black', $
                   sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=col,sym_transparency=20,fill_transparency=20,/ov)
        to = text(0.06,0.85,target=pxm,/relative,'X-ray non-det. sources',font_name='Times')
        po = plot([1.5e-1],[3.65e-15],'o',col='black',sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=col,sym_transparency=20,/ov)
        to = text(0.15,0.75,target=pxm,/relative,'Model',font_name='Times')
        po = plot([1.5e-1],[2.32e-15],'s',col='black',sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
        to = text(0.15,0.65,target=pxd,/relative,'X-ray stack',font_name='Times')

        ;; plot NH distribution
        pct = plot([24.,26.],[1.,1.]*enh.yra[1],_extra=enh,pos=pos[*,3], $
                   fill_background=1,fill_color='light grey',fill_level=0.,fill_transparency=0)
        pct.stairstep = 0   ;; CT shading
        pnh = errorplot(nhm.xhoff,nhm.yh,nhm.sig,'-',thick=4,errorbar_thick=4,_extra=enh, $
                        pos=pos[*,3],fill_background=1,fill_color=col,fill_transparency=30,name='This work')
        pan = plot(nh_ana_lox.xh,nh_ana_lox.yh*frac_lox+nh_ana_hix.yh*frac_hix,'--',thick=2,_extra=enh,/ov,name='$Ananna+2019 (weighted)$')
        tct = text(24.1,0.06,'CT',target=pnh,/data,font_name='Times',font_size=14)
        ctad = text(25.,(pnh.yra[1]+max(nhm.yh+nhm.sig))/2.,'$!8f!7_{CT} = '+ctf+'$',target=pnh,/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=0.5)
        leg = legend(target=[pnh,pan],position=[0.125,0.60],/normal,horizontal_alignment=0.,font_size=14,font_name='Times')

        ;; figure reference for Nature caption
        t = text(0.03,0.97,"a",/normal,font_name='Times',font_style='Bold',font_size=16)
        t = text(0.03,0.62,"b",/normal,font_name='Times',font_style='Bold',font_size=16)
        
        ;; save image
        if keyword_set(sav) then begin
            print, '    SAVING PLOT'
            if (i eq 0) then savfile = 'nh_fixed' else $
                             savfile = 'nh_free'
            if (strupcase(strtrim(sav,2)) eq 'EPS') then pnh.save,'figures/'+savfile+'.eps',/BITMAP else $
                                                         pnh.save,'figures/'+savfile+'.png',resolution=res
        endif    
    endfor
endif


;;----------------------------------------------------------------------------------------
;; DISTRIBUTION OF BEST-FIT RX FROM MODELING
;;----------------------------------------------------------------------------------------
if keyword_set(rxdist) then begin

    data = ['RXDV','RXDV1']
    models = ['RXMV','RXMV1']
    indices = ['IIMV','IIMV1']
    ;models = ['RX_MODV','RX_MODV_']
    ;indices = ['IIMODV','IIMODV_']
    letter = ['a','b']
    file = ['rx_fixed','rx_free']
    
    rxbin = 0.2d
    
    for i = 0,1 do begin
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
             font_name:'Times',font_size:16, $
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
        t = text(0.03,0.90,letter[i],/normal,font_name='Times',font_style='Bold',font_size=16)

        if keyword_set(sav) then begin
            print, '    SAVING PLOT'
            if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/'+file[i]+'.eps',/BITMAP else $
                                                         p.save,'figures/'+file[i]+'.png',resolution=res
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
    line = ['-','--','-.',':','__']
    ;col = [[213,94,0],[204,121,167],[0,114,178],[240,228,66],[0,158,115]]

    ;;      TOTAL       PL         SCPL        BORUS
    ;;      GREEN       ORANGE     BLUE        MAGENTA
    col = [[0,158,115],[213,94,0],[0,114,178],[204,121,167]]

    e = {xra:[0.1,200],yra:[1e-4,50], $
         xlog:1,ylog:1, $
         xtickname:['0.1','1','10','100'],ytickname:'10!U'+strtrim([-4:1],2), $
         ytickinterval:1, $
         xtitle:'Energy  [keV]',ytitle:'!8EF$_E$!7  [keV Photons cm$^{-2}$ s$^{-1}$]', $
         font_name:'Times',font_size:16}
    p = plot(energy,tot21,_extra=e,/nodata)
    thick = 2
    ;; total SED vs. NH
    re = execute('pt21 = plot(energy,'+tot[0]+',linestyle=line[0],col=col[*,0],thick=thick,/ov,name="Total log !8N!7$_H$ = 21")')
    re = execute('pt23 = plot(energy,'+tot[1]+',linestyle=line[1],col=col[*,0],thick=thick,/ov,name="Total log !8N!7$_H$ = 23")')
    re = execute('pt25 = plot(energy,'+tot[2]+',linestyle=line[2],col=col[*,0],thick=thick,/ov,name="Total log !8N!7$_H$ = 25")')
    ;for i = 0,nnh-1 do re = execute('pt = plot(energy,'+tot[i]+',linestyle=line[i],col=col[*,0],/ov)')
    re = execute('pp = plot(energy,'+pl[0]+',linestyle=line[0],col=col[*,1],thick=thick,/ov,name="Direct")')
    for i = 1,nnh-1 do re = execute('p = plot(energy,'+pl[i]+',linestyle=line[i],col=col[*,1],thick=thick,/ov)')
    re = execute('ps = plot(energy,'+scpl[0]+',linestyle=line[0],col=col[*,2],thick=thick,/ov,name="Scattered")')
    for i = 1,nnh-1 do re = execute('p = plot(energy,'+scpl[i]+',linestyle=line[i],col=col[*,2],thick=thick,/ov)')
    re = execute('pb = plot(energy,'+bor[0]+',linestyle=line[0],col=col[*,3],thick=thick,/ov,name="Torus")')
    for i = 1,nnh-1 do re = execute('p = plot(energy,'+bor[i]+',linestyle=line[i],col=col[*,3],thick=thick,/ov)')
    ;; legend
    leg = legend(target=[pt21,pt23,pt25],position=[0.43,0.87],horizontal_spacing=0.06,font_name='Times',font_size=14)
    leg = legend(target=[pp,ps,pb],position=[0.86,0.30],horizontal_spacing=0.06,font_name='Times',font_size=14)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/xspec_model.eps',/BITMAP else $
                                                     p.save,'figures/xspec_model.png',resolution=res
    endif        
endif







END










