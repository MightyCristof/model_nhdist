PRO plot_model_nhdist, PROPERTIES = properties, $
                       NHDIST = nhdist, $
                       RXDIST = rxdist, $
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

    iix = iiwac
    iid = iiwd
    iin = iiwn
    ;iix = bytarr(n_elements(z))+1b
    ;iid = iiad
    ;iin = iian
    
    bn = 0.05
    yh = histogram(z[where(iix)],bin=bn,location=xh,min=0.,max=0.8)
    yhd = histogram(z[where(iid)],bin=bn,location=xhd,min=0.,max=0.8)
    yhn = histogram(z[where(iin)],bin=bn,location=xhn,min=0.,max=0.8)

    e = {xra:[0.,0.85],yra:[0.,110.],$
         stairstep:1,thick:2,fill_background:1,fill_transparency:25, $
         xtitle:'$!8z!7$',ytitle:'Frequency', $
         font_name:'Times',font_size:16, $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1
    
    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources')
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.')
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected')
    leg = legend(target=[pd,pn,p],position=[0.61,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=12,font_name='Times')
    t = text(0.03,0.90,"a",/normal,font_name='Times',font_style='Bold',font_size=16)
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/z_dist.eps',/BITMAP else $
                                                     p.save,'figures/z_dist.png',resolution=res
    endif    

    bn = 0.25
    yh = histogram(loglir[where(iix)],bin=bn,location=xh,min=41.,max=46.)
    yhd = histogram(loglir[where(iid)],bin=bn,location=xhd,min=41.,max=46.)
    yhn = histogram(loglir[where(iin)],bin=bn,location=xhn,min=41.,max=46.)
    
    e = {xra:[42.,46.],yra:[0.,150.],$
         stairstep:1,thick:2,fill_background:1,fill_transparency:25, $
         xtitle:'$log !8L!7_{MIR} [erg s^{-1} cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times', $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources')
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.')
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected')
    leg = legend(target=[pd,pn,p],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=12,font_name='Times')
    t = text(0.03,0.90,"b",/normal,font_name='Times',font_style='Bold',font_size=16)

    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/lir_dist.eps',/BITMAP else $
                                                     p.save,'figures/lir_dist.png',resolution=res
    endif    

    bn = 0.25
    yh = histogram(loglxir[where(iix)],bin=bn,location=xh,min=41.,max=46.)
    yhd = histogram(loglxir[where(iid)],bin=bn,location=xhd,min=41.,max=46.)
    yhn = histogram(loglxir[where(iin)],bin=bn,location=xhn,min=41.,max=46.)
    
    e = {xra:[42.,46.],yra:[0.,150.],$
         stairstep:1,thick:2,fill_background:1,fill_transparency:25, $
         xtitle:'$log !8L!7_{X}(!8L!7_{MIR}) [erg s^{-1} cm^{-2}]$',ytitle:'Frequency', $
         font_name:'Times', $
         buffer:0}
    if keyword_set(hide) then e.buffer = 1

    p = plot(xh+bn/2.,yh,':',_extra=e,fill_color='light grey',name='All sources')
    pn = plot(xhn+bn/2.,yhn,'-',_extra=e,fill_color='orange',/ov,name='X-ray non-det.')
    pd = plot(xhd+bn/2.,yhd,'--',_extra=e,fill_color='dodger blue',/ov,name='X-ray detected')
    leg = legend(target=[pd,pn,p],position=[0.15,0.86],/normal,sample_width=0.16,horizontal_spacing=0.06,horizontal_alignment=0.,font_size=12,font_name='Times')
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

    ;; set up NH model for plotting
    var = ['XNH','XNHO','YNH','ENH','MNH','FCT']
    ;; FIXED
    xnh = nh_mod.xh
    xnho = nh_mod.xhoff
    ynh = nh_mod.yh
    enh = nh_mod.sig
    ;mnh = (hist2d_avg(nhmv,1.0d,normalize='where(xh lt 24.,/null)')).mad
    i20 = where(xnh lt 20.,n20)
    i24 = where(xnh ge 24.,n24)
    if (n20 gt 0) then begin
        ynh[i20[-1]+1] = total([ynh[i20],ynh[i20[-1]+1]])
        ynh[i20] = 0.
        enh[i20[-1]+1] = sqrt(total(([enh[i20],enh[i20[-1]+1]])^2.))
        enh[i20] = 0.
    endif
    ;; FREE
    xnh_ = nh_mod_.xh
    xnho_ = nh_mod_.xhoff
    ynh_ = nh_mod_.yh
    enh_ = nh_mod_.sig
    ;mnh_ = (hist2d_avg(nhmv2,1.0d,normalize='where(xh lt 24.,/null)')).mad
    i20 = where(xnh_ lt 20.,n20)
    i24 = where(xnh_ ge 24.,n24)
    if (n20 gt 0) then begin
        ynh_[i20[-1]+1] = total([ynh_[i20],ynh_[i20[-1]+1]])
        ynh_[i20] = 0.
        enh_[i20[-1]+1] = sqrt(total(([enh_[i20],enh_[i20[-1]+1]])^2.))
        enh_[i20] = 0.
    endif

    ;; X-ray fluxes    
    energy_center = [1.25,4.5]
    energy_range = [0.75,2.5]
    modoff = 0.85*[1.,1.]
    stack_flux = fxstak[1,0:1]
    stack_err = e_fxstak[1,0:1]
    xy = findgen(35,start=-17)
    ra = minmax(xy)

    ;; flux estimates
    fx_mod = [mode(fx_non.soft,bin=scott(fx_non.soft)),mode(fx_non.hard,kde=kde_bandwidth(fx_non.hard))]
    e_fx_mod = sqrt([mode(fx_non.e_soft,kde=kde_bandwidth(fx_non.e_soft)),mode(fx_non.e_hard,kde=kde_bandwidth(fx_non.e_hard))]^2.+([stddev(fx_non.soft),stddev(fx_non.hard)])^2.)

    dim = [880,780]
    sq = 180
    gap = 40
    pos = [[100,560,100+sq,560+sq],[100+sq+gap,560,100+2.*sq+gap,560+sq],[600,560,840,560+sq],$
                                   [100,70,840,480]]
    e = {xra:[20.,26.],yra:[0.,rnd(max(ynh)>max(ynh_),1)+0.1>1.1],$
         stairstep:1, $
         xtitle:'$log !8N!7_{H} [cm^{-2}]$',ytitle:'Frequency (normalized)', $
         font_name:'Times',font_size:14, $
         dim:dim,device:1,buffer:0}
    if keyword_set(hide) then e.buffer = 1
    
    nhcol = [28,144,153]
            ;; purple     [84,39,143]
            ;; teal       [28,144,153]
            ;; mint       [127,205,187]
            ;; aquamarine [28,144,153]
    
    for i = 0,1 do begin
        if (i eq 0) then begin
            xxnh = xnho
            yynh = ynh
            eenh = enh
            ;mmnh = mnh
            cctf = string(fct,format='(d4.2)')+'\pm'+string(e_fct,format='(d4.2)')
        endif else begin
            xxnh = xnho_
            yynh = ynh_
            eenh = enh_
            ;mmnh = mnh_
            cctf = string(fct_,format='(d4.2)')+'\pm'+string(e_fct_,format='(d4.2)')
        endelse
        ;pad = errorplot(xxnh[i24],yynh[i24],mmnh[i24],linestyle='',_extra=e, $
        ;                pos=pos[*,3],errorbar_color='black',errorbar_linestyle='-.')
        pad = errorplot(xxnh,yynh,eenh,'-',thick=4,errorbar_thick=4,_extra=e, $
                        pos=pos[*,3],fill_background=1,fill_color=nhcol,fill_transparency=40,name='This work')
        pavg = plot(nh_ana_lox.xh,nh_ana_lox.yh*frac_lox+nh_ana_hix.yh*frac_hix,'--',thick=2,_extra=e,/ov,name='$Ananna+2019 (weighted)$')
        ;; CT line
        pct = plot([24.,24.],e.yra,_extra=e,thick=2,linestyle=':',/ov)
        ;ctl = text(alog10(1.5e24)+0.1,1.0,'CT',/data,font_size=16,font_name='Times')
        ;; CT fraction
        ctad = text(25.,(e.yra[1]+max(yynh+eenh))/2.,'$!8f!7_{CT} = '+cctf+'$',/data,font_size=16,font_name='Times',alignment=0.5,vertical_alignment=1.)
        ;; add legend
        leg = legend(target=[pavg,pad],position=[0.125,0.60],/normal,horizontal_alignment=0.,font_size=14,font_name='Times')
        ;; add X-ray stacked images
        psoft = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,0],/current,/device,ytitle='offset in Decl. [arcsec.]',font_name='Times',font_size=14)
        imsoft = image(soft,pos=pos[*,0],/current,/device)    
        lsoft = text(target=psoft,0.5,0.82,'  0.5$-$2 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
        phard = plot(xy,xy,/nodata,xra=ra,yra=ra,pos=pos[*,1],/current,/device);,xtitle='offset in R.A. [arcsec.]                     ')
        imhard = image(hard,pos=pos[*,1],/current,/device)
        lhard = text(target=phard,0.5,0.82,'  2$-$7 keV  ',font_name='Times',font_style='Bold',font_size=16,fill_background=1,alignment=0.5,/relative)
        xt = text(target=plo,(pos[2,0]+pos[0,1])/2.,pos[1,0]-42,/device,'offset in R.A. [arcsec.]',alignment=0.5,font_name='Times',font_size=14)
        
        ;; add X-ray flux estimates
        p = errorplot(energy_center,stack_flux,energy_range,stack_err, $
                       linestyle='',/xlog,/ylog,pos=pos[*,2],/current,/device, $
                       xtitle='$energy [keV]$',ytitle='$!8F!7_{X}  [erg s^{-1} cm^{-2}]$', $
                       font_name='Times',font_size=14)
                       ;xtickunit='exponent'
        px = plot(energy_center,stack_flux,'s',col='black', $
                  sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',/ov,name='Stack')
        p.xra=[1e-1,1e1]
        p.yra=[1e-16,1e-14]
        p = errorplot(energy_center*modoff,fx_mod,e_fx_mod, $
                      linestyle='',errorbar_capsize=0.1,/ov)
        pad = plot(energy_center*modoff,fx_mod,'o',col='black', $
                   sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
        pad = plot(energy_center*modoff,fx_mod,'o',col='black', $
                   sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=nhcol,sym_transparency=20,/ov)
        ;; add "legend"  
        to = text(0.06,0.85,target=pad,/relative,'X-ray non-det. sources',font_name='Times')
        po = plot([1.5e-1],[3.65e-15],'o',col='black',sym_size=1.4,sym_thick=1.5,sym_filled=1,sym_fill_color=nhcol,sym_transparency=20,/ov)
        to = text(0.15,0.75,target=pad,/relative,'Model',font_name='Times')
        po = plot([1.5e-1],[2.32e-15],'s',col='black',sym_size=1.2,sym_thick=1.5,sym_filled=1,sym_fill_color='white',sym_transparency=0,/ov)
        to = text(0.15,0.65,target=px,/relative,'X-ray stack',font_name='Times')
        ;; figure reference for Nature caption
        t = text(0.03,0.97,"a",/normal,font_name='Times',font_style='Bold',font_size=16)
        t = text(0.03,0.62,"b",/normal,font_name='Times',font_style='Bold',font_size=16)
        
        if keyword_set(sav) then begin
            print, '    SAVING PLOT'
            if (i eq 0) then savfile = 'nh_fixed' else $
                             savfile = 'nh_free'
            if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/'+savfile+'.eps',/BITMAP else $
                                                         p.save,'figures/'+savfile+'.png',resolution=res
        endif    
    endfor
endif


;;----------------------------------------------------------------------------------------
;; DISTRIBUTION OF BEST-FIT RX FROM MODELING
;;----------------------------------------------------------------------------------------
if keyword_set(rxdist) then begin

    models = ['RXMV','RXMV1']
    indices = ['IIMV','IIMV1']
    letter = ['a','b']
    file = ['rx_fixed','rx_free']
    for i = 0,1 do begin
        re = execute('model = '+models[i])
        re = execute('index = '+indices[i])
        ;; models
        rxm = hist2d_avg(model,0.2d,iidet=index,normalize=1,confidence=1)
        ind = minmax(where(rxm.yh_det gt 0.))
        ind = [ind[0]:ind[1]:1]
        ;; data
        yhd = histogram(rxd,bin=0.2d,location=xhd,min=min(rxm.xh),max=max(rxm.xh))
        indd = minmax(where(yhd gt 0.))
        indd = [indd[0]:indd[1]:1]
        ;; arbitrary normalization of data to match model peak
        ipeak = where(rxm.yh_det eq max(rxm[where(rxm.xh gt -1)].yh_det))
        normd = (rxm[ipeak].yh_det/yhd[ipeak])[0]

        e = {xra:[-2.6,1.],yra:[0.,0.20],$
             stairstep:1,fill_background:1, $
             xtitle:'$!8R!7_{!8L!7_X}$',ytitle:'Frequency (normalized)', $
             font_name:'Times',font_size:16, $
             buffer:0}
        if keyword_set(hide) then e.buffer = 1

        col1 = [4,90,141]
        col2 = [54,144,192]
        col3 = [116,169,207]
    
        p3u = plot(rxm[ind].xhoff,rxm[ind].sig3u_det,_extra=e,color='black',transparency=75,fill_color=col3)
        p2u = plot(rxm[ind].xhoff,rxm[ind].sig2u_det,_extra=e,color='black',transparency=75,fill_color=col2,/ov)
        p1u = plot(rxm[ind].xhoff,rxm[ind].sig1u_det,_extra=e,color='black',transparency=75,fill_color=col1,/ov)
        p1l = plot(rxm[ind].xhoff,rxm[ind].sig1l_det,_extra=e,color='black',transparency=75,fill_color=col2,/ov)
        p2l = plot(rxm[ind].xhoff,rxm[ind].sig2l_det,_extra=e,color='black',transparency=75,fill_color=col3,/ov)
        p3l = plot(rxm[ind].xhoff,rxm[ind].sig3l_det,_extra=e,color='black',transparency=75,fill_color='white',/ov)
        p = plot(rxm.xhoff,rxm.yh_det,/stairstep,thick=2,/ov,name='$!8R!7_X^{model}$')
        pd = plot(xhd[indd]+0.1,yhd[indd]*normd,'o',sym_size=1.,sym_filled=1,/ov,name='$!8R!7_{!8L!7_X} (data)$')
        plim = plot([1.,1.]#rx2nh(23.,/rx_out),e.yra,'--',thick=2,/ov)    
        psym3 = plot([0.,0.],[-1.,-1.],'s',sym_size=1.5,sym_filled=1,col='black',sym_fill_color=col3,name='3$\sigma$ conf.',/ov)
        psym2 = plot([0.,0.],[-1.,-1.],'s',sym_size=1.5,sym_filled=1,col='black',sym_fill_color=col2,name='2$\sigma$ conf.',/ov)
        psym1 = plot([0.,0.],[-1.,-1.],'s',sym_size=1.5,sym_filled=1,col='black',sym_fill_color=col1,name='1$\sigma$ conf.',/ov)
        leg = legend(target=[psym1,psym2,psym3,pd],position=[0.35,0.85],/normal,sample_width=0.,horizontal_spacing=0.06,font_name='Times',font_size=14)
        ;; NH limit
        t = text(0.63,0.8,/normal,'$!8N!7_H \geq$ 10$^{23}$ cm$^{-2}$',alignment=1.,font_name='Times',font_size=14)
        p = text(0.63,0.75,/normal,'$\leftarrow$',alignment=1.,font_size=26)
        ;; figure reference for Nature caption
        t = text(0.03,0.90,letter[i],/normal,font_name='Times',font_style='Bold',font_size=16)

        if keyword_set(sav) then begin
            print, '    SAVING PLOT'
            if (strupcase(strtrim(sav,2)) eq 'EPS') then p.save,'figures/'+file[i]+'.eps',/BITMAP else $
                                                         p.save,'figures/'+file[i]+'.png',resolution=res
        endif        
    endfor    

endif



END










