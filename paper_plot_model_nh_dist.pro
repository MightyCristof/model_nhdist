PRO paper_plot_model_nh_dist, LL_DIST = ll_dist, $
                              NH_DIST = nh_dist, $
                              SHOW = show, $
                              SAV = sav
                               


;; load data
common _fits    
common _inf_cha 
common _inf_xmm 
common _inf_nst 
common _det_cha 
common _det_xmm 
common _det_nst 
common _wac 
common _softx   
common _fluxlim 
common _comp    
common _agnlum 
common _clean_cha
common _clean_xmm
common _clean_nst
common _quality  
common _lum_ratio
common _nhdist


;;----------------------------------------------------------------------------------------
;; LUMINOSITY RATIO DISTRIBUTION
;;----------------------------------------------------------------------------------------
if keyword_set(ll_dist) then begin
    print, '    PREPARING PAPER LUMINOSITY RATIO DISTRIBUTION'

    e = {, $
         buffer:1}
    if keyword_set(show) then e.buffer = 0
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strtrim(sav,2) eq 'EPS') then p.save,'ll_distribution.eps',/BITMAP else $
                                          p.save,'ll_distribution.png'
    endif
endif


;;----------------------------------------------------------------------------------------
;; NH DISTRIBUTION
;;----------------------------------------------------------------------------------------
if keyword_set(nh_dist) then begin
    print, '    PREPARING PAPER NH DISTRIBUTION'

    e = {, $
         buffer:1}
    if keyword_set(show) then e.buffer = 0
    
    if keyword_set(sav) then begin
        print, '    SAVING PLOT'
        if (strtrim(sav,2) eq 'EPS') then p.save,'nh_distribution.eps',/BITMAP else $
                                          p.save,'nh_distribution.png'
    endif
endif


END








    ;; MODEL L/L COMPARISON TO LANSBURY+15
    e = {xtitle:'$log( L_{10-40 keV} / L_{10-40 keV}(L_{IR}) )$',peak:1,buffer:1}
    h0 = histplot(xllsamp,hllsamp,_extra=e,name='Sampled Lansbury+15')
    h1 = histplot(xlldet_hb,hlldet_hb,e_hlldet_hb,_extra=e,sym='o',col='dodger blue',ov=h0,name='AGN + X-ray det.')
    h2 = histplot(xlllim_nox,hlllim_nox,_extra=e,col='orange',ov=h0,name='AGN X-ray non-det.')
    ;h3 = histplot(xllmod,hllmod,/peak,col='teal',ov=h0,name='$Corrected L15 model$')
    leg = legend(target=[h0,h1,h2],position=[abs(width(h0.xra))/2.+h0.xra[0],h0.yra[1]*0.95],/data,/auto_text_color,sample_width=0)
    h0.save,'model_ll_nox.png'

    ;;  MODEL NH COMPARISON TO LANSBURY+15
    e = {xtitle:'$Log(N_H/cm^{-2})$',peak:1,buffer:1}
    h0 = histplot(xlan,hlan,_extra=e,name='Lansbury+15')
    h1 = histplot(xnhdet_hb,hnhdet_hb,_extra=e,col='dodger blue',ov=h0,name='AGN + X-ray det.')
    h2 = histplot(xnhlim_nox,hnhlim_nox,_extra=e,col='orange',ov=h0,name='AGN X-ray non-det.')
    ;h3 = histplot(xnhmod,hnhmod,/peak,col='teal',ov=h0,name='$Corrected N_H model$')    
    leg = legend(target=[h0,h1,h2],position=[width(h0.xra)/2.+h0.xra[0],h0.yra[1]*0.95],/data,/auto_text_color,sample_width=0)
    h0.save,'model_nh_nox.png'





