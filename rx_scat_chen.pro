PRO rx_scat_chen, shortest = shortest, chen = chen



if keyword_set(chen) then begin
    file = '~/Research/collaborations/yan+21/sed_ages.fits'
    ;; AGES
    r = mrdfits(file,1)
    ii = where(r.lir gt 0. and r.lx gt 0.)
    lir = r[ii].lir
    lx = r[ii].lx
endif else begin
    files = '/Users/ccarroll/Research/projects/xray_lack_agn/workspace/run_20201008_final_updated/'+['src_luminosities.sav','combined_lum.sav']
    ;; Carroll+21
    for i = 0,1 do restore,files[i]
    ii = where(loglir ne -9999. and loglx ne -9999.)
    lir = loglir[ii]
    lx = loglx[ii]
endelse

;; relationship
xrel = [40.:50.:0.01]
yrel = dblarr(n_elements(xrel))
iichen = xrel lt 44.79
yrel[where(iichen)] = 0.84*(xrel[where(iichen)]-45.)+44.60
yrel[where(~iichen)] = 0.40*(xrel[where(~iichen)]-45.)+44.51


if keyword_set(shortest) then begin

    ;; or shortest distance to line
    x = lir
    y = lx
    mm = minmax(x)
    ii = where(xrel ge mm[0] and xrel le mm[1])
    xr = xrel[ii]
    yr = yrel[ii]
    rx = dblarr(n_elements(x))
    ry = dblarr(n_elements(x))
    d = dblarr(n_elements(x))
    ind = dblarr(n_elements(x))
    for i = 0,n_elements(x)-1 do begin
        d[i] = min(sqrt((xr-x[i])^2. + (yr-y[i])^2.),imin)
        rx[i] = xr[imin]-x[i]
        ry[i] = yr[imin]-y[i]
        ind[i] = imin
    endfor
    dsig = stddev(d)
    dd = [-d,d]
    plothist, dd
    ddsig = stddev(dd)
    print, ddsig
     ;0.29251033
     
    p = plot(x,y,'.')
    p = plot(xr,yr,'-r',/ov)
    p = plot(xr-dsig,yr+dsig,'.b',/ov)
    p = plot(xr+dsig,yr-dsig,'.b',/ov)

stop


endif else begin



    ;; plot
    p = plot(lir, lx, '.', xra=[41.,47.], yra=[41.,47.],/ov)
    p = plot(xrel, yrel, thick=2, col='red',/ov)

    ;; sync axes
    ;; what is RELATIONSHIP LX for a given LIR
    rlx = interpol(yrel,xrel,lir)
    ;; what is RELATIONSHIP LIR for a given LX
    rlir = interpol(xrel,yrel,lx)

    ;; scatter
    sig_lx = stddev(sqrt((rlx-lx)^2.))
    sig_lir = stddev(sqrt((rlir-lir)^2.))
    sig_rel = stddev(sqrt((rlx-lx)^2. + (rlir-lir)^2.))

    ;; test lx scatter by bins of LIR
    lirb = [40.:48.:1.]
    nlx = lonarr(n_elements(lirb))
    dlx = dblarr(n_elements(lirb))

    for i = 0,n_elements(lirb)-2 do begin
        ig = where(lir gt lirb[i] and lir le lirb[i+1],ct)
        if (ct eq 0) then continue
        nlx[i] = ct
        dlx[i] = stddev(rlx[ig]-lx[ig])
    endfor

    print, sig_lx
     ;0.22540535
    print, sig_lir
     ;0.35255237
    print, sig_rel
     ;0.41288108

endelse







END