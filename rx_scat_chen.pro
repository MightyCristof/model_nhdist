PRO rx_scat_chen



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



;; or shortest distance to line
x = lir
y = lx
nobj = n_elements(x)
mm = minmax(x)
ii = where(xrel ge mm[0] and xrel le mm[1])
xr = xrel[ii]
yr = yrel[ii]
d = dblarr(nobj)
dx = dblarr(nobj)
dy = dblarr(nobj)
ind = dblarr(nobj)
for i = 0,nobj-1 do begin
    d[i] = min(sqrt((xr-x[i])^2. + (yr-y[i])^2.),imin)
    dx[i] = xr[imin]-x[i]
    dy[i] = yr[imin]-y[i]
    ind[i] = imin
endfor
;; combined uncertainty
dsig = sqrt(stddev(dx)^2. + stddev(dy)^2.)

;; simulate spread about zero
irand = randomi(nobj,nobj,/nodup)
dd = [-d[irand[0:nobj/2-1]],d[irand[nobj/2:-1]]]
plothist, dd
ddsig = stddev(dd)

print, dsig, ddsig
 
p = plot(x,y,'.')
p = plot(xr,yr,'-r',/ov)
p = plot(xr-dsig,yr+dsig,'.b',/ov)
p = plot(xr+dsig,yr-dsig,'.b',/ov)
p = plot(xr-ddsig,yr+ddsig,'.b',/ov)
p = plot(xr+ddsig,yr-ddsig,'.b',/ov)


END












