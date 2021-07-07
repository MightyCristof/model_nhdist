FUNCTION fscat_erf, nh_in, $
                    PLT = plt, $
                    LOW = low


;; pull data
readcol,'Brightman*.csv',nh,fs,e_fs,e_nh
;; add bounds
if keyword_set(low) then begin
    nh_pref = 24.
    npref = n_elements(nh_pref)
    nh_suff = [24.5,25.0,25.5]
    nsuff = n_elements(nh_suff)
    nhv = [nh_pref,nh[-1],nh_suff]
    e_nhv = [dblarr(npref)+e_nh[-1],e_nh[-1],dblarr(nsuff)+e_nh[-1]]
    fsv = [dblarr(npref)+0.01,fs[-1],0.008,0.0065,0.005]
    e_fsv = [dblarr(npref)+e_fs[-1],e_fs[-1],dblarr(nsuff)+e_fs[-1]]
endif else begin
    nh_pref = 22.;[21.:nh[0]:1.]
    npref = n_elements(nh_pref)
    nh_suff = 25.5;reverse([25.5:nh[-1]:-1.])
    nsuff = n_elements(nh_suff)
    nhv = [nh_pref,nh,nh_suff]
    e_nhv = [dblarr(npref)+median(e_nh,/even),e_nh,dblarr(nsuff)+median(e_nh,/even)]
    fsv = [dblarr(npref)+0.08,fs,dblarr(nsuff)+0.005]
    e_fsv = [dblarr(npref)+median(e_fs,/even)/5.,e_fs,dblarr(nsuff)+median(e_fs,/even)/5.]
endelse

;; FUNCTION: ERF
x = [-5.:5.:0.1]
y = erfc(x)
c1 = [0.05:2.:0.05]
c2 = [22.:24.:0.1]
c3 = [0.01:0.10:0.005]
c4 = [-0.01:0.01:0.001]

;; FUNCTION: 1/x
;x = [-9.9d:9.9d:0.2d]
;y = 1./x
;c1 = [2.:8.:0.05]
;c2 = [21.:23.:0.1]
;c3 = [0.001:0.02:0.001]
;c4 = [-0.05:0.05:0.001]

nx = n_elements(x)
nc1 = n_elements(c1)
nc2 = n_elements(c2)
nc3 = n_elements(c3)
nc4 = n_elements(c4)

xx = dblarr(nx,nc1,nc2,nc3,nc4)
yy = dblarr(nx,nc1,nc2,nc3,nc4)

for i = 0,nc1-1 do begin
    for j = 0,nc2-1 do begin
        for k = 0,nc3-1 do begin
            for l = 0,nc4-1 do begin
                xx[*,i,j,k,l] = c1[i]*x+c2[j]
                yy[*,i,j,k,l] = c3[k]*y+c4[l]
            endfor
        endfor
    endfor
endfor

chiy = dblarr(nc1,nc2,nc3,nc4)
for i = 0,nc1-1 do begin
    for j = 0,nc2-1 do begin
        for k = 0,nc3-1 do begin
            for l = 0,nc4-1 do begin
                chiy[i,j,k,l] = ndx2(fsv,e_fsv,interpol(yy[*,i,j,k,l],xx[*,i,j,k,l],nhv))
            endfor
        endfor
    endfor
endfor

chiym = min(chiy,imin)
ind = array_indices(chiy,imin)

if keyword_set(plt) then begin
    p = errorplot(nhv,fsv,e_nhv,e_fsv,'o',sym_thick=4,xra=[21.,26.],yra=[-0.005,0.105])
    p = errorplot(nh,fs,e_nh,e_fs,'ob',sym_thick=4,sym_filled=1,/ov)
    p = plot(xx[*,ind[0],ind[1],ind[2],ind[3]],yy[*,ind[0],ind[1],ind[2],ind[3]],col='red',/ov)
    p.save,'vscat.png'
endif

fs_out = interpol(yy[*,ind[0],ind[1],ind[2],ind[3]],xx[*,ind[0],ind[1],ind[2],ind[3]],nh_in)
return, fs_out


END


