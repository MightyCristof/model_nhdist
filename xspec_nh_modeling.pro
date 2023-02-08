PRO xspec_nh_modeling, XSPEC = xspec, $
                       CONV = conv, $
                       RUN = run


;; NH array
nh = [21.d:25.d:0.25d]
nnh = n_elements(nh)
nhstr = string(nh,format='(f5.2)')
nhlo = strtrim(10.^nhstr[0:3]/1e22,2)
nnhlo = n_elements(nhlo)
nhv = [nhlo,nhstr[4:-1]]

;; DIMENSIONS
;; (looped in reverse order)
;; 0 energy
;; 1 photon index
;; 2 scattering fracxtion
;; 3 reflection strength
;; 4 opening angle
;; 5 redshift
;; 6 compton reflection model (if using)

;; 0 energy
;; energy ranges and string prefix
eng = ['2 10','0.5 2','2 7']
neng = n_elements(eng)

;; 1 photon index
phot = string([1.8,1.9,2.0],format='(f3.1)')
nphot = n_elements(phot)

;; 2 scattering fraction
;scat = 'LR'
scat = 'mc'
case scat of
    'lr': readcol,'Gupta+2021_Fig01_model_LR.csv',nhr,fsr,format='d,d'
    'mc': readcol,'Gupta+2021_Fig01_model_MC.csv',nhm,fsm,format='d,d'
endcase 
;; interpolate, convert to decimal, string
fsm = 10.^interpol(alog10(fsm),nhm,nh)
fsig = [-2.:0.5:0.5d]
nscat = n_elements(fsig)
fsv = reform(string(fsm#10d^fsig,format='(d9.4)'),nnh,n_elements(fsig))
fsig = string(fsig,format='(f4.1)')

;; 3 reflection strength
refl = [0.:2.:0.2d]
rrefl = string(refl,format='(f3.1)')
nrefl = n_elements(refl)

;; 4 opening angle
;; stuff for later
oang = [10.:80.:10.d]
oadeg = string(oang,format='(f4.1)')
oarad = strtrim(oang*!DTOR,2)
cosoa = strtrim(cos(oang*!DTOR),2)
noang = n_elements(oang)

;; 5 redshift
;; redshift bins for model
zm = string([0.:0.8:0.02d],format='(f4.2)')
nzm = n_elements(zm)

;; 6 compton reflection model
;; compton reflection
model = ['borus','buchner']
model = model[0]

if keyword_set(xspec) then begin    
    ;; for the above dimensions to work out properly with REFORM() down the line,
    ;; you must loop over everything *backwards* in the next step.

    ;; loop over...
    ;;; ...energy
    ;for ie = 0,n_elements(eng)-1 do begin
    ;    ;; if 2-10 keV, restframe only
    ;    if (eng[ie] eq '2 10') then z = zm[0] else z = zm
    ;    nz = n_elements(z)

        ;; XSPEC script setup
        ;openw,1,'model_script_'+scat+'_'+model+'_'+strcompress(eng[ie],/rem)+'kev.xcm'
        openw,1,'model_script_'+scat+'_'+model+'.xcm'
        printf,1,'abund angr'
        printf,1,'cosmo 70 0 0.73'
        printf,1,'method leven 10 0.01'
        printf,1,'xsect vern'
        printf,1,'xset delta 0.01'    
        printf,1,'set id [open model_flux_'+scat+'_'+model+'.dat a]'
        ;; modeling begins
        printf,1,'model  constant*phabs*zphabs(atable{/Users/ccarroll/heasoft-6.25/localmodels/borus/borus02_v170323a.fits} + zphabs*cabs*pexmon + constant*cutoffpl)'
        printf,1,'              1         -1          0          0      1e+10      1e+10'   
        printf,1,'           0.01         -1          0          0     100000      1e+06'
        printf,1,'           0.01         -1          0          0     100000      1e+06'
        printf,1,'           '+zm[0]+'      -0.01     -0.999     -0.999         10         10'
        printf,1,'            '+phot[0]+'       0.05        1.4       1.45       2.55        2.6'
        printf,1,'            300         -1         20         50        500       2000'
        printf,1,'             22        0.1         22      22.25         25       25.5'
        printf,1,'       '+cosoa[0]+'         -1'
        printf,1,'      0.0871557         -1'
        printf,1,'              1         -1'
        printf,1,'= p4'
        printf,1,'              1       0.01          0          0      1e+20      1e+24'
        printf,1,'              1      0.001          0          0     100000      1e+06'
        printf,1,'= p4'
        printf,1,'= p13'
        printf,1,'= p5'
        printf,1,'= p6'
        printf,1,'            '+rrefl[0]+'         -1'
        printf,1,'= p4'
        printf,1,'/'
        printf,1,'/'
        printf,1,'             85         -1'                        
        printf,1,'= p12'
        printf,1,'      '+fsv[0,0]+'         -1          0          0      1e+10      1e+10'
        printf,1,'= p5'
        printf,1,'= p6'
        printf,1,'= p12'
        ;;========================================================================
        ;;Model Model Component  Parameter  Unit     Value
        ;; par  comp
        ;;   1    1   constant   factor              1.00000      frozen
        ;;   2    2   phabs      nH         10^22    1.00000E-02  frozen
        ;;   3    3   zphabs     nH         10^22    1.00000E-02  frozen
        ;;   4    3   zphabs     Redshift            0.0          frozen
        ;;   5    4   borus02    PhoIndex            1.80000      +/-  0.0
        ;;   6    4   borus02    Ecut       keV      300.000      frozen
        ;;   7    4   borus02    logNHtor            22.0000      +/-  0.0
        ;;   8    4   borus02    CFtor               0.174533     frozen
        ;;   9    4   borus02    cos(thInc)          8.71557E-02  frozen      (85 degrees)
        ;;  10    4   borus02    A_Fe       A_Sun    1.00000      frozen
        ;;  11    4   borus02    z                   0.0          = p4
        ;;  12    4   borus02    norm                1.00000      +/-  0.0
        ;;  13    5   zphabs     nH         10^22    1.00000      +/-  0.0
        ;;  14    5   zphabs     Redshift            0.0          = p4
        ;;  15    6   cabs       nH         10^22    1.00000      = p13
        ;;  16    7   pexmon     PhoIndex            1.80000      = p5
        ;;  17    7   pexmon     foldE      keV      300.000      = p6
        ;;  18    7   pexmon     rel_refl            0.0          frozen
        ;;  19    7   pexmon     redshift            0.0          = p4
        ;;  20    7   pexmon     abund               1.00000      frozen
        ;;  21    7   pexmon     Fe_abund            1.00000      frozen
        ;;  22    7   pexmon     Incl       deg      85.0000      frozen
        ;;  23    7   pexmon     norm                1.00000      = p12
        ;;  24    8   constant   factor              0.102700     frozen
        ;;  25    9   cutoffpl   PhoIndex            1.80000      = p5
        ;;  26    9   cutoffpl   HighECut   keV      300.000      = p6
        ;;  27    9   cutoffpl   norm                1.00000      = p12
        ;;________________________________________________________________________
       
        ;; ...redshift (for 0.5-2 and 2-7 keV conversions only)
        ;; par 4
        for iz = 0,nzm-1 do begin
            printf,1,'newpar 4 '+zm[iz]
            printf,1,'puts $id "'+zm[iz]+' -1 redshift"'
            
            ;; opening angle
            ;; par 8
            for ia = 0,noang-1 do begin
                printf,1,'newpar 8 '+cosoa[ia]
                printf,1,'puts $id "'+oadeg[ia]+' -1 angle"'
                
                ;; reflection strength   
                ;; par 18 
                for ir = 0,nrefl-1 do begin
                    printf,1,'newpar 18 '+rrefl[ir]
                    printf,1,'puts $id "'+rrefl[ir]+' -1 refl"'
                    
                    ;; ... scattering fraction variance
                    ;; par 24
                    for is = 0,nscat-1 do begin
                        ;; fsv changes with NH (below), not here
                        printf,1,'puts $id "'+strtrim(fsig[is],2)+' -1 fscat"'

                        ;; ...photon index
                        ;; par 5
                        for ip = 0,nphot-1 do begin
                            printf,1,'newpar 5 '+phot[ip]
                            printf,1,'puts $id "'+strtrim(phot[ip],2)+' -1 gamma"'

                            ;; ... low NH
                            ;; low NH par 13
                            for n = 0,nnhlo-1 do begin
                                printf,1,'newpar 13 '+strtrim(nhv[n],2)
                                printf,1,'newpar 24 '+strtrim(fsv[n,is],2)
                                
                                ;; ... energy
                                ;; print 2-10 keV for z=0 only
                                if (iz eq 0) then begin
                                    printf,1,'flux '+eng[iz]
                                    printf,1,'puts $id "[tcloutr flux 2]"'
                                endif
                                for ie = 1,neng-1 do begin
                                    printf,1,'flux '+eng[ie]
                                    printf,1,'puts $id "[tcloutr flux 2]"'
                                ;; end energy
                                endfor
                            ;; end low NH                    
                            endfor
                            
                            ;; ... high NH
                            ;; high NH par 24
                            ;; link obscuration to torus at log NH = 22
                            printf,1,'newpar 13 = 1.0e-22*10.0^p7'
                            for n = nnhlo,nnh-1 do begin
                                printf,1,'newpar 7 '+strtrim(nhv[n],2)
                                printf,1,'newpar 24 '+strtrim(fsv[n,is],2)
                                
                                ;; ... energy
                                ;; print 2-10 keV for z=0 only
                                if (iz eq 0) then begin
                                    printf,1,'flux '+eng[iz]
                                    printf,1,'puts $id "[tcloutr flux 2]"'
                                endif
                                for ie = 1,neng-1 do begin
                                    printf,1,'flux '+eng[ie]
                                    printf,1,'puts $id "[tcloutr flux 2]"'
                                ;; end energy
                                endfor
                            ;; end high NH
                            endfor      
                        ;; end redshift
                        endfor
                    ;; end scattering fraction
                    endfor
                ;; end opening angle
                endfor
            ;; end reflection strength
            endfor
        ;; end photon index          
        endfor
        printf,1,'close $id'
        close,1
    ;;; energy
    ;endfor
    
    if keyword_set(run) then begin
        ;; XSPEC run file for entire project
        openw,1,'model_script_'+scat+'_'+model+'_run.xcm'
        for ie = 0,n_elements(eng)-1 do printf,1,'@model_script_'+scat+'_'+model+'_'+strcompress(eng[ie],/rem)+'kev.xcm'
        close,1
    endif
endif

;; STOP AND RUN XSPEC
;; ctrl+c to continue
;stop

if keyword_set(conv) then begin
    ;; DIMENSIONS
    ;; 0 energy
    ;; 1 photon index
    ;; 2 scattering fracxtion
    ;; 3 reflection strength
    ;; 4 opening angle
    ;; 5 redshift
    ;; 6 compton reflection model (if using)
    ;ffile = 'model_flux210kev_'+scat+'_'+model+'.dat'
    
    ffile = 'model_flux_'+scat+'_'+model+'.dat'
    readcol,ffile,modfx,flag,descr,format='d,d,a'
    iz = where(descr eq 'redshift',/null)
    re = execute('ifull = where(flag['+strjoin(strtrim(iz[0:1],2),':')+'] eq 0.0d,nfull)')
    if (nfull eq 0) then stop
    full = reform(modfx[ifull],neng,nnh,nphot,nscat,nrefl,noang,1)
    flx_full = reform(full[0,*,*,*,*,*])
    ;; create restframe 2-10 keV flx_modfx ratio
    flx0 = rebin(flx_full[0,*,*,*,*],nnh,nphot,nscat,nrefl,noang)
    rx = alog10(flx_full/flx0)
    
    ;; observed frame for soft 0.5-2 kev and hard 2-7 keV fluxes
    re = execute('irest = iz[1]+where(flag['+strjoin([strtrim(iz[1],2),"*"],":")+'] eq 0.0,nrest)')
    if (nrest eq 0) then stop
    rest = reform(modfx[irest],neng-1,nnh,nphot,nscat,nrefl,noang,nzm-1)
    flx_soft = dblarr(nnh,nphot,nscat,nrefl,noang,nzm)
    flx_soft[*,*,*,*,*,0] = reform(full[1,*,*,*,*,*,*])
    flx_soft[*,*,*,*,*,1:-1] = reform(rest[0,*,*,*,*,*,*])
    flx_hard = dblarr(nnh,nphot,nscat,nrefl,noang,nzm)
    flx_hard[*,*,*,*,*,0] = reform(full[2,*,*,*,*,*,*])
    flx_hard[*,*,*,*,*,1:-1] = reform(rest[1,*,*,*,*,*,*])
    ;; make conversion arrays
    flx0z = rebin(flx_full,nnh,nphot,nscat,nrefl,noang,nzm)
    c_soft = flx_soft/flx0z
    c_hard = flx_hard/flx0z

    ;;; fine grid (must loop due to interpolation)
    ;nh_fine = [nh[0]:nh[-1]:0.002d]
    ;rx_fine = dblarr([size(nh_fine,/dim),nphot,nscat])
    ;c_hard_fine = dblarr([size(rx_fine,/dim),nzm])
    ;c_soft_fine = dblarr([size(rx_fine,/dim),nzm])
    ;for ip = 0,nphot-1 do begin
    ;    for is = 0,nscat-1 do begin
    ;        rx_fine[*,ip,is] = interpol(rx[*,ip,is],nh,nh_fine)
    ;        for iz = 0,nzm-1 do begin
    ;            c_hard_fine[*,ip,is,iz] = interpol(c_hard[*,ip,is,iz],rx[*,ip,is],rx_fine[*,ip,is])
    ;            c_soft_fine[*,ip,is,iz] = interpol(c_soft[*,ip,is,iz],rx[*,ip,is],rx_fine[*,ip,is])
    ;        endfor
    ;    endfor
    ;endfor
    
    ;; output FITS file with header comments
    mkhdr,hdr,rx
    sxaddhist,['EXT 0: RX','EXT 1: SOFT FX CONVERSION','EXT 2: HARD FX CONVERSION','EXT 3: LOG NH','EXT 4: PHOTON INDEX','EXT 5: SIGMA FSCAT','EXT 6: R REFLECTION','EXT 7: OPENING ANGLE [DEGREES]','EXT 8: REDSHIFT'],hdr,/comment
    mwrfits,rx,'rx.fits',hdr,/create
    mwrfits,c_soft,'rx.fits'
    mwrfits,c_hard,'rx.fits'
    mwrfits,nh,'rx.fits'
    mwrfits,double(phot),'rx.fits'
    mwrfits,double(fsig),'rx.fits'
    mwrfits,refl,'rx.fits'
    mwrfits,oang,'rx.fits'
    mwrfits,double(zm),'rx.fits'
    
    save,phot,fsm,fsig,fsv,refl,oang,zm,eng, $
         flx_full,flx_soft,flx_hard, $
         nh, $;nh_fine, $
         rx, $; rx_fine, $
         c_soft,c_hard, $; c_hard_fine,c_soft_fine, $
         file='rx_'+scat+'_'+model+'.sav'

endif


END
























;; old modeling

                            ;printf,1,'model  constant*phabs*zphabs(atable{/Users/ccarroll/heasoft-6.25/localmodels/borus/borus02_v170323a.fits} + zphabs*cabs*cutoffpl + constant*cutoffpl)'
                            ;printf,1,'              1         -1          0          0      1e+10      1e+10'   
                            ;printf,1,'           0.01         -1          0          0     100000      1e+06'
                            ;printf,1,'           0.01         -1          0          0     100000      1e+06'
                            ;printf,1,'           '+z[0]+'      -0.01     -0.999     -0.999         10         10'
                            ;printf,1,'            '+phot[ip]+'       0.05        1.4       1.45       2.55        2.6'
                            ;printf,1,'            300         -1         20         50        500       2000'
                            ;printf,1,'             22        0.1         22      22.25         25       25.5'
                            ;printf,1,'/'
                            ;printf,1,'           0.05         -1'
                            ;printf,1,'              1         -1'
                            ;printf,1,'= p4'
                            ;printf,1,'              1       0.01          0          0      1e+20      1e+24'
                            ;printf,1,'              1      0.001          0          0     100000      1e+06'
                            ;printf,1,'= p4'
                            ;printf,1,'= p13'
                            ;printf,1,'= p5'
                            ;printf,1,'= p6'
                            ;printf,1,'= p12'
                            ;printf,1,'      '+fsv[0,is]+'         -1          0          0      1e+10      1e+10'
                            ;printf,1,'= p5'
                            ;printf,1,'= p6'
                            ;printf,1,'= p12'
                            ;;  Model Model Component  Parameter  Unit     Value
                            ;;   par  comp
                            ;;     1    1   constant   factor              1.00000      frozen
                            ;;     2    2   phabs      nH         10^22    1.00000E-02  frozen
                            ;;     3    3   zphabs     nH         10^22    1.00000E-02  frozen
                            ;;     4    3   zphabs     Redshift            0.0          frozen
                            ;;     5    4   borus02    PhoIndex            1.80000      +/-  0.0
                            ;;     6    4   borus02    Ecut       keV      300.000      frozen
                            ;;     7    4   borus02    logNHtor            22.0000      +/-  0.0
                            ;;     8    4   borus02    CFtor               0.500000     +/-  0.0
                            ;;     9    4   borus02    cos(thInc)          5.00000E-02  frozen
                            ;;    10    4   borus02    A_Fe       A_Sun    1.00000      frozen
                            ;;    11    4   borus02    z                   0.0          = p4
                            ;;    12    4   borus02    norm                1.00000      +/-  0.0
                            ;;    13    5   zphabs     nH         10^22    1.00000      +/-  0.0
                            ;;    14    5   zphabs     Redshift            0.0          = p4
                            ;;    15    6   cabs       nH         10^22    1.00000      = p13
                            ;;    16    7   cutoffpl   PhoIndex            1.80000      = p5
                            ;;    17    7   cutoffpl   HighECut   keV      300.000      = p6
                            ;;    18    7   cutoffpl   norm                1.00000      = p12
                            ;;    19    8   constant   factor              0.177200     +/-  0.0
                            ;;    20    9   cutoffpl   PhoIndex            1.80000      = p5
                            ;;    21    9   cutoffpl   HighECut   keV      300.000      = p6
                            ;;    22    9   cutoffpl   norm                1.00000      = p12   



