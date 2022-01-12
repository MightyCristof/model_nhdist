PRO xspec_nh_modeling


;; NH array
nh = [21.d:25.d:0.25d]
nnh = n_elements(nh)
nhstr = string(nh,format='(f5.2)')
nhlo = strtrim(10.^nhstr[0:3]/1e22,2)
nnhlo = n_elements(nhlo)
nhv = [nhlo,nhstr[4:-1]]

;; DIMENSIONS
;; 0 energy
;; 1 photon index
;; 2 scattering fracxtion
;; --skip for now 3 opening angle
;; 4 redshift
;; 5 compton reflection model (if using)

;; 0 energy
;; energy ranges and string prefix
eng = ['2 10','0.5 2','2 7']
;; 1 photon index
phot = string([1.4,1.6,1.8,2.0],format='(f3.1)')
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
fsig = [-2.:1.:0.5]
nscat = n_elements(fsig)
fsv = reform(string(fsm#10d^fsig,format='(d9.4)'),nnh,n_elements(fsig))

;; 3 redshift
;; redshift bins for model
zm = string([0.:0.8:0.01],format='(f4.2)')
nzm = n_elements(zm)

;; 4 opening angle
;; stuff for later

;; 5 compton reflection model
;; compton reflection
model = ['borus','buchner']
model = model[0]

;; XSPEC script setup
openw,1,'model_script_'+scat+'_'+model+'.xcm'
printf,1,'abund angr'
printf,1,'cosmo 70 0 0.73'
printf,1,'method leven 10 0.01'
printf,1,'xsect vern'
printf,1,'xset delta 0.01'

;; for the above dimensions to work out properly with REFORM() down the line,
;; you must loop over everything *backwards* in the next step.

;; loop over energy...
for ie = 0,n_elements(eng)-1 do begin
    printf,1,'set id [open model_flux'+strcompress(eng[ie],/rem)+'kev_'+scat+'_'+model+'.dat a]'
    ;; if 2-10 keV, restframe only
    if (eng[ie] eq '2 10') then z = zm[0] else z = zm
    ;; ...redshift (for 0.5-2 and 2-7 keV conversions only)
    for iz = 0,n_elements(z)-1 do begin
        ;; ... scattering fraction
        for is = 0,nscat-1 do begin
            ;; flag fscat variance
            printf,1,'puts $id "'+strtrim(fsig[is],2)+' -1 fscat"'
            ;; ...photon index
            for ip = 0,nphot-1 do begin
                ;; flag photon index change
                printf,1,'puts $id "'+strtrim(phot[ip],2)+' 0.05 gamma"'
                ;; modeling begins
                printf,1,'model  constant*phabs*zphabs(atable{/Users/ccarroll/heasoft-6.25/localmodels/borus/borus02_v170323a.fits} + zphabs*cabs*cutoffpl + constant*cutoffpl)'
                printf,1,'              1         -1          0          0      1e+10      1e+10'   
                printf,1,'           0.01         -1          0          0     100000      1e+06'
                printf,1,'           0.01         -1          0          0     100000      1e+06'
                printf,1,'           '+z[0]+'      -0.01     -0.999     -0.999         10         10'
                printf,1,'            '+phot[ip]+'       0.05        1.4       1.45       2.55        2.6'
                printf,1,'            300         -1         20         50        500       2000'
                printf,1,'             22        0.1         22      22.25         25       25.5'
                printf,1,'/'
                printf,1,'           0.05         -1'
                printf,1,'              1         -1'
                printf,1,'= p4'
                printf,1,'              1       0.01          0          0      1e+20      1e+24'
                printf,1,'              1      0.001          0          0     100000      1e+06'
                printf,1,'= p4'
                printf,1,'= p13'
                printf,1,'= p5'
                printf,1,'= p6'
                printf,1,'= p12'
                printf,1,'      '+fsv[0,is]+'         -1          0          0      1e+10      1e+10'
                printf,1,'= p5'
                printf,1,'= p6'
                printf,1,'= p12'
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
                ;; flag redshift change 
                printf,1,'new 4 '+z[iz]
                printf,1,'puts $id "[tcloutr param 4]"'
                for n = 0,nnhlo-1 do begin
                    printf,1,'new 13 '+nhv[n]
                    printf,1,'new 19 '+fsv[n,is]
                    printf,1,'flux '+eng[ie]
                    printf,1,'puts $id "[tcloutr flux 2]"'
                endfor
                ;; link obscuration to torus at log NH = 22
                printf,1,'new 13 = 1.0e-22*10.0^p7'
                for n = nnhlo,nnh-1 do begin
                    printf,1,'new 7 '+nhv[n]
                    printf,1,'new 19 '+fsv[n,is]
                    printf,1,'flux '+eng[ie]
                    printf,1,'puts $id "[tcloutr flux 2]"'
                endfor      
            ;; redshift
            endfor
        ;; scattering fraction
        endfor
    ;; photon index          
    endfor
    printf,1,'close $id'
;; energy
endfor
close,1


;; STOP AND RUN XSPEC
;; ctrl+c to continue
stop


;; DIMENSIONS
;; 0 energy
;; 1 photon index
;; 2 scattering fracxtion
;; --skip for now 3 opening angle
;; 4 redshift
;; 5 compton reflection model (if using)
file = 'model_flux210kev_'+scat+'_borus.dat'
readcol,file,full,dfull,sfull,format='d,d,a'
;; index dimensions
;ip = where(strmatch(supp,'gamma',/fold),plen)
;is = where(strmatch(supp,'fscat',/fold),slen)
;ia = where(dval eq 1.0,alen)
;iz = where(dval eq -0.01d,zlen)
;; reform flx_full values only    
ifull = where(dfull eq 0.0d,nfull)
if (nfull eq 0) then stop
flx_full = reform(val[inh],nnh,nphot,nscat)
;; create restframe 2-10 keV flx_full ratio
flx0 = rebin(flx_full[0,*,*],nnh,nphot,nscat)    
rx = alog10(flx_full/flx0)

;; observed frame conversion to 0.5-2 kev
file = 'model_flux0.52kev_mc_borus.dat'
readcol,file,val,dval,supp,format='d,d,a'
;; index dimensions
;ip = where(strmatch(supp,'gamma',/fold),plen)
;is = where(strmatch(supp,'fscat',/fold),slen)
;ia = where(dval eq 1.0,alen)
;iz = where(dval eq -0.01d,zlen)
;; reform flx_full values only    
inh = where(dval eq 0.0d,nval)
flx_soft = reform(val[inh],nnh,nphot,nscat,nzm)    

;; observed frame conversion to 2-7 kev
file = 'model_flux27kev_gupta_borus.dat'
readcol,file,val,dval,supp,format='d,d,a'
;; index dimensions
;ip = where(strmatch(supp,'gamma',/fold),plen)
;is = where(strmatch(supp,'fscat',/fold),slen)
;ia = where(dval eq 1.0,alen)
;iz = where(dval eq -0.01d,zlen)
;; reform flx_full values only    
inh = where(dval eq 0.0d,nval)
flx_hard = reform(val[inh],nnh,nphot,nscat,nzm)

;; make conversion arrays
c_soft = dblarr(size(flx_soft,/dim))
c_hard = dblarr(size(flx_hard,/dim))
for i = 0,zlen/slen-1 do begin
    c_soft[*,*,*,i] = flx_soft[*,*,*,i]/flx_full
    c_hard[*,*,*,i] = flx_hard[*,*,*,i]/flx_full
endfor






save,zm,nh,fsv, $
     rx, $
     flx_full,flx_soft,flx_hard, $
     c_soft,c_hard, $
     rx_fine, $
     c_hard_fine,c_soft_fine, $
     file='rx_'+scat+'_'+model+'.sav'


END







