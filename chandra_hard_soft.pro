PRO chandra_hard_soft, sband, $
                       mband, $
                       hband, $
                       softx, $
                       hardx


nin = n_elements(sband[*,0])
softx = dblarr(nin,2)
hardx = dblarr(nin,2)

;; Flux conversions
;; WebPIMMS parameters: Galactic NH=2E20, power law photon index=1.8
;; https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3pimms/w3pimms.pl
;; Chandra *1.1 to convert 90% aperture flux to 100%
;; energy band (keV)
;; Chandra          
;; b  0.5-7         
;; h  2-7           
;; m  1.2-2         2.414E+00 ergs/cm/cm/s
;; s  0.5-1.2       1.707E+00 ergs/cm/cm/s
;; u  0.2-0.5       
;; w  0.1-10        

;; stack soft band x-rays according to least fractional conversion factor to 0.5-2 keV
;; chandra m
iim = mband[*,0] gt 0. and mband[*,1] gt 0. and finite(mband[*,0]) and finite(mband[*,1])
im = where(iim,mct)
if (mct gt 0) then softx[im,*] = mband[im,*]*2.414E+00*1.1
;; chandra s
iis = sband[*,0] gt 0. and sband[*,1] gt 0. and finite(sband[*,0]) and finite(sband[*,1])
is = where(iis,sct)
if (sct gt 0) then softx[is,*] = sband[is,*]*1.707E+00*1.1

;; chandra h
iih = hband[*,0] gt 0. and hband[*,1] gt 0. and finite(hband[*,0]) and finite(hband[*,1])
ih = where(iih,hct)
if (hct gt 0) then hardx[ih,*] = hband[ih,*]*1.1


END