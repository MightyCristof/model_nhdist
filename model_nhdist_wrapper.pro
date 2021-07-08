
;model_nhdist,/data

dir = 'final_model/'
;;; WISE AGN AD TEST
model_nhdist,dir+'lan_ad_wac_pm',setnh='nh_lan_nst',/setrx,group='WAC',test='AD',/fixed,/free,/model,/postmod
;popd
;file_copy,dir+'lan_ad_wac_mc',dir+'lan_ad_wac_pm_mc',/recursive,/overwrite
;pushd,dir+'lan_ad_wac_pm'
;model_rxdist,/postmod


;dir = 'variant_models/plus_c17/scat01/'
;; WISE AGN AD TEST
;model_nhdist,dir+'lan_ad_all',setnh='nh_lan_nst',/setrx,group='ALL',test='AD',/fixed,/free,/split,/model
;popd
;file_copy,dir+'lan_ad_all',dir+'lan_ad_all_pm',/recursive,/overwrite
;pushd,dir+'lan_ad_all_pm'
;model_rxdist,/postmod

;; ALL SOURCES AD TEST
;model_nhdist,'variant_models/lan_ad_all',setnh='nh_lan_nst',/setrx,group='ALL',test='AD',/fixed,/free,/split,/model
;popd
;file_copy,'variant_models/lan_ad_all','variant_models/lan_ad_all_pm',/recursive,/overwrite
;pushd,'variant_models/lan_ad_all_pm'
;model_rxdist,/postmod
;popd

;; WISE AGN JOINT TEST
;model_nhdist,'variant_models/lan_joint_wac',setnh='nh_lan_nst',/setrx,group='WAC',test='JOINT',/fixed,/free,/split,/model
;popd
;file_copy,'variant_models/lan_joint_wac','variant_models/lan_joint_wac_pm',/recursive,/overwrite
;pushd,'variant_models/lan_joint_wac_pm'
;model_rxdist,/postmod
;popd

;; ALL SOURCES JOINT TEST
;model_nhdist,'variant_models/lan_joint_all',setnh='nh_lan_nst',/setrx,group='ALL',test='JOINT',/fixed,/free,/split,/model
;popd
;file_copy,'variant_models/lan_joint_all','variant_models/lan_joint_all_pm',/recursive,/overwrite
;pushd,'variant_models/lan_joint_all_pm'
;model_rxdist,/postmod
;popd


;if keyword_set(batch) then begin
;;; BATCH!!
;
;;; WISE AGN JOINT TEST - change directory
;for i = 2,10 do begin
;    model_nhdist,'variant_models/lan_ad_all',/free
;    file_move,'ctf_free.sav','ctf_free'+string(i,format='(i02)')+'.sav'
;    popd
;endfor
;file = file_search('ctf_free*')
;vars = ['FCTV1','F24V1','F25V1','STATV1','NREJ1','NHMV1','RXMV1','IIMV1']
;for i = 0,n_elements(file)-1 do begin
;    restore,file[i]
;    if (i eq 0) then begin
;        for j = 0,n_elements(vars)-1 do re = execute('temp_'+vars[j]+' = '+vars[j])
;    endif else begin
;        temp_FCTV1  = [temp_FCTV1   ,FCTV1   ]
;        temp_F24V1  = [temp_F24V1   ,F24V1   ]
;        temp_F25V1  = [temp_F25V1   ,F25V1   ]
;        temp_STATV1 = [[temp_STATV1],[STATV1]]
;        temp_NREJ1  = [temp_NREJ1   ,NREJ1   ]
;        temp_NHMV1  = [[temp_NHMV1] ,[NHMV1] ]
;        temp_RXMV1  = [[temp_RXMV1] ,[RXMV1] ]
;        temp_IIMV1  = [[temp_IIMV1] ,[IIMV1] ]
;    endelse
;endfor
;for j = 0,n_elements(vars)-1 do re = execute(vars[j]+' = temp_'+vars[j])
;sav_str = strjoin(vars,',')
;re = execute('save,'+sav_str+',file="ctf_free.sav"')
;
;;; WISE AGN JOINT TEST - change directory
;for i = 2,10 do begin
;    model_nhdist,'variant_models/lan_ad_all',/split
;    file_move,'nh_split.sav','nh_split'+string(i,format='(i02)')+'.sav'
;    popd
;endfor
;file = file_search('nh_split*')
;vars = ['F24V2','F25V2','STATV2','NREJ2','NHMV2','RXMV2','IIMV2']
;for i = 0,n_elements(file)-1 do begin
;    restore,file[i]
;    if (i eq 0) then begin
;        for j = 0,n_elements(vars)-1 do re = execute('temp_'+vars[j]+' = '+vars[j])
;    endif else begin
;        temp_F24V2  = [temp_F24V2   ,F24V2   ]
;        temp_F25V2  = [temp_F25V2   ,F25V2   ]
;        temp_STATV2 = [[temp_STATV2],[STATV2]]
;        temp_NREJ2  = [temp_NREJ2   ,NREJ2   ]
;        temp_NHMV2  = [[temp_NHMV2] ,[NHMV2] ]
;        temp_RXMV2  = [[temp_RXMV2] ,[RXMV2] ]
;        temp_IIMV2  = [[temp_IIMV2] ,[IIMV2] ]
;    endelse
;endfor
;for j = 0,n_elements(vars)-1 do re = execute(vars[j]+' = temp_'+vars[j])
;sav_str = strjoin(vars,',')
;re = execute('save,'+sav_str+',file="nh_split.sav"')
;endif


END






