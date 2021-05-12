

;; WISE AGN AD TEST
;model_nhdist,'variant_models/lan_ad_wac',setnh='nh_lan_nst',/setrx,group='WAC',test='AD',/fixed,/free,/split,/model
;popd
;file_copy,'variant_models/lan_ad_wac','variant_models/lan_ad_wac_pm',/recursive,/overwrite
;pushd,'variant_models/lan_ad_wac_pm'
;model_rxdist,/postmod
;popd

;; ALL SOURCES AD TEST
;model_nhdist,'variant_models/lan_ad_all',setnh='nh_lan_nst',/setrx,group='ALL',test='AD',/fixed,/free,/split,/model
;popd
;file_copy,'variant_models/lan_ad_all','variant_models/lan_ad_all_pm',/recursive,/overwrite
;pushd,'variant_models/lan_ad_all_pm'
;model_rxdist,/postmod
;popd

;; WISE AGN JOINT TEST
model_nhdist,'variant_models/lan_joint_wac',setnh='nh_lan_nst',/setrx,group='WAC',test='JOINT',/fixed,/free,/split,/model
popd
file_copy,'variant_models/lan_joint_wac','variant_models/lan_joint_wac_pm',/recursive,/overwrite
pushd,'variant_models/lan_joint_wac_pm'
model_rxdist,/postmod
popd


;; ALL SOURCES JOINT TEST
;model_nhdist,'variant_models/lan_joint_all',setnh='nh_lan_nst',/setrx,group='ALL',test='JOINT',/fixed,/free,/split,/model
;popd
;file_copy,'variant_models/lan_joint_all','variant_models/lan_joint_all_pm',/recursive,/overwrite
;pushd,'variant_models/lan_joint_all_pm'
;model_rxdist,/postmod
;popd



END






