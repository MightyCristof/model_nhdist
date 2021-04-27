

model_nhdist,'variant_models/lan_ad',setnh='nh_lan_nst',/setrx,group='wac',test='AD',/fixed
popd
model_nhdist,'variant_models/lan_ad_postmod',setnh='nh_lan_nst',/setrx,group='wac',test='AD',/postmod,/fixed
popd
model_nhdist,'variant_models/lan_joint',setnh='nh_lan_nst',/setrx,group='wac',test='JOINT',/fixed
popd
model_nhdist,'variant_models/lan_joint_postmod',setnh='nh_lan_nst',/setrx,group='wac',test='JOINT',/postmod,/fixed
popd
model_nhdist,'variant_models/lan_all_ad',setnh='nh_lan_nst',/setrx,group='all',test='AD',/fixed
popd
model_nhdist,'variant_models/lan_all_ad_postmod',setnh='nh_lan_nst',/setrx,group='all',test='AD',/postmod,/fixed
popd
model_nhdist,'variant_models/lan_all_joint',setnh='nh_lan_nst',/setrx,group='all',test='JOINT',/fixed
popd
model_nhdist,'variant_models/lan_all_joint_postmod',setnh='nh_lan_nst',/setrx,group='all',test='JOINT',/postmod,/fixed
popd




END
