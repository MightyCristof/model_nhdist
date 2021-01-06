FUNCTION mcmc_testmodel,x,param,_EXTRA=ex
;a very simple model of NEE
  nee = (param[0] * exp(param[1]*(x[0,*]-15.0))) -(( param[2]*x[1,*]) / (param[3]+x[1,*]))
  return,nee
END