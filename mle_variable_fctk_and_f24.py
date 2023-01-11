#!/usr/bin/env python
import sys
import os
sys.path.append(os.getcwd())
import numpy as np
#applies to all three functions
import gc
gc.enable()
from scipy.interpolate import interp1d, RegularGridInterpolator
from astropy.cosmology import WMAP9 as cosmo
import gzip
from astropy.cosmology import WMAP9 as cosmo
from astropy.io import ascii, fits
import scipy.stats as ss
#from numba import jit
import copy
import emcee
import time
from scipy.integrate import quad,dblquad
#import matplotlib as mpl
#mpl.use('agg')
#import pylab as plt
#import corner
import multiprocessing
from multiprocessing import Pool
import glob
import math
import gzip
import contextlib
import pylab as plt
from chainconsumer import ChainConsumer
import numpy as np
#import matplotlib.pyplot as plt
from astropy.io import ascii
#from sklearn.decomposition import PCA
import time
from scipy import stats
from astropy.cosmology import WMAP9 as cosmo
import pandas as pd
from scipy.integrate import quad
#import twodkstest_peacock as ks2d
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline, RegularGridInterpolator
from scipy.integrate import quad
from scipy.stats import anderson_ksamp,chisquare

'''
to run in grace:
. $HEADAS/headas-init.sh
module load miniconda
conda activate run_marx 
sbatch -p day -t 23:00:00 --mail-type=ALL --cpus-per-task 40 mle_new_for_yale.sh 
sbatch -p week -t 160:00:00 --mail-type=ALL --cpus-per-task 36 mle_new_for_yale.sh 
'''

#define a results dir in your google drive
tt=time.time()
global results_dir
results_dir=""

mpc_to_cm=3.086e+24
w_AD = 0.5
rx_grid_nums=600
non_detected_chandra_samp_size=91
histbins=20

t0=time.time()
#defining MCMC parameters:
#define variables and parameters, for example, priors:
nwalkers,ndim,nsteps,nburn=36,3,1050,10 #oa,r,fscat,gam,frac_24_25

cpu_count=nwalkers
#prior_info=np.array([[1.5, 2.3],[0.1,0.6],[0.1,0.7]])
prior_info="_"+str(ndim)+"_params_"#+str(fctk_fixed)+"_"
#param_bounds=np.array([[0.1,0.6],[0.1,0.7]])
if ndim==6:
	latex_param=["OA","R","fscat",r"$\Gamma$",r"$f_{ctk}$",r"$f_{24-25}$"]#
	param_bounds=np.array([[30,70],[0,2.],[-2,0.5],[1.8,2.0],[0.2,0.7],[0.05,0.95]])

	param_mean=[-99,0.99,-99.,1.9,-99.,-99.]
	sigma=[-99,0.185,-99.,0.085,-99.,-99.]

	priors=[ss.uniform(param_bounds[0][0],param_bounds[0][1]-param_bounds[0][0]),
	        ss.norm(param_mean[1],sigma[1]),
	        ss.uniform(param_bounds[2][0],param_bounds[2][1]-param_bounds[2][0]),
	        ss.norm(param_mean[3],sigma[3]),
	        ss.uniform(param_bounds[4][0],param_bounds[4][1]-param_bounds[4][0]),
	        ss.uniform(param_bounds[5][0],param_bounds[5][1]-param_bounds[5][0])]

	#log prior:
	def log_prior(params):
	        prob=0
	        #flat prior
	        #print("G1",params)
	        if param_bounds[0][0] <= params[0] <= param_bounds[0][1]\
	           and param_bounds[1][0] <= params[1] <= param_bounds[1][1]\
	           and param_bounds[2][0] <= params[2] <= param_bounds[2][1]\
	           and param_bounds[3][0] <= params[3] <= param_bounds[3][1]\
	           and param_bounds[4][0] <= params[4] <= param_bounds[4][1]\
	           and param_bounds[5][0] <= params[5] <= param_bounds[5][1]:
	             for i in np.arange(ndim):
	               prob+=np.log(priors[i].pdf(params[i]))
	               return prob#+np.log(len(obj_info))#0#prob
	        else:
	             return -np.inf

	def log_likelihood(params):

		oa,r,fscat,gam,fctk_fixed,frac_24_25=params
		oa,r,fscat,gam=np.full(model_size,oa),np.full(model_size,r),np.full(model_size,fscat),np.full(model_size,gam)
		#produce log nH distribution in bins of size 0.01
		all_y_bins,x=all_y_bins_and_x_given_fctk_variable(fctk_fixed,frac_24_25)
		
		#model_fscat_sig_vals=np.random.uniform(-2,0.5,model_size)
		model_nH_vals=np.random.choice(all_x_bins,model_size,p=all_y_bins/np.sum(all_y_bins))
		#model_nH_vals=np.clip(model_nH_vals,None,25.)

		put_in_one_bin=np.array([oa,r,fscat,gam,model_nH_vals]).T
		#print("put_in_one_bin",put_in_one_bin)
		rx_vals_model=rx(put_in_one_bin)
		#print("rlx",rx_vals_model)
		#mrx:
		#make a histogram:
		prob=plt.hist(rx_vals_model,bins=int(model_size/10))

		#turn the rlx distribution into a 1d interpolation for integration
		Mrx=interp1d((prob[1][:-1]+prob[1][1:])/2.,prob[0],bounds_error=False,fill_value="extrapolate")

		#non-normalized Mrx values
		mrlx_vals=Mrx(grid_of_rlx)
		mrlx_vals=np.clip(mrlx_vals,0,None)

		#normalize
		norm=np.trapz(mrlx_vals,grid_of_rlx)
		mrlx_vals=mrlx_vals/norm
		#print("mrlx_vals",mrlx_vals,params)

		#normal distributions in all_rlx_grid
		#first term
		first_term_integs=mrlx_vals*first_term_rlx_grid
		all_integs=x_diff*(0.5*(first_term_integs[:,0]+first_term_integs[:,-1])+np.sum(first_term_integs[:,1:-1],axis=1))
		#print("all_integs",all_integs)
		#full first term:
		first_term=np.sum(np.log(all_integs))

		#rlx limit comparison
		#generate 540 Rlx for these parameters? No.
		#add up Mrlx for detected sources
		#second term - should it be a sum
		second_term_integs=mrlx_vals*second_term_rlx_grid
		all_integs=x_diff*(0.5*(second_term_integs[0]+second_term_integs[-1])+np.sum(second_term_integs[1:-1]))
		#mrlx_data=np.clip(mrlx_data,0,None)
		second_term=np.sum(all_integs) #sum of Mrx of all the detected sources as the sensitivity function is a heaviside step function

		print("first and second",first_term,-second_term)
		return first_term - second_term
elif ndim==5:
	latex_param=["OA","R","fscat",r"f$_{ctk}$",r"f$_{24-25}$"]#
	param_bounds=np.array([[30,70],[0,2.],[-2,0.5],[0.2,0.6],[0.05,0.95]])

	param_mean=[-99,0.99,-99.,-99.,-99.]
	sigma=[-99,0.185,-99.,-99.,-99.]

	priors=[ss.uniform(param_bounds[0][0],param_bounds[0][1]-param_bounds[0][0]),
	        ss.norm(param_mean[1],sigma[1]),
	        ss.uniform(param_bounds[2][0],param_bounds[2][1]-param_bounds[2][0]),
	        ss.uniform(param_bounds[3][0],param_bounds[3][1]-param_bounds[3][0]),
	        ss.uniform(param_bounds[4][0],param_bounds[4][1]-param_bounds[4][0])]

	#log prior:
	def log_prior(params):
	        prob=0
	        #flat prior
	        #print("G1",params)
	        if param_bounds[0][0] <= params[0] <= param_bounds[0][1]\
	           and param_bounds[1][0] <= params[1] <= param_bounds[1][1]\
	           and param_bounds[2][0] <= params[2] <= param_bounds[2][1]\
	           and param_bounds[3][0] <= params[3] <= param_bounds[3][1]\
	           and param_bounds[4][0] <= params[4] <= param_bounds[4][1]:
	             for i in np.arange(ndim):
	               prob+=np.log(priors[i].pdf(params[i]))
	               return prob#+np.log(len(obj_info))#0#prob
	        else:
	             return -np.inf

	def log_likelihood(params):

		oa,r,fscat,fctk_fixed,frac_24_25=params
		gam = 1.9
		oa,r,fscat,gam=np.full(model_size,oa),np.full(model_size,r),np.full(model_size,fscat),np.full(model_size,gam)
		#produce log nH distribution in bins of size 0.01
		all_y_bins,x=all_y_bins_and_x_given_fctk_variable(fctk_fixed,frac_24_25)
		
		#model_fscat_sig_vals=np.random.uniform(-2,0.5,model_size)
		model_nH_vals=np.random.choice(all_x_bins,model_size,p=all_y_bins/np.sum(all_y_bins))
		#model_nH_vals=np.clip(model_nH_vals,None,25.)

		put_in_one_bin=np.array([oa,r,fscat,gam,model_nH_vals]).T
		#print("put_in_one_bin",put_in_one_bin)
		rx_vals_model=rx(put_in_one_bin)
		#print("rlx",rx_vals_model)
		#mrx:
		#make a histogram:
		prob=plt.hist(rx_vals_model,bins=int(model_size/10))

		#turn the rlx distribution into a 1d interpolation for integration
		Mrx=interp1d((prob[1][:-1]+prob[1][1:])/2.,prob[0],bounds_error=False,fill_value="extrapolate")

		#non-normalized Mrx values
		mrlx_vals=Mrx(grid_of_rlx)
		mrlx_vals=np.clip(mrlx_vals,0,None)

		#normalize
		norm=np.trapz(mrlx_vals,grid_of_rlx)
		mrlx_vals=mrlx_vals/norm
		#print("mrlx_vals",mrlx_vals,params)

		#normal distributions in all_rlx_grid
		#first term
		first_term_integs=mrlx_vals*first_term_rlx_grid
		all_integs=x_diff*(0.5*(first_term_integs[:,0]+first_term_integs[:,-1])+np.sum(first_term_integs[:,1:-1],axis=1))
		#print("all_integs",all_integs)
		#full first term:
		first_term=np.sum(np.log(all_integs))

		#rlx limit comparison
		#generate 540 Rlx for these parameters? No.
		#add up Mrlx for detected sources
		#second term - should it be a sum
		second_term_integs=mrlx_vals*second_term_rlx_grid
		all_integs=x_diff*(0.5*(second_term_integs[0]+second_term_integs[-1])+np.sum(second_term_integs[1:-1]))
		#mrlx_data=np.clip(mrlx_data,0,None)
		second_term=np.sum(all_integs) #sum of Mrx of all the detected sources as the sensitivity function is a heaviside step function

		print("first and second",first_term,-second_term)
		return first_term - second_term

elif ndim==4:
	latex_param=["OA","fscat",r"f$_{ctk}$",r"f$_{24-25}$"]#
	param_bounds=np.array([[30,65],[-2,0.5],[0.2,0.6],[0.1,0.7]])

	param_mean=[-99,-99.,-99.,-99.]
	sigma=[-99,-99.,-99.,-99.]

	priors=[ss.uniform(param_bounds[0][0],param_bounds[0][1]-param_bounds[0][0]),
	        ss.uniform(param_bounds[1][0],param_bounds[1][1]-param_bounds[1][0]),
	        ss.uniform(param_bounds[2][0],param_bounds[2][1]-param_bounds[2][0]),
	        ss.uniform(param_bounds[3][0],param_bounds[3][1]-param_bounds[3][0])]

	#log prior:
	def log_prior(params):
	        prob=0
	        #flat prior
	        #print("G1",params)
	        if param_bounds[0][0] <= params[0] <= param_bounds[0][1]\
	           and param_bounds[1][0] <= params[1] <= param_bounds[1][1]\
	           and param_bounds[2][0] <= params[2] <= param_bounds[2][1]\
	           and param_bounds[3][0] <= params[3] <= param_bounds[3][1]:
	             for i in np.arange(ndim):
	               prob+=np.log(priors[i].pdf(params[i]))
	               return prob#+np.log(len(obj_info))#0#prob
	        else:
	             return -np.inf

	def log_likelihood(params):

		oa,fscat,fctk_fixed, frac_24_25=params
		r,gam = 0.99,1.9
		oa,r,fscat,gam=np.full(model_size,oa),np.full(model_size,r),np.full(model_size,fscat),np.full(model_size,gam)
		#produce log nH distribution in bins of size 0.01
		all_y_bins,x=all_y_bins_and_x_given_fctk_variable(fctk_fixed,frac_24_25)
		
		#model_fscat_sig_vals=np.random.uniform(-2,0.5,model_size)
		model_nH_vals=np.random.choice(all_x_bins,model_size,p=all_y_bins/np.sum(all_y_bins))
		#model_nH_vals=np.clip(model_nH_vals,None,25.)

		put_in_one_bin=np.array([oa,r,fscat,gam,model_nH_vals]).T
		#print("put_in_one_bin",put_in_one_bin)
		rx_vals_model=rx(put_in_one_bin)
		#print("rlx",rx_vals_model)
		#mrx:
		#make a histogram:
		prob=plt.hist(rx_vals_model,bins=int(model_size/10))

		#turn the rlx distribution into a 1d interpolation for integration
		Mrx=interp1d((prob[1][:-1]+prob[1][1:])/2.,prob[0],bounds_error=False,fill_value="extrapolate")

		#non-normalized Mrx values
		mrlx_vals=Mrx(grid_of_rlx)
		mrlx_vals=np.clip(mrlx_vals,0,None)

		#normalize
		norm=np.trapz(mrlx_vals,grid_of_rlx)
		mrlx_vals=mrlx_vals/norm
		#print("mrlx_vals",mrlx_vals,params)

		#normal distributions in all_rlx_grid
		#first term
		first_term_integs=mrlx_vals*first_term_rlx_grid
		all_integs=x_diff*(0.5*(first_term_integs[:,0]+first_term_integs[:,-1])+np.sum(first_term_integs[:,1:-1],axis=1))
		#print("all_integs",all_integs)
		#full first term:
		first_term=np.sum(np.log(all_integs))

		#rlx limit comparison
		#generate 540 Rlx for these parameters? No.
		#add up Mrlx for detected sources
		#second term - should it be a sum
		second_term_integs=mrlx_vals*second_term_rlx_grid
		all_integs=x_diff*(0.5*(second_term_integs[0]+second_term_integs[-1])+np.sum(second_term_integs[1:-1]))
		#mrlx_data=np.clip(mrlx_data,0,None)
		second_term=np.sum(all_integs) #sum of Mrx of all the detected sources as the sensitivity function is a heaviside step function

		#print("first and second",first_term,-second_term)
		return first_term - second_term
elif ndim==3:
	latex_param=["fscat",r"f$_{ctk}$",r"f$_{24-25}$"]#
	param_bounds=np.array([[-2,0.5],[0.2,0.7],[0.1,0.8]])

	param_mean=[-99.,-99.,-99.]
	sigma=[-99.,-99.,-99.]

	priors=[ss.uniform(param_bounds[0][0],param_bounds[0][1]-param_bounds[0][0]),
			ss.uniform(param_bounds[1][0],param_bounds[1][1]-param_bounds[1][0]),
			ss.uniform(param_bounds[2][0],param_bounds[2][1]-param_bounds[2][0])]

	#log prior:
	def log_prior(params):
	        prob=0
	        #flat prior
	        #print("G1",params)
	        if param_bounds[0][0] <= params[0] <= param_bounds[0][1]\
	           and param_bounds[1][0] <= params[1] <= param_bounds[1][1]\
	           and param_bounds[2][0] <= params[2] <= param_bounds[2][1]:
	             for i in np.arange(ndim):
	               prob+=np.log(priors[i].pdf(params[i]))
	               return prob#+np.log(len(obj_info))#0#prob
	        else:
	             return -np.inf

	def log_likelihood(params):

		fscat,fctk_fixed, frac_24_25=params
		oa,r,gam = 60.,0.99,1.9
		oa,r,fscat,gam=np.full(model_size,oa),np.full(model_size,r),np.full(model_size,fscat),np.full(model_size,gam)
		#produce log nH distribution in bins of size 0.01
		all_y_bins,x=all_y_bins_and_x_given_fctk_variable(fctk_fixed,frac_24_25)
		
		#model_fscat_sig_vals=np.random.uniform(-2,0.5,model_size)
		model_nH_vals=np.random.choice(all_x_bins,model_size,p=all_y_bins/np.sum(all_y_bins))
		#model_nH_vals=np.clip(model_nH_vals,None,25.)

		put_in_one_bin=np.array([oa,r,fscat,gam,model_nH_vals]).T
		#print("put_in_one_bin",put_in_one_bin)
		rx_vals_model=rx(put_in_one_bin)
		#print("rlx",rx_vals_model)
		#mrx:
		#make a histogram:
		prob=plt.hist(rx_vals_model,bins=int(model_size/10))

		#turn the rlx distribution into a 1d interpolation for integration
		Mrx=interp1d((prob[1][:-1]+prob[1][1:])/2.,prob[0],bounds_error=False,fill_value="extrapolate")

		#non-normalized Mrx values
		mrlx_vals=Mrx(grid_of_rlx)
		mrlx_vals=np.clip(mrlx_vals,0,None)

		#normalize
		norm=np.trapz(mrlx_vals,grid_of_rlx)
		mrlx_vals=mrlx_vals/norm
		#print("mrlx_vals",mrlx_vals,params)

		#normal distributions in all_rlx_grid
		#first term
		first_term_integs=mrlx_vals*first_term_rlx_grid
		all_integs=x_diff*(0.5*(first_term_integs[:,0]+first_term_integs[:,-1])+np.sum(first_term_integs[:,1:-1],axis=1))
		#print("all_integs",all_integs)
		#full first term:
		first_term=np.sum(np.log(all_integs))

		#rlx limit comparison
		#generate 540 Rlx for these parameters? No.
		#add up Mrlx for detected sources
		#second term - should it be a sum
		second_term_integs=mrlx_vals*second_term_rlx_grid
		all_integs=x_diff*(0.5*(second_term_integs[0]+second_term_integs[-1])+np.sum(second_term_integs[1:-1]))
		#mrlx_data=np.clip(mrlx_data,0,None)
		second_term=np.sum(all_integs) #sum of Mrx of all the detected sources as the sensitivity function is a heaviside step function

		#print("first and second",first_term,-second_term)
		return first_term - second_term


#importing specific files
#ricci:
lans=ascii.read(results_dir+"Ricci+2017_Fig23_data.csv")
#lansbury:
#lans=ascii.read(results_dir+"lansbury_2015_nh_distribution.txt")

lans_x=np.array([21.0,21.999,22.0,22.999,23.0,24.])
#np.concatenate((np.array([21.]),lans["NH"].data[2:]))
lans_y=np.array([lans["FREQ"].data[2],lans["FREQ"].data[2],lans["FREQ"].data[3],
								 lans["FREQ"].data[3],lans["FREQ"].data[4],lans["FREQ"].data[4]])
#np.array([2.92,2.92,6.57,6.57,15.81,15.81])

lans_interp=interp1d(lans_x,lans_y)
area_current=quad(lans_interp,21,24)[0]
correct_dist=np.arange(21.,24,.01)
ctk_bins=np.arange(24.,26.005,0.01)
first_norm_lans_dist=lans_interp(correct_dist)/area_current
all_x_bins=np.concatenate((correct_dist,ctk_bins))
all_x_bins=np.clip(all_x_bins,None,26.0)
model_size=2000
sigma_for_uncertainty=0.23

def all_y_bins_and_x_given_fctk(fctk):
	x = fctk/(1.-fctk)
	each_bin=x/2.
	all_y_bins=np.concatenate((first_norm_lans_dist,np.full(len(ctk_bins),each_bin)))
	all_y_bins=all_y_bins/(1.+x)
	return all_y_bins,x


def all_y_bins_and_x_given_fctk_variable(fctk,frac_24_26):
	x = fctk/(1.-fctk)
	each_bin=frac_24_26*x
	#append probability for the log nH 24-25 bin
	all_y_bins=np.concatenate((first_norm_lans_dist,
							   np.full(int(len(ctk_bins)/2),each_bin)))
	each_bin=(1.-frac_24_26)*x
	#append probability for the log nH 24-25 bin
	all_y_bins=np.concatenate((all_y_bins,
							   np.full(int(len(ctk_bins)/2)+1,each_bin)))
	all_y_bins=all_y_bins/(1.+x)
	return all_y_bins,x


#now associate model to observations:
new_rx=fits.open("rx.fits")
#new_rx_dat=new_rx[1].data
#new_rx_dat=Table(new_rx_dat)
rx=RegularGridInterpolator((new_rx[7].data,new_rx[6].data,new_rx[5].data,
                                 new_rx[4].data, new_rx[3].data), new_rx[0].data, method='linear',
                                bounds_error=False,fill_value=None)

hardx_=RegularGridInterpolator((new_rx[8].data,new_rx[7].data,new_rx[6].data,new_rx[5].data,
                                 new_rx[4].data, new_rx[3].data), new_rx[2].data, method='linear',
                                bounds_error=True)#,fill_value=0)


softx_=RegularGridInterpolator((new_rx[8].data,new_rx[7].data,new_rx[6].data,new_rx[5].data,
                                 new_rx[4].data, new_rx[3].data), new_rx[1].data, method='linear',
                                bounds_error=True)

rxl=ascii.read(results_dir+"rx_data.csv")
rx_data=rxl["RX"].data
rx_data_lim=rxl["RXL"].data
#detected sources:
detected_data=np.where(np.greater_equal(rx_data,rxl["RXL"].data))[0]
detected_rx_values=rx_data[detected_data]

grid_of_rlx=np.linspace(-3,1,rx_grid_nums)
first_term_rlx_grid=[]
for mu,mu_lim in zip(rx_data,rx_data_lim):
  if mu > -4.2:
    blah=ss.norm.pdf(grid_of_rlx,mu,0.23)
    norm=np.trapz(blah,grid_of_rlx)
    first_term_rlx_grid.append(blah/norm)


first_term_rlx_grid=np.array(first_term_rlx_grid)
#first_term_rlx_grid=np.sum(first_term_rlx_grid,axis=0)
#print("first term",first_term_rlx_grid)

second_term_rlx_grid=[]
for mu_lim in rx_data_lim:
	inds=np.where(grid_of_rlx>mu_lim)[0]
	blah=np.zeros(len(grid_of_rlx))
	blah[inds] = 1.
	second_term_rlx_grid.append(blah)

second_term_rlx_grid=np.array(second_term_rlx_grid)
second_term_rlx_grid=np.sum(second_term_rlx_grid,axis=0)
#print("second term",second_term_rlx_grid)

x_diff=grid_of_rlx[1]-grid_of_rlx[0]
#print("all_rlx_grid",all_rlx_grid)
#add_uncertainty

def lnprob(params):
	#add log prior and log likelihood
	tt1=time.time()
	#print("G0 starting prior",params)
	lp=log_prior(params)
	if not np.isfinite(lp):
		return -np.inf
	else:
		ll=log_likelihood(params)
		total=lp+ll
		gc.collect()
		#if not np.isinf(ll):
		#	print("PARAM OUTPUT:", params, ll, lp, total, 'Time:',time.time()-tt1)
		if np.isnan(total):
			return -np.inf
		else:
			return total


with contextlib.closing(Pool(processes=nwalkers)) as pool:
		pos=np.full(nwalkers*ndim, 1.25).reshape(nwalkers,ndim)
		for i in range(0,ndim):
			pos[:,i]=priors[i].rvs(nwalkers)

			#print("G")
		'''
		filename = results_dir+"tutorial_fctk_"+str(fctk_fixed)+".h5"
		if glob.glob(filename):
			os.remove(filename)
		backend = emcee.backends.HDFBackend(filename)
		backend.reset(nwalkers, ndim)
		
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, backend=backend)#, threads=nwalkers)#, pool=pool)
		'''
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)#, threads=nwalkers)#, pool=pool)

		#print("H")

		#burn in:
		pos = sampler.run_mcmc(pos, nburn)#[0]
		print("done with burns in", time.time()-t0)
		#reset
		sampler.reset()

		#sampling:
		result = sampler.run_mcmc(pos, nsteps, progress=True)#[0]
		pool.close()
		np.save(results_dir+"chris_mcmc_chain_"+prior_info+str(nwalkers)+"_walkers_"+str(nsteps)+"_steps.npy", sampler.chain)

os.system("gzip "+results_dir+"chris_mcmc_chain_"+prior_info+str(nwalkers)+"_walkers_"+str(nsteps)+"_steps.npy")
print("total time for "+str(int(nsteps+nburn))+" runs:", time.time()-t0)
#print_all_results(sampler.chain, False)

sampler_chain=sampler.chain
figu, ax = plt.subplots(ndim,figsize=(8,4*ndim))

shapeyy=sampler_chain.shape[1]

#for i in range(0,ndim):
for i in range(0,ndim):
    ax[i].plot(sampler_chain[:,:,i].T, '-', color='k', alpha=0.3)
    #ax[i].set_title("weight # "+str(i+1))
    ax[i].set_xlabel(latex_param[i])
plt.subplots_adjust(hspace=0.2)
plt.savefig("chris_mcmc_chain_"+prior_info+str(nwalkers)+"_walkers_"+str(nsteps)+"_steps.png", dpi=100)
plt.close()

sampler_chain=sampler_chain.reshape(-1,ndim)

c = ChainConsumer().add_chain(sampler_chain, parameters=latex_param, 
		name=prior_info)

table = c.analysis.get_latex_table(caption="Summary", label="tab:summary_mcmc")

print(table)

cov_table=c.analysis.get_covariance_table(caption='Parameter Covariance', 
										  label='tab:_parameter_covariance')

print("Covariance table:")
print(cov_table)

c.configure(colors=["blue"], shade=True, shade_alpha=0.7, bar_shade=True,#"#D50B53""#1DCFD6"#BD3E85
			flip=False, sigma2d=True, summary=True, diagonal_tick_labels=False,
			tick_font_size=7, sigmas=[1,2], max_ticks=5,label_font_size=12)# 

fig = c.plotter.plot(filename="chain_consumer_chris_mcmc_chain_"+prior_info+str(nwalkers)+"_walkers_"+str(nsteps)+"_steps.png")

print("Total time:",time.time()-tt)