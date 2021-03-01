#!/cvmfs/icecube.opensciencegrid.org/py3-v4/RHEL_7_x86_64/bin/python

import sys, os, time

#CSKY IMPORTS
import csky as cy
from csky.hyp import PowerLawFlux
from csky.utils import Arrays

#OTHER IMPORTS
import copy
import socket
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import histlite as hl
import datetime
from astropy.time import Time
import pickle  
import data_rebinning
import argparse


################################################################################
# DEFINITIONS
################################################################################

def Write2File(msg,logfile):

    if logfile!="":
        log = open(logfile, "a")
        log.write(msg+"\n")
        log.close

    return
#################################################################################


################################################################################
# INPUT SECTION
################################################################################

p.add_argument("--pivot_energy", default=100.0, type=float,
               help="Pivot energy in TeV (E0) (default=100.0)")
p.add_argument("--wrkdir", default="/data/user/adesai/csky/mojave/mojave_time_integrated_analysis_gthub_rep/n_bias_test", type=str,
               help="Output directory (default:/data/user/adesai/csky/mojave/mojave_time_integrated_analysis_gthub_rep/n_bias_test)")
p.add_argument("--index", default=2.0, type=float,
               help="Spectral Index (default=2.0)")
p.add_argument("--ana_dir_path", default="/data/user/adesai/csky/mojave/time_integrated_analysis/cache/", type=str,
               help="Location of analysis directory. (Default:/data/user/adesai/csky/mojave/time_integrated_analysis/cache/)")




args = p.parse_args()

pivot_energy = args.pivot_energy
index = args.index
wrkdir = args.wrkdir
ana_dir_path = args.ana_dir_path


if(os.path.isdir(wrkdir)!=True):
    print("Input ERROR, directory_to_make_new_folders does not exist")
    sys.exit()



GeV = 1.        # base SkyLab energy units are GeV                        
TeV = 1000*GeV    
E0    = pivot_energy*TeV                # pivot energy 




################################################################################
# GET MOJAVE DATA SET (NAME RA DEC av_flux pk_flux)SET Directory
################################################################################

input_file="/data/user/adesai/radio_stacking/mojavexv/1A_wt_rd_flux/MOJAVE_multiyear_complete_sample.txt"


array = np.loadtxt(input_file,dtype='str')
src_name=np.asarray(array[:,0])
src_ra=np.asarray(array[:,1],dtype=float)*np.pi/180
src_dec=np.asarray(array[:,2],dtype=float)*np.pi/180
flux_used=np.asarray(array[:,3],dtype=float)


########################################
#  Sample Checks and outdir
########################################

print("---------------------------------------------------------\n")
print("WARNING!!!! CHECK BELOW IF AVG FLUX IS USED (STATEMENT IS TRUE)")
print("Max Avg Flux of sample %s; Flux used:%s"%(21433.22727,np.amax(flux_used)))
if np.amax(flux_used)!=21433.22727:
      print("ERROR: Fluxes not same, Check column number used!!!!")
      #exit() #IN IPYNB THIS WILL KILL THE KERNEL
print("---------------------------------------------------------\n\n\n")


#if len(flux_used)!=no_of_sources:                                                                       
print("---------------------------------------------------------\n")
print("Number of sources used in stacking: %s"%(len(flux_used)))
print("---------------------------------------------------------\n")
no_of_sources=len(flux_used)






    
################################################################################
# Weights!
################################################################################

equal_weights=False

wt=[]
for j in range(len(flux_used)):
    if equal_weights:
        wt.append(1.0)
    else:
        wt.append(flux_used[j])

wt=wt/np.sum(wt)
print("\n WEIGHTS USED FOR THIS RUN:\n %s \n"%(wt))


################################################################################
# CSKY SETUP
################################################################################

ana_dir = cy.utils.ensure_dir(ana_dir_path)
repo = cy.selections.repo

ana = cy.get_analysis(cy.selections.repo, cy.selections.PSDataSpecs.ps_10yr, dir=ana_dir)
#ana.save(ana_dir) #Un-comment if you are running the analysis for the first time

cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = 2





################################################################################
# CSKY SOURCE DESCRIPTION
################################################################################

src = cy.utils.Sources(ra=src_ra,dec=src_dec,weight=wt,deg=False)

flux    = cy.hyp.PowerLawFlux(gamma=index,norm=E0)




################################################################################
# BACKGROUND
################################################################################

#conf = dict(fitter_args = _fmin_method='bfgs')

timer = cy.timing.Timer()
time = timer.time
with time('trial runner construction'):
    tr = cy.get_trial_runner(src=src, ana=ana,flux=flux, mp_cpus=2, sindec_bandwidth=np.radians(.1))

          
                              

################################################################################
# N Bias test
################################################################################

n_sig_arr =[1,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,437]


#line="n_inj\tn_sig\ttrial_arr"
#Write2File(line ,wrkdir+"results.txt")

trials = [tr.get_many_fits(100, n_sig=n_sig, logging=False, seed=n_sig*np.random.randint(1,1000)) for n_sig in n_sig_arr]

for (n_sig, t) in zip(n_sig_arr, trials):
    t['ntrue'] = np.repeat(n_sig, len(t))

time = timer.time

allt = cy.utils.Arrays.concatenate(trials)
#allt
np.savez(wrkdir+"trials_with_sig_inj_"+str(index)+".npz",allt.as_array)



fig, axs = plt.subplots(1, 2, figsize=(6,3))
dns = np.mean(np.diff(n_sig_arr))
ns_bins = np.r_[n_sig_arr - 0.5*dns, n_sig_arr[-1] + 0.5*dns]
expect_kw = dict(color='C0', ls='--', lw=1, zorder=-10)
expect_gamma = tr.sig_injs[0].flux[0].gamma

ax = axs[0]
h = hl.hist((allt.ntrue, allt.ns), bins=(ns_bins, 100))
hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')

lim = ns_bins[[0, -1]]
ax.set_xlim(ax.set_ylim(lim))
ax.plot(lim, lim, **expect_kw)
ax.set_aspect('equal')


ax = axs[1]
h = hl.hist((allt.ntrue, allt.gamma), bins=(ns_bins, 100))
hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
ax.axhline(expect_gamma, **expect_kw)
ax.set_xlim(axs[0].get_xlim())

for ax in axs:
    ax.set_xlabel(r'$n_{inj}$')
    ax.grid()
axs[0].set_ylabel(r'$n_s$')
axs[1].set_ylabel(r'$\gamma$')

plt.tight_layout()
plt.savefig(wrkdir+"trials_with_sig_inj_"+str(index)+".png")   
    
