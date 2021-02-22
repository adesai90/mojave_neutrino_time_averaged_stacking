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
import random
import pickle 


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

command_line = sys.argv
if len(command_line) < 3:
    print("WRONG NUMBER OF INPUTS, CHECK CODE")
    sys.exit(1)
else:
    pivot_energy   = float(command_line[1])
    wrkdir      = "/data/user/adesai/csky/mojave/time_integrated_analysis/sensitivity_with_diff_inj_enj/"   
    energy_inj    = command_line[2] #TeV
    index = float(command_line[3])      # spectral index   

# Check if input directory exists        
if(os.path.isdir(wrkdir)!=True):
    print("Input ERROR, directory_to_make_new_folders does not exist")
    sys.exit()



GeV = 1.        # base SkyLab energy units are GeV                        
TeV = 1000*GeV    

E0    = pivot_energy*TeV                # pivot energy 
inj_enj_lo=(float(energy_inj[:3]))*TeV   # inj energy LL
inj_enj_hi=(float(energy_inj[4:]))*TeV   # inj energy UL 

print("Min max inj energy")
print(inj_enj_lo,inj_enj_hi)



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




outdir="%s/avg_flux_weight_stacking_%s_sources_gamma_%s_pvenergy_%s_injenj_%s/"%(wrkdir,len(flux_used),index,E0,pivot_energy)

if(os.path.isdir(outdir)!=True):
    os.system('mkdir %s' %outdir)
os.chdir(outdir)
print("Current WORKING DIRECTORY is: %s"%os.getcwd())


  
    
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
# CSKY SETUP and LOAD ICECUBE DATA
################################################################################



hostname = socket.gethostname ()
username = 'adesai'
ana_dir = cy.utils.ensure_dir('/data/user/adesai/csky/mojave/time_integrated_analysis/cache/')
repo = cy.selections.repo

ana = cy.get_analysis(cy.selections.repo, cy.selections.PSDataSpecs.ps_10yr, dir=ana_dir)
#ana.save(ana_dir)
job_basedir = outdir

cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = 20





################################################################################
# CSKY SOURCE DESCRIPTION
################################################################################

src = cy.utils.Sources(ra=src_ra,dec=src_dec,weight=wt,deg=False)

flux    = cy.hyp.PowerLawFlux(gamma=index,norm=E0,energy_range=(inj_enj_lo,inj_enj_hi))




################################################################################
# SET UP TRIAL RUNNER
################################################################################

timer = cy.timing.Timer()
time = timer.time
with time('trial runner construction'):
    tr = cy.get_trial_runner(src=src, ana=ana,flux=flux, sindec_bandwidth=np.radians(.1), mp_cpus=20)


################################################################################
# BACKGROUND
################################################################################

bkgfile_comb= "/data/user/adesai/csky/mojave/time_integrated_analysis/bkg_files/background_for_index_%s_energy_inj_%s.npy"%(index,energy_inj)


if(os.path.isfile(bkgfile_comb)==True):
    bg_arr=np.load(bkgfile_comb)
    bg = cy.dists.Chi2TSD(cy.utils.Arrays(init=bg_arr))
else:
    print("computing background:")
    bg_arr = tr.get_many_fits(1000, seed=random.randint(1, 1000))
    bg = cy.dists.Chi2TSD(bg_arr)
    np.save(bkgfile_comb,bg_arr.as_array)
    
print(bg.description)
################################################################################
# BACKGROUND PLOTTING (OPTIONAL)
################################################################################
                              
"""
fig, ax = plt.subplots()

h = bg.get_hist(bins=30)
hl.plot1d(ax, h, crosses=True, label='{} bg trials'.format(bg.n_total))

x = h.centers[0]
norm = h.integrate().values
ax.semilogy(x, norm * bg.pdf(x), lw=1, ls='--',label=r'$\chi^2[{:.2f}\ dof,\ \eta={:.3f}]$'.format(bg.ndof, bg.eta))

ax.set_xlabel(r'TS')
ax.set_ylabel(r'number of trials')
ax.legend()
plt.tight_layout()
plt.savefig("background_plot_for_index_%s_energy_bin_%s.pickle"%(index,energy_inj))
"""


################################################################################
# SENSITIVITY
################################################################################

with time('sensitivity'):
    sens = tr.find_n_sig(
        bg.median(), 0.9,
        n_sig_step=5,
        batch_size=100,
        max_batch_size=100,
        tol=0.05
        )
    disc = tr.find_n_sig(bg.isf_nsigma(5), 0.5, n_sig_step=5, batch_size=100, tol=0.05)
    
e_mid=10**((math.log((inj_enj_lo/1000),10)+math.log((inj_enj_hi/1000),10))/2)

print("\nSensitivity Flux in TeV/cm2/s  @ %s TeV:"%(E0/TeV))
print(tr.to_E2dNdE(sens, E0=e_mid, unit=1e3))   # TeV/cm2/s  @  100TeV
e2dnde=tr.to_E2dNdE(sens['n_sig'], E0=e_mid, unit=1e3)
e2dnde_err=tr.to_E2dNdE(sens['n_sig_error'], E0=e_mid, unit=1e3)

#print(tr.to_E2dNdE(disc, E0=e_mid, unit=1e3))   # TeV/cm2/s  @  100TeV
#e2dnde=tr.to_E2dNdE(disc['n_sig'], E0=e_mid, unit=1e3)
#e2dnde_err=tr.to_E2dNdE(disc['n_sig_error'], E0=e_mid, unit=1e3)


print(e2dnde,"+-",e2dnde_err)

################################################################################
# SAVE RESULTS TO COMMON FILE
################################################################################

#E_range gamma e_norm sensflux(e2dnde)calc sensflux(e2dnde)err
out_line="%s\t%s\t%s\t%.2e\t%.2e\t%s"%(energy_inj,index,e_mid,e2dnde,e2dnde_err,"TeV/cm2/s")

combresult_file= "/data/user/adesai/csky/mojave/time_integrated_analysis/sens_combined_result_file_index_2_0.txt"
if(os.path.isfile(combresult_file)!=True):
    Write2File("E_range\tgamma\te_norm\tsensflux(E2dnde)calc\tsensflux(E2dnde)err",combresult_file)

Write2File(out_line,combresult_file)


################################################################################    
# SAVE RESULTS                                                                                
################################################################################    

print("\n DP Flux in TeV/cm2/s  @ %s TeV:"%(E0/TeV))
print(tr.to_E2dNdE(disc, E0=e_mid, unit=1e3))   # TeV/cm2/s  @  100TeV
e2dnde=tr.to_E2dNdE(disc['n_sig'], E0=e_mid, unit=1e3)
e2dnde_err=tr.to_E2dNdE(disc['n_sig_error'], E0=e_mid, unit=1e3)


#E_range gamma e_norm sensflux(e2dnde)calc sensflux(e2dnde)err                                       
out_line="%s\t%s\t%s\t%.2e\t%.2e\t%s"%(energy_inj,index,e_mid,e2dnde,e2dnde_err,"TeV/cm2/s")

combresult_file= "/data/user/adesai/csky/mojave/time_integrated_analysis/dp_combined_result_file_i\
ndex_2_0.txt"
if(os.path.isfile(combresult_file)!=True):
    Write2File("E_range\tgamma\te_norm\tdpflux(E2dnde)calc\tdpflux(E2dnde)err",combresult_file)

Write2File(out_line,combresult_file)

    
