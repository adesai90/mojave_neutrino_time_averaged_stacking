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
import random

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
    wrkdir      = "/data/user/adesai/csky/mojave/time_integrated_analysis/results_with_fixed_energy_inj/"   
    #energy_inj    = command_line[2] #TeV
    index = float(command_line[2])      # spectral index   


# Check if input directory exists        
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




outdir="%s/avg_flux_weight_stacking_%s_sources_gamma_%s_pvenergy_%s/"%(wrkdir,len(flux_used),index,pivot_energy)

if(os.path.isdir(outdir)!=True):
    os.system('mkdir %s' %outdir)
os.chdir(outdir)
print("Current WORKING DIRECTORY is: %s"%os.getcwd())


"""

tr_file="trialrunner_sensitivity_results_for_index_%s_energy_bin_%s.pickle"%(index,energy_inj)
results_file="sensitivity_results_for_index_%s_energy_bin_%s.npy"%(index,energy_inj)

if(os.path.isfile(results_file)==True):
    print("\n \n WARNING: \n RESULT TAKEN FROM PRECOMPUTED ANALYSIS")
    filehandler = open(tr_file, 'rb') 
    tr=pickle.load(filehandler)
    filehandler.close
    sens_res=np.load(results_file,allow_pickle=True)[()]
    print("\nSensitivity Flux in TeV/cm2/s  @ %s TeV:"%(E0/TeV))
    print(tr.to_E2dNdE(sens_res, E0=E0/TeV, unit=1e3))   # TeV/cm2/s  @  100TeV
else:
    print("No precomputed result so running analysis")
    
"""    
    
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

flux    = cy.hyp.PowerLawFlux(gamma=index,norm=E0)
#,energy_range=(inj_enj_lo,inj_enj_hi))



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

bkgfile_comb= "/data/user/adesai/csky/mojave/time_integrated_analysis/bkg_files/nondefault_mostsenseng_background_for_index_%s.npy"%(index)


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
        n_sig_step=2,
        batch_size=500,
        max_batch_size=500,
        # 10% tolerance -- let's get an estimate quick!
        tol=0.05
        )
    #disc = tr.find_n_sig(bg.isf_nsigma(5), 0.5, n_sig_step=5, batch_size=100, tol=0.05)



print("\nSensitivity Flux in TeV/cm2/s  @ %s TeV:"%(E0/TeV))

ednde=tr.to_dNdE(sens['n_sig'], E0=pivot_energy, unit=1e3)
ednde_err=tr.to_dNdE(sens['n_sig_error'], E0=pivot_energy, unit=1e3)
e2dnde=tr.to_E2dNdE(sens['n_sig'], E0=pivot_energy, unit=1e3)
e2dnde_err=tr.to_E2dNdE(sens['n_sig_error'], E0=pivot_energy, unit=1e3)


print(e2dnde,"+-",e2dnde_err)
print(timer)

################################################################################
# SAVE RESULTS TO COMMON FILE
################################################################################

combresult_file= "/data/user/adesai/csky/mojave/time_integrated_analysis/results_with_fixed_energy_inj/nondefault_mostsenseng_sensitivity_fixed_einj_result_file_index_%s.txt"%(index)
if(os.path.isfile(combresult_file)!=True):
    Write2File("#Gamma\tE_norm\tN_sig\tdnde\tdnde_err\tsensflux(E2dnde)calc\tsensflux(E2dnde)err",combresult_file)
    
pv_eng=pivot_energy

ednde=tr.to_dNdE(sens['n_sig'], E0=pv_eng, unit=1e3)
ednde_err=tr.to_dNdE(sens['n_sig_error'], E0=pv_eng, unit=1e3)
e2dnde=tr.to_E2dNdE(sens['n_sig'], E0=pv_eng, unit=1e3)
e2dnde_err=tr.to_E2dNdE(sens['n_sig_error'], E0=pv_eng, unit=1e3)

out_line="%s\t%s\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e"%(index,pv_eng,sens['n_sig'],ednde,ednde_err,e2dnde,e2dnde_err)
    
Write2File(out_line,combresult_file)


################################################################################                            
# DISCOVERY POTENTIAL                                                                                   
################################################################################                                    
disc = tr.find_n_sig(bg.isf_nsigma(5), 0.5, n_sig_step=2, batch_size=500, tol=0.05)

combresult_file2= "/data/user/adesai/csky/mojave/time_integrated_analysis/results_with_fixed_energy_inj/nondefault_mostsenseng_discovery_fixed_einj_result_file_%s.txt"%(index)
if(os.path.isfile(combresult_file2)!=True):
    Write2File("#Gamma\tE_norm\tN_sig\tdnde\tdnde_err\tdiscflux(E2dnde)calc\discflux(E2dnde)err",combresult_file2)

   
    
ednde=tr.to_dNdE(disc['n_sig'], E0=pv_eng, unit=1e3)
ednde_err=tr.to_dNdE(disc['n_sig_error'], E0=pv_eng, unit=1e3)
e2dnde=tr.to_E2dNdE(disc['n_sig'], E0=pv_eng, unit=1e3)
e2dnde_err=tr.to_E2dNdE(disc['n_sig_error'], E0=pv_eng, unit=1e3)

print("\nDiscovery Potential in TeV/cm2/s  @ %s TeV:"%(E0/TeV))
print(e2dnde,"+-",e2dnde_err)


out_line="%s\t%s\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e"%(index,pv_eng,disc['n_sig'],ednde,ednde_err,e2dnde,e2dnde_err)
    
Write2File(out_line,combresult_file2)                                                                                                                                                 
                                                                                                                                                 
                                                                                                                                                 
