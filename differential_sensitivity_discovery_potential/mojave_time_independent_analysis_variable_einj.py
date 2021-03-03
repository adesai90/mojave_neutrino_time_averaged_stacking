#!/cvmfs/icecube.opensciencegrid.org/py3-v4/RHEL_7_x86_64/bin/python

import os, time
#CSKY IMPORTS
import csky as cy
from csky.hyp import PowerLawFlux
from csky.utils import Arrays
#OTHER IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import histlite as hl
from astropy.time import Time
import random
import argparse

################################################################################
# DEFINITIONS
################################################################################
def Write2File(msg,logfile):
    if logfile!="":
        log = open(logfile, "a")
        log.write(msg+"\n")
        log.close()
    else:
        print("WARNING: File path empty, nothing is written")
    return

def computing_fluxes_from_ns(trials,sens_disc_arr,pivot_energy,units):
    ednde = trials.to_dNdE(sens_disc_arr['n_sig'], E0=pivot_energy, unit=units)
    ednde_err = trials.to_dNdE(sens_disc_arr['n_sig_error'], E0=pivot_energy, unit=units)
    e2dnde = trials.to_E2dNdE(sens_disc_arr['n_sig'], E0=pivot_energy, unit=units)
    e2dnde_err = trials.to_E2dNdE(sens_disc_arr['n_sig_error'], E0=pivot_energy, unit=units)
    return ednde,ednde_err,e2dnde,e2dnde_err
################################################################################
# INPUT SECTION
################################################################################
#############
# ARGUMENTS #
#############
p = argparse.ArgumentParser(description="Calculates Sensitivity and Discovery"
                            " Potential Fluxes for the MOJAVE AGN Radio stacking analysis",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--pivot_energy", default=100.0, type=float,
               help="Pivot energy in TeV (E0) (default=100.0)")
p.add_argument("--index", default=2.0, type=float,
               help="Spectral Index (default=2.0)")
p.add_argument("--wrkdir", default="result_dir", type=str,
               help="Output directory (default:result_dir)")
p.add_argument("--equal_weights", default=False, 
               help="Equal weights for stacking? If no then average fluxes are used. (default:False)")
p.add_argument("--nscramble", default=1000, type=int,
               help="Number of background only scrambles used "
               "to measure TS distribution (default=1000)")
p.add_argument("--nsample", default=100, type=int,
               help="Batch size for sensitivity and discovery potential calculation ")
p.add_argument("--make_background_plot", default=False,
               help="Make plot of background trials (defaut:False)")
p.add_argument("--inj_energy", default="1e0_1e1", type=str, 
               help="Energy bin in TeV to use for computing differential sensitivity in form of string 1e0_1e1 (default)")
p.add_argument("--nstep", default=5, type=int,
               help="Number of signal injection steps (default=5)")
p.add_argument("--input_file", default="../MOJAVE_multiyear_complete_sample.txt", type=str,
               help="Location of input file: MOJAVE_multiyear_complete_sample.txt")
p.add_argument("--ana_dir_path", default="../cache/", type=str,
               help="Location of analysis directory. (Default:../cache/)")
p.add_argument("--discovery_thresh", default=5, type=float,
               help="Discovery threshold in sigma (default=5)")
p.add_argument("--cpus_used",default = 2, type=int,help = 'number of cpus to be used (default 2)')

args = p.parse_args()
#############
# ARGUMENTS #
#############    
pivot_energy = args.pivot_energy
index = args.index
wrkdir = args.wrkdir
equal_weights = args.equal_weights
nscramble = args.nscramble
nsample = args.nsample
make_background_plot =args.make_background_plot
inj_energy = args.inj_energy
nstep = args.nstep
input_file = args.input_file
ana_dir_path = args.ana_dir_path
discovery_thresh = args.discovery_thresh
compute_sig_trials = args.compute_sig_trials
cpus_used = args.cpus_used

wrkdir = args.wrkdir
if not os.path.exists(wrkdir):
    os.makedirs(wrkdir)
GeV = 1.        # base SkyLab and csky energy units are in GeV mostly
TeV = 1000*GeV    
E0    = pivot_energy*TeV                # pivot energy 
inj_enj_lo = (float(inj_energy[:3]))*TeV   # inj energy LL
inj_enj_hi = (float(inj_energy[4:]))*TeV   # inj energy UL 
print("Min max inj energy in GeV")
print(inj_enj_lo,inj_enj_hi)

################################################################################
# GET MOJAVE DATA SET (NAME RA DEC av_flux pk_flux)SET Directory
################################################################################
array = np.loadtxt(input_file,dtype='str')
src_name = np.asarray(array[:,0])
src_ra = np.asarray(array[:,1],dtype=float)*np.pi/180
src_dec = np.asarray(array[:,2],dtype=float)*np.pi/180
flux_used = np.asarray(array[:,3],dtype=float)
print("---------------------------------------------------------\n")
print("Number of sources used in stacking: %s"%(len(flux_used)))
print("---------------------------------------------------------\n")
no_of_sources = len(flux_used)

################################################################################
# Weights!
################################################################################
if equal_weights:
    wt = np.ones(len(flux_used))
else:
    wt = np.copy(flux_used)
wt=wt/np.sum(wt)
print("\n WEIGHTS USED FOR THIS RUN:\n %s \n"%(wt))

################################################################################
# CSKY SETUP and LOAD ICECUBE DATA
################################################################################
ana_dir = cy.utils.ensure_dir(ana_dir_path)
repo = cy.selections.repo
ana = cy.get_analysis(cy.selections.repo, cy.selections.PSDataSpecs.ps_10yr, dir=ana_dir)
if(os.path.isdir(ana_dir)==True):
    print("Not saving analysis cache as folder exists")
else:
    print("Saving analysis cache to ",ana_dir)
    ana.save(ana_dir)cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = cpus_used

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
    tr = cy.get_trial_runner(src=src, ana=ana,flux=flux, sindec_bandwidth=np.radians(.1), mp_cpus=cpus_used)

################################################################################
# BACKGROUND
################################################################################
bkgfile_comb= wrkdir+"/bkg_files/background_for_index_%s_energy_inj_%s.npy"%(index,energy_inj)
if(os.path.isfile(bkgfile_comb)==True):
    bg_arr=np.load(bkgfile_comb)
    bg = cy.dists.Chi2TSD(cy.utils.Arrays(init=bg_arr))
else:
    print("computing background:")
    bg_arr = tr.get_many_fits(nscramble, seed=random.randint(1, 1000))
    bg = cy.dists.Chi2TSD(bg_arr)
    np.save(bkgfile_comb,bg_arr.as_array)
print(bg.description)
                         
if make_background_plot:
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
    plt.savefig(wrkdir+"background_plot_for_index_%s_energy_bin_%s.png"%(index,energy_inj))

################################################################################
# SENSITIVITY
################################################################################
sens = tr.find_n_sig(
        bg.median(), 0.9, #90%
        n_sig_step=nstep,
        batch_size=nsample,
        max_batch_size=nsample,
        tol=0.05 #Change tolerance if required
        )

    
e_mid=10**((math.log((inj_enj_lo/1000),10)+math.log((inj_enj_hi/1000),10))/2)

ednde,ednde_err,e2dnde,e2dnde_err = computing_fluxes_from_ns(tr,sens,e_mid,1e3)
print("\nSensitivity Flux in TeV/cm2/s  @ %s TeV:"%(E0/TeV))
print(e2dnde,"+-",e2dnde_err)


# SAVE RESULTS TO COMMON FILE
# Column names E_range gamma e_norm sensflux(e2dnde)calc sensflux(e2dnde)err
out_line="%s\t%s\t%s\t%.2e\t%.2e\t%s"%(energy_inj,index,e_mid,e2dnde,e2dnde_err,"TeV/cm2/s")
combresult_file= wrkdir+"/sens_combined_result_file_index_%s.txt"%(index)
if(os.path.isfile(combresult_file)!=True):
    Write2File("E_range\tgamma\te_norm\tsensflux(E2dnde)calc\tsensflux(E2dnde)err",combresult_file)
Write2File(out_line,combresult_file)

################################################################################
# DISCOVERY POTENTIAL
################################################################################                                    
disc = tr.find_n_sig(bg.isf_nsigma(discovery_thresh), 
                     0.5, #50%
                     n_sig_step=nstep, 
                     batch_size=nsample, 
                     tol=0.05) #Change tolerance if required


ednde,ednde_err,e2dnde,e2dnde_err = computing_fluxes_from_ns(tr,disc,e_mid,1e3)
print("\n DP Flux in TeV/cm2/s  @ %s TeV:"%(E0/TeV))
print(e2dnde,"+-",e2dnde_err)


# SAVE RESULTS TO COMMON FILE
# Column names E_range gamma e_norm dpflux(e2dnde)calc dpflux(e2dnde)err
                
out_line="%s\t%s\t%s\t%.2e\t%.2e\t%s"%(energy_inj,index,e_mid,e2dnde,e2dnde_err,"TeV/cm2/s")
combresult_file=wrkdir+ "/dp_combined_result_file_index_%s.txt"%(index)
if(os.path.isfile(combresult_file)!=True):
    Write2File("E_range\tgamma\te_norm\tdpflux(E2dnde)calc\tdpflux(E2dnde)err",combresult_file)
Write2File(out_line,combresult_file)
