#!/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import cosmolopy
import math
import scipy
from scipy import stats
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
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
####################################################
def doppler_factor(tau,theta):
    return (tau-math.cos(np.radians(theta))*((tau)**2.-1.)**(1./2.))**(-1.)
####################################################
cosmology = {'omega_M_0': 0.308, 'omega_lambda_0': 0.692,'omega_k_0': 0.0, 'h': 0.678}
####################################################
def lum_v(S,z):
    s_lum_a=S*4*np.pi*((3e24*cosmolopy.distance.luminosity_distance(z, **cosmology))**2)
    p=2+alpha
    return s_lum_a/(((doppler_factor(tau,theta))**p)*((1+z)**(1+alpha)))
####################################################
def e(z):
    return ((1+z)**k)*math.exp(z/eta)
####################################################
def dVdz(z):
    return  cosmolopy.distance.diff_comoving_volume(z, **cosmology)/((180/np.pi)**2)
####################################################
def phi(Lp,z):    
    A=1#dVdz(z)/math.log(Lp,10)#*1e5**1e-13/(4*3.141512) #mpc-3erg-1s ---> s/(erg deg2)
    L_star=1#2.61e46 #erg/sec
    L_by_z=Lp/e(z)
    #print(math.log(L_by_z,10))
    return (A/math.log(L_by_z,10))*(L_by_z/L_star)**gamma  #Max lum 1e33W/s --> 1e26*1e7 erg/s
####################################################
def dn(s_e):
    #dn_ds=scipy.integrate.quad(lambda z: (phi(z,lum_s(s_e,z))*dVdz(z)*lum_s(s_e,z)/3282.8 ),0.01,6) 
    dn_ds=scipy.integrate.quad(lambda z: (phi(lum_v(s_e,z),z)),0.01,6) 
    return dn_ds[0]
####################################################
def inner_integral(L1,z1):
    #print(L1)
    return scipy.integrate.quad(lambda lum: phi(lum,z1),L1,lum_max)[0] 
####################################################
def fitting_func(x, a,b):
    return x*a+b
####################################################
# PARAMETERS REQUIRED FOR INTEGRALS (DERIVED USING THE MOJAVE CATALOG)
k=8.0
eta=-0.35
gamma=-3.1
alpha=0.0
tau=10
theta=2

#############
# ARGUMENTS #
#############
p = argparse.ArgumentParser(description="Calculates COMPLETENESS fOR MOJAVE CATALOG",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--save_to_file", default=False,
               help="Save trial to file (default=False)")
p.add_argument("--plot_trial", default=False,
               help="Plot trial (default=False)")
p.add_argument("--input_file", default="../MOJAVE_multiyear_complete_sample.txt", type=str,
               help="Location of input file: MOJAVE_multiyear_complete_sample.txt")
args = p.parse_args()
save_to_file = args.save_to_file
plot_trial = args.plot_trial
input_file = args.input_file
    
########################################################################################################
# GET LORENTZ FACTOR AND VIEWING ANGLE DISTRIBUTIONS
########################################################################################################
lorentz_input_file=np.loadtxt("Lorentz_factor_dist_mojaveXVII.txt")
lorentz_dist_x=lorentz_input_file[:,0]
lorentz_dist_y=lorentz_input_file[:,1]

viewing_angle_input_file=np.loadtxt("viewing_angle_dist_mojaveXVII.txt")
viewing_ang_dist_x=viewing_angle_input_file[:,0]
viewing_ang_dist_y=viewing_angle_input_file[:,1]

lorentz_dist_y=lorentz_dist_y/sum(lorentz_dist_y)
viewing_ang_dist_y=viewing_ang_dist_y/sum(viewing_ang_dist_y)

#SAMPLE LORENTZ FACTOR
f = interp1d(lorentz_dist_x,lorentz_dist_y)
xnew = np.linspace(1.55, 40., num=25, endpoint=True) ##PARAMETES FIXED TO MAKE SURE THE DISTRIBUTION IS REPRODECED, Can be edited as required
ynew=[]
for i in range(len(xnew)):
    ynew.append(f(xnew[i]))
ynew=np.asarray(ynew)
ynew=ynew/np.sum(ynew)
gamma_samples=(np.random.choice(xnew,100, p=ynew))

#SAMPLE VIEWING ANGLE    
f2 = interp1d(viewing_ang_dist_x,viewing_ang_dist_y)
xnew2 = np.linspace(0.3, 30.00, num=55, endpoint=True) ##PARAMETES FIXED TO MAKE SURE THE DISTRIBUTION IS REPRODECED, Can be edited as required
ynew2=[]
for i in range(len(xnew2)):
    ynew2.append(f2(xnew2[i]))
ynew2=np.asarray(ynew2)
ynew2=ynew2/np.sum(ynew2)
viewing_ang_sample=(np.random.choice(xnew2,100, p=ynew2))

########################################################################################################
# GET MOJAVE SAMPLE FOR COMPARISION 
########################################################################################################
array = np.loadtxt(input_file,dtype='str')
flux_used=[]
name=[]
for i in range(len(array[:,3])):
    if float(array[i,3])<=45.0:   #REMOVING LOWEST FLUX AS IT IS OUT OF RANGE TO BE TESTED
        print("removing flux value of",array[i,3])
        continue
    flux_used.append(array[i,3])
    name.append(array[i,0])
flux_used=np.asarray(flux_used,dtype="float")
#print(np.sort(flux_used))
mojave_flux_eg_prcm_s=flux_used*(1e-3*1e-23*15*1e9) #Convert mJy to erg sec-1 cm-2 (15GHz used)
mojave_flux_jy=flux_used*(1e-3) #Convert mJy to Jy (15GHz used)
max_flux=np.amax(mojave_flux_eg_prcm_s)#50*(1e-3*1e-23*15*1e9)

print("Maximum MOJAVE Flux:",max_flux)
lum_max=lum_v(max_flux,6) #1e31*1e7*15*1e9 # W/Hz --> erg/s/Hz---> erg/sec
print("Maximum MOJAVE Luminosity:",lum_max)
min_flux=np.amin(mojave_flux_eg_prcm_s)
print("Minimum MOJAVE Flux:",min_flux)
lum_min=lum_v(min_flux,6) #1e31*1e7*15*1e9 # W/Hz --> erg/s/Hz---> erg/sec
print("Minimum MOJAVE Luminosity:",lum_min)

########################################################################################################
# COMPUTE LOGN-LOGS 
########################################################################################################
dn_arr_fin=[]
dn_arr_all=[]
flux_arr=np.logspace(math.log(max_flux,10)-0.001,-22,100) #Length of x axis (log S)

#For EACH SKY SIMULATION which runs len(viewing_ang_sample) times, COMPUTE THE dN/dS and save it in dn_arr_all. 
for j in range(len(viewing_ang_sample)):    
    tau=gamma_samples[j]
    theta=viewing_ang_sample[j]
    dn_arr=[]
    print("running",j)
    for i in range(len(flux_arr)):            
        dn_arr.append(scipy.integrate.quad(lambda z: inner_integral(lum_v(flux_arr[i],z),z),0.01,6)[0]) 
    dn_arr_all.append(np.asarray(dn_arr))
dn_arr_all=np.asarray(dn_arr_all)
dn_arr_fin=[]
dn_arr_std=[]
for i in range(len(flux_arr)):  
    dn_arr_fin.append(np.mean(dn_arr_all[:,i]))
    dn_arr_std.append(np.std(dn_arr_all[:,i]))
dn_arr_fin=np.asarray(dn_arr_fin)

#Get Limits based on the minimum and maximum of the lorentz distribution 
dn_array_hi=[]
dn_arr_low=[]
tau=np.amin(lorentz_dist_x)#lorentz_sample[i]
theta=np.amin(viewing_ang_dist_x)#viewing_ang_sample[i]
dn_arr=[]
for i in range(len(flux_arr)):            
        dn_arr=(scipy.integrate.quad(lambda z: inner_integral(lum_v(flux_arr[i],z),z),0.01,6)[0]) 
        dn_arr_low.append(dn_arr)
tau=np.amax(lorentz_dist_x)#lorentz_sample[i]
theta=np.amax(viewing_ang_dist_x)#viewing_ang_sample[i]
dn_arr=[]
for i in range(len(flux_arr)):            
        dn_arr=(scipy.integrate.quad(lambda z: inner_integral(lum_v(flux_arr[i],z),z),0.01,6)[0]) 
        dn_array_hi.append(dn_arr)

#Compute N>S for the y axis of logN-logS
n_gr_than_S=[]
flux_bins=np.logspace(math.log(max_flux,10)-0.001,math.log(min_flux,10),40)
for i in range(len(flux_bins)):
    num = mojave_flux_eg_prcm_s > flux_bins[i]
    n_gr_than_S.append(len(mojave_flux_jy[num]))    
n_gr_than_S=np.asarray(n_gr_than_S,dtype="float")
print("N>S: ",n_gr_than_S)
print("Flux bins: ",flux_bins)

n_gr_than_S_up=n_gr_than_S*(1/(4*3.14151)) #Dividing by 4pi will take into account the completeness due to sky position also

#NORMALIZE
constant=(dn_arr_fin[0])/(1/(4*3.141512))#/n_gr_than_S_up[0]
#print(constant)
dn_arr_fin_up=np.asarray(dn_arr_fin)/constant 
dn_arr_std_up=np.asarray(dn_arr_std)/((dn_arr_std[0])/(1/(4*3.141512)))

########################################################################################################
# MODIFY LOGN-LOGS to run to a point when no new "low flux" sources are simulated 
########################################################################################################
new_flux_arr=[]
new_dn_arr=[]
new_dn_arr_std_up=[]
for i in range(len(flux_arr)):
    if flux_arr[i]<=1e-19 and dn_arr_fin_up[i]-dn_arr_fin_up[i-1] <=0.7:
        print(flux_arr[i])
        continue
    new_flux_arr.append(flux_arr[i])
    new_dn_arr.append(dn_arr_fin_up[i])
    new_dn_arr_std_up.append(dn_arr_std_up[i])    
new_dn_arr=np.asarray(new_dn_arr)
new_flux_arr=np.asarray(new_flux_arr)
new_dn_arr_std_up=np.asarray(new_dn_arr_std_up)

########################################################################################################
# Compute Completeness 
########################################################################################################
mojave=scipy.integrate.simps(n_gr_than_S_up,x=-flux_bins)
sim=scipy.integrate.simps(new_dn_arr,x=-new_flux_arr)
print("Area under the curve for MOJAVE and simulated sources: ",mojave,sim)
print("Ratio: ",mojave/sim)
result_final=mojave/sim

########################################################################################################
# Compute FIT log(y)=m*log(x)+b 
########################################################################################################
fitting_x=[]
fitting_y=[]
for i in range(len(flux_bins)):
    if flux_bins[i]>=1e-13:
        fitting_x.append(math.log(flux_bins[i],10))
        fitting_y.append(math.log(n_gr_than_S_up[i],10))
x=np.logspace(math.log(flux_bins[1],10),-14,100)
popt1_1, pcov1_1 = curve_fit(fitting_func, fitting_x,fitting_y)
plt.plot( fitting_x,fitting_y)
y=[]
for i in range(len(x)):
    logy=fitting_func(math.log(x[i],10),popt1_1[0],popt1_1[1])
    y.append(10**logy)
perr1_1 = np.sqrt(np.diag(pcov1_1)) #1 sigma errors
print("FITTING RESULTS:")
print(x,y)
print(popt1_1, pcov1_1,perr1_1)

########################################################################################################
# PLOT
########################################################################################################
if plot_trial:
    plt.plot(flux_bins,n_gr_than_S_up,".",label="Mojave Sources")
    plt.plot(new_flux_arr,new_dn_arr,".",label="Simulated Source Count")
    plt.fill_between(new_flux_arr,new_dn_arr-new_dn_arr_std_up,new_dn_arr+new_dn_arr_std_up,color="orange",alpha=0.2)
    plt.plot(x,y,label="logN=(%.2f$\pm$%.2f)*log$_{10}$(S)+(%.1f$\pm$%.1f)"%(popt1_1[0],perr1_1[0],popt1_1[1],perr1_1[1])) 
    plt.legend(loc="upper left",fontsize=12)
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(1e-1,2e5)
    plt.grid()
    plt.xlabel("Flux log$_{10}$S erg/sec/cm$^2$",fontsize=14)
    plt.ylabel("N>S str$^{-1}$",fontsize=14)
    plt.savefig("LogNlogS_Mojave_complete.png",bbox_inches='tight')

########################################################################################################
# RUN ONE TRIAL AND SAVE TO FILE
########################################################################################################
if save_to_file:
    filename= "result_file_all.txt"
    Write2File(str(result_final),filename)

