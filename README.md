This folder includes scripts for the Radio Neutrino Correlation analysis using MOJAVE sources (time averaged case)
> Analysis wiki:https://wiki.icecube.wisc.edu/index.php/Radio_neutrino_correlation

> contact: Abhishek Desai (desai25@wisc.edu)

> IceCube environment  : py3-v4.1.0 (python 3)

> Csky  version        : v1.1.0

> Pycondor version     : v0.5.0

> Datasets             : 10 yr PointSourceTracks_v003p02

> Method               : Likelihood stacking (with weights)

The python 3 environment can be setup using:
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`

The required csky software can be installed and setup using the instructions given on:
https://user-web.icecube.wisc.edu/~mrichman/docs/csky/installation.html


# There are 4 main folders in the analysis.

## 1. calc_sensitivity_discovery_potential:

  This folder contains the code required to run the sensitivity and discovery potential calculations and their results for the complete energy range case.
  The "mojave_time_independent_analysis_fixed_einj.py" is the main code used for running the calculation while the result_dir is the directory where the results are stored.
  An example to run the main code is:
  > python mojave_time_independent_analysis_fixed_einj.py --index=2.0 --pivot_energy=100.0 --cpus_used=3

  Other options include:
  --wrkdir,--equal_weights,--nscramble,--nsample,--make_background_plot,--inj_energy,--nstep,--input_file,--ana_dir_path,--discovery_thresh,--cpus_used
  

  The Sensitivity and Discovery Potential for different indices plots shown on the wiki page are made using the results that are outputted in the textfile derived using the above code. 

  The pycondor_jobsubmit_fixed_einj.py is used to submit multiple jobs to the submitter. (To make changes to the path of the code and/or other input parameters modify the pycondor code directly)

  
## 2. differential_sensitivity_discovery_potential:

  This folder contains the code required to run the differential sensitivity and discovery potential calculations and their results for the complete energy range case.
  The "mojave_time_independent_analysis_fixed_einj.py" is the main code used for running the calculation while the result_dir is the directory where the results are stored.
  An example to run the main code is:
  > python mojave_time_independent_analysis_variable_einj.py --inj_energy=1e1_1e2 --index=2.0 --pivot_energy=100.0 --cpus_used=3

  Other options include:
  --wrkdir,--equal_weights,--nscramble,--nsample,--make_background_plot,--inj_energy,--nstep,--input_file,--ana_dir_path,--discovery_thresh,--cpus_used
  
  The Differential Sensitivity and Discovery Potential for different indices plots shown on the wiki page are made using the results that are outputted in the textfile derived using the above code. 

  The pycondor_jobsubmit_fixed_einj.py is used to submit multiple jobs to the submitter. (To make changes to the path of the code and/or other input parameters modify the pycondor code directly)

## 3. n_bias_test:

  This folder contains the code required to run the bias test.
  The "mojave_time_independent_analysis_n_bias_test.py" is the main code used for running the calculation 
  An example to run the code is:
  > python mojave_time_independent_analysis_n_bias_test.py --index=2.0 --pivot_energy=100.0

  Other options include:
  --wrkdir,--ana_dir_path
  
  The bias plots for different indices using csky and variable gamma, shown on the wiki page, are made using the above code. 

## 4. completeness_calculation :

  This folder contains the code required to run the completeness calculation.
  This code simulates the sky multiple times to get the completeness value for one trial. This can be plotted using the code: 
  > python simulate_source_count_one_time.py  --plot_trial=True

  To derive the final completeness value, you run the following to get one trial that is saved in a result file. Running this multiple times gives you a distribution of the completeness values which can be then used to get the mean value along with the standard deviation. To see how the resultfile will look like after multiple runs check "sample_result_file.txt". To see how the histogram of the sample_result_file looks like using plt.hist see "sample_result_dist.png" 
  > python simulate_source_count_one_time.py  --save_to_file=True
  
  Other options include:
  --input_file
