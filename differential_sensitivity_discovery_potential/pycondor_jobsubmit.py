import pycondor
from pycondor import *
import os,sys
from random import randint
import numpy as np

command_line = sys.argv


################################################################################
# SOURCE INPUT
################################################################################

index=3.0
E0=100 #TeV
e_range_arr=["0.1_1e0","1e0_1e1","1e1_1e2","1e2_1e3","1e3_1e4","1e4_1e5","1e5_1e6","1e6_1e7"]
#e_range="const"
index_arr=[2.0,2.5,3.0]
for e_range in e_range_arr:
    #for index in index_arr:
    mem_req=30000
    cpu_req=10
    file_loc="/scratch/adesai/"
    error=file_loc
    output=file_loc
    log=file_loc
    submit=file_loc


    job = pycondor.Job ('csky_dp_t_indep_%s_index_%s'%(e_range,index),#name='csky_tindep_index_%s'%(index), 
                executable='/data/user/adesai/csky/mojave/mojave_time_integrated_analysis_gthub_rep/differential_sensitivity_discovery_potential/mojave_time_independent_analysis_variable_einj.py',
                #executable='/data/user/adesai/csky/mojave/time_integrated_analysis/mojave_time_independent_analysis_fixed_einj.py',
                error=error,
                output=output,
                log=log,
                submit=submit,
                getenv=True,
                universe='vanilla',
                verbose=2,
                request_memory= mem_req,  
                request_cpus = cpu_req,
                extra_lines=['priority = 20','should_transfer_files = YES',  'Requirements =  (Machine != "node128.icecube.wisc.edu")'])
                #extra_lines= ['request_memory = (NumJobStarts is undefined) ? 2 * pow(2, 10) : 2048 * pow(2, NumJobStarts + 1) periodic_release = (HoldReasonCode =?= 21 && HoldReasonSubCode =?= 1001) || HoldReasonCode =?= 21 periodic_remove = (JobStatus =?= 5 && (HoldReasonCode =!= 34 && HoldReasonCode =!= 21)) || (RequestMemory > 13192)','should_transfer_files = YES',  'Requirements =  (Machine != "node128.icecube.wisc.edu")']   )

    
    job.add_arg('--pivot_energy=%s --inj_energy=%s --index=%s --cpus_used=10'%(E0,e_range,index))
    #job.add_arg('%s %s'%(E0,index))

    job.build_submit()



    
