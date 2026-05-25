#!/bin/bash -l

#PBS -N Run_PolarRes_SUMMARY_MARCESM
#PBS -l select=1:ncpus=4:mem=300GB:cpu_id=any
#PBS -l walltime=24:00:00
#PBS -j oe 
#PBS -J 1-15


###############################################
#
#
#  Display PBS info
#
#
###############################################
print_pbs_info(){
    echo ------------------------------------------------------
    echo -n 'Job is running on node '; cat $PBS_NODEFILE
    echo ------------------------------------------------------
    echo PBS: qsub is running on $PBS_O_HOST
    echo PBS: originating queue is $PBS_O_QUEUE
    echo PBS: executing queue is $PBS_QUEUE
    echo PBS: working directory is $PBS_O_WORKDIR
    echo PBS: execution mode is $PBS_ENVIRONMENT
    echo PBS: job identifier is $PBS_JOBID
    echo PBS: job name is $PBS_JOBNAME
    echo PBS: node file is $PBS_NODEFILE
    echo PBS: current home directory is $PBS_O_HOME
    echo PBS: PATH = $PBS_O_PATH
	echo PBS: sub-job is $PBS_ARRAY_INDEX
    echo ------------------------------------------------------

    # displaying some additional node info
    # is handy for debugging some things or know if you are
    # encountering any problematic nodes
    echo ""
    echo ------------------------------------------------------
}

###############################################
#
#
#  Helper/Setup Functions
#
#
###############################################


load_modules(){
        #module load mamba
		conda activate RLibs_V4.4.1D
		# module load r/4.4.1-gfbf-2023a
		# export R_LIBS="/mnt/home/n11222026/.conda/envs/RLibs_V4.4.1D/lib/R/library"
}
copy_in(){
    #copy some data to  your input directory
    #nothing to copy in on this script
    #For empty bash functions, must put a colon in them,
    #otherwise it will throw an error
    :
}


copy_out(){
    #nothing to copy out on this script
    #For empty bash functions, must put a colon in them,
    #otherwise it will throw an error
    :
}



run_program(){
    cd $PBS_O_WORKDIR
	Rscript PolarRes_Summarising_MAR_CESM2_ALL_VARIABLES_v2.R $PBS_ARRAY_INDEX
	
}


run_clean(){
    #nothing to clean for this script
    #For empty bash functions, must put a colon in them,
    #otherwise it will throw an error
    :
}

###############################################
#
#
#  Running everything
#
#
###############################################

print_pbs_info
load_modules
copy_in
copy_out
run_program
run_clean
