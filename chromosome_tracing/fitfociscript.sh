#!/bin/bash

#SBATCH --job-name=fitfoci                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=100:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=myrun.o%j     #These are the STDOUT and STDERR files
#SBATCH --error=myrun.e%j

#You selected an array of jobs from 1 to 10 with 10 simultaneous jobs
#SBATCH --array=1-10
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ahilya.sawh@unibas.ch        #You will be notified via email when your task ends or fails




cd /scicore/home/mangos/finatu15/200625_pgl/

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.



module load MATLAB




matlab -nodisplay -r "ANSrun_fitfoci_nofig_denoised(${SLURM_ARRAY_TASK_ID})"


