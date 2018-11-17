#!/bin/csh

#$ -M anematba@nd.edu	 # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  gpu 	 # Specify queue
#$ -l gpu_card=1
#$ -pe smp 4 
#$ -N  run_May11	 # Specify job name
#$ -t  1-8:1 

set param = ( N01_0 N01_1 N01_2 N01_3 N02_0 N02_1 N02_2 N02_3 )
module load gcc/4.9.2
module load cuda/7.0
module load bertini 
module load boost/1.58
echo -n "It is currently: ";date
echo -n "I am logged on as ";who am i
echo -n "This computer is called ";hostname
echo -n "I am currently in the directory ";pwd
#setenv PATH /afs/crc.nd.edu/user/a/anematba/Public/2015/Oct/11th/SceCells/bin:$PATH
./bin/runDiscSimulation_M -slurm $param[${SGE_TASK_ID}]
