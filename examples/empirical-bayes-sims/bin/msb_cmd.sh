#! /bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=2G
#$ -l vf=2G
#$ -q all.q
#$ -pe orte 10

if [ -n "$SGE_O_WORKDIR" ]
then
    source ~/.bash_profile
    cd /share/work1
    cd $SGE_O_WORKDIR
fi

nprocs=10
nreps=100
nprior=10000
npost=100

msb.py --np $nprocs -o ../configs/m5.cfg -p ../configs/*.cfg -r $nreps \
    -n $nprior --num-posterior-samples $npost --output-dir ../results \
    --merge-priors --keep-priors --seed 4849390

