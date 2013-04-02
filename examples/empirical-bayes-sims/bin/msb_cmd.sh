#! /bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=16G
#$ -l vf=16G
#$ -q all.q
#$ -pe orte 8

if [ -n "$SGE_O_WORKDIR" ]
then
    source ~/.bash_profile
    cd /share/work1
    cd $SGE_O_WORKDIR
    tmp_dir=$(mktemp -d /tmp/output.XXXXXXXXX)
fi

nprocs=8
nreps=10
nprior=1000
npost=100

msb.py --np $nprocs \
    -o ../configs/m5.cfg \
    -p ../configs/*.cfg \
    -r $nreps \
    -n $nprior \
    --num-posterior-samples $npost \
    --output-dir ../results \
    --staging-dir $tmp_dir \
    --merge-priors --keep-priors --seed 4849390

echo 'Here are the contents of the local temp directory '${$tmp_dir}':'
ls -Fla $tmp_dir
echo 'Removing the local temp directory...'
rm -r $tmp_dir
echo 'Done!'

