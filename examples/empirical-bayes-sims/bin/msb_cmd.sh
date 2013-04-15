#! /bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=16G
#$ -l vf=16G
#$ -q all.q
#$ -pe orte 8

tmp_dir="../results"
if [ -n "$SGE_O_WORKDIR" ]
then
    source ~/.bash_profile
    cd /share/work1
    cd $SGE_O_WORKDIR
    tmp_dir=$(mktemp -d /tmp/output.XXXXXXXXX)
fi

nprocs=4
nreps=1000
nprior=1000000
npost=1000
seed=48492341

msb.py --np $nprocs \
    -o ../configs/m5.cfg \
    -p ../configs/m[123478].cfg \
    -r $nreps \
    -n $nprior \
    --num-posterior-samples $npost \
    --rejection-tool msreject \
    --regression-method glm \
    --output-dir ../results \
    --staging-dir $tmp_dir \
    --merge-priors --keep-priors --seed $seed

if [ -n "$SGE_O_WORKDIR" ]
then
    echo 'Here are the contents of the local temp directory '${$tmp_dir}':'
    ls -Fla $tmp_dir
    echo 'Removing the local temp directory...'
    rm -r $tmp_dir
fi
echo 'Done!'

