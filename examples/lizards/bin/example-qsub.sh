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
fi

staging_dir=$(mktemp -d /tmp/output.XXXXXXXXX)


nprocs=8
nprior=10000000
batch_size=12500
nsums=100000
npost=1000
nquantiles=1000
reporting_freq=100000
sort_index=0
seed=845225390

output_dir="../results/sample-output/lizard-analysis"
if [ ! -d "$output_dir" ]
then
    mkdir -p $output_dir
fi

dmc.py --np $nprocs \
    -o ../configs/dpp-simple.cfg \
    -p ../configs/dpp-simple.cfg ../configs/msbayes.cfg \
    -n $nprior \
    --prior-batch-size $batch_size \
    --num-posterior-samples $npost \
    --num-standardizing-samples $nsums \
    -q $nquantiles \
    --reporting-frequency $reporting_freq \
    --sort-index $sort_index \
    --output-dir $output_dir \
    --staging-dir $staging_dir \
    --temp-dir $staging_dir \
    --compress \
    --seed $seed

echo "Here are the contents of the local temp directory '${staging_dir}':"
ls -Fla $staging_dir
echo 'Removing the local temp directory...'
rm -r $staging_dir
echo 'Done!'

