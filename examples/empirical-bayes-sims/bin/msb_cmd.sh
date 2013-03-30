#!/bin/sh

nprocs = 10
nreps = 100
nprior = 10000
npost = 100

msb.py --np $nprocs -o ../configs/m5.cfg -p ../configs/*.cfg -r $nreps \
    -n $nprior --num-posterior-samples $npost --output-dir ../results \
    --merge-priors --keep-priors --seed 4849390

