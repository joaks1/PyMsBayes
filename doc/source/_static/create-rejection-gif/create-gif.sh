#! /bin/bash

python reject-sampling-plot.py

convert -loop 0 -delay 50 rejection-sampling-observed.png rejection-sampling-tolerance.png rejection-sampling-10.png rejection-sampling-200.png rejection-sampling-500.png rejection-sampling-1000.png rejection-sampling-10000.png rejection-sampling-20000.png -delay 250 rejection-sampling-20000-post.png ../rejection-sampling.gif

rm rejection-sampling*.png
