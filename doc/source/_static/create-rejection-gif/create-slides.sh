#! /bin/bash

tmp_py_file="reject-sampling-plot-slides.py"
sed "s/\.png/\.pdf/g" reject-sampling-plot.py > $tmp_py_file
python "$tmp_py_file" 

pdftk rejection-sampling-observed.pdf rejection-sampling-tolerance.pdf rejection-sampling-10.pdf rejection-sampling-200.pdf rejection-sampling-500.pdf rejection-sampling-1000.pdf rejection-sampling-10000.pdf rejection-sampling-20000.pdf cat output ../rejection-sampling.pdf

rm "$tmp_py_file"
rm rejection-sampling*.pdf
