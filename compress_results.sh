#!/usr/bin/bash

# In parallel, using 60 cores, compress all the results. Use the
# "notebooks/compress.py" script; it take two arguments: input (original) and
# output (compressed) file paths. Apply this script to every file with the .nc
# extension in the results directory recursively, moving results/foo/bar.nc to
# results_compressed/foo/bar.nc.

find -L results -name "*.nc" | xargs -P 60 -I {} \
    bash -c 'mkdir -p "results_compressed/$(dirname {})" && \
    python notebooks/compress.py {} "results_compressed/{}"'