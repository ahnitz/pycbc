`which python` -m cProfile -o log /home/ahnitz/projects/pycbc/bin/inference/pycbc_inference \
--config-file `dirname "$0"`/single.ini \
--nprocesses=1 \
--output-file single.hdf \
--seed 0 \
--force \
--verbose
