source activate pdx-subtract-env
snakemake -p -j 4 --nolock --snakefile Snakefile_fusions.py
