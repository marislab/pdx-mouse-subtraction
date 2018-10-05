source activate pdx-subtract-env
snakemake -p -j 16 --nolock --snakefile Snakefile_phase1
