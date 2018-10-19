source activate pdx-subtract-env
snakemake --rerun-incomplete -p -j 16 --nolock --snakefile Snakefile_phase2.py
