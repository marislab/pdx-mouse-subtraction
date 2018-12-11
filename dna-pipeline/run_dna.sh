source activate pdx-subtract-env
snakemake --rerun-incomplete -p -j 16 --nolock --snakefile Snakefile_dna.py
