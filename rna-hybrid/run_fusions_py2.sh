source activate star-fusion-env
snakemake --rerun-incomplete -p -j 16 --nolock --snakefile Snakefile_py2.py
