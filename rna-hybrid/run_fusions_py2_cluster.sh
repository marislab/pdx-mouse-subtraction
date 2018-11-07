source activate star-fusion-env
snakemake -s Snakefile_py2.py -p -j 3 --cluster-config cluster.yaml -c "qsub -cwd -e error.txt -o output.txt -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -l m_mem_free={cluster.m_mem_free} -pe smp {threads}" &
