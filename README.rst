.. |date| date::

******************************
PDX Mouse Subtraction Pipeline
******************************

:authors: Komal Rathi
:contact: rathik@email.chop.edu
:organization: DBHi, CHOP
:status: Work in progress
:date: |date|

.. meta::
   :keywords: pdx, mouse, 2016
   :description: pdx mouse subtraction pipeline.

Introduction
============

The goal of this repo is to make the Mouse subtraction pipeline from BCM (Wheeler Lab) reproducible.

Installation
============

.. code-block:: bash

	# for all publicly available tools 
	conda create --name pdx-subtract-env
	conda activate pdx-subtract-env
	conda install -c bioconda samtools
	conda install -c bioconda sambamba
	conda install -c bioconda picard
	conda install -c bioconda cufflinks
	conda install -c bioconda rna-seqc
	conda install -c bioconda star
	conda install -c bioconda star-fusion
	conda install -c bioconda bwa
	conda install -c bioconda alignstats

	# soapfuse has to be installed separately
	# not available on conda
	wget https://sourceforge.net/projects/soapfuse/files/SOAPfuse_Package/SOAPfuse-v1.26.tar.gz
	tar -xzf SOAPfuse-v1.26.tar.gz
	cd SOAPfuse-v1.26

	# get soapfuse database
	wget http://public.genomics.org.cn/BGI/soap/SOAPfuse/hg19-GRCh37.59.for.SOAPfuse.tar.gz
	
	# update soapfuse config file according to http://soap.genomics.org.cn/soapfuse.html

Prepare reference fasta and gtf:
================================

.. code-block:: bash

	# Code to prepare reference fasta and gtf (But as a backup I have all those reference files from Maria as well):
	bash scripts/generate_ref.sh

BCM-specific scripts and software:
==================================

.. code-block:: bash

    1. call_htseq.sh
    2. run-defuse.sh
    3. pindel_0.2.5b5_tdonly
    4. ERCCPlot.jar
    5. RnaSeqLimsData.pl

Steps to run the pipeline:
==========================

1. Create config file (see config.yaml)
2. Use the command below to run the pipeline:

.. code-block:: bash

	source activate pdx-subtract-env
	snakemake -p -j 16 --nolock --snakefile Snakefile_phase1

