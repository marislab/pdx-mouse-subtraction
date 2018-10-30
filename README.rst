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
	conda install -c bioconda htslib
	conda install -c bioconda sambamba
	conda install -c bioconda picard
	conda install -c bioconda cufflinks
	conda install -c anaconda java-1.7.0-openjdk-cos6-x86_64 # required by rna-seqc
	conda install -c bioconda rna-seqc
	conda install -c bioconda htseq
	conda install -c bioconda star=2.5.3a
	conda install -c bioconda trinity=2.5.1 # required by star-fusion
	conda install -c bioconda star-fusion=1.1.0
	conda install -c bioconda bwa
	conda install -c bioconda alignstats
	conda config --add channels https://conda.anaconda.org/dranew
	conda install defuse

SOAPfuse
========

.. code-block:: bash

	# SOAPfuse has to be installed separately as it is not available on conda
	wget https://sourceforge.net/projects/soapfuse/files/SOAPfuse_Package/SOAPfuse-v1.26.tar.gz
	tar -xzf SOAPfuse-v1.26.tar.gz
	cd SOAPfuse-v1.26

	# get SOAPfuse database
	cd /mnt/isilon/cbmi/variome/reference/soapfuse_db
	wget http://public.genomics.org.cn/BGI/soap/SOAPfuse/hg19-GRCh37.59.for.SOAPfuse.tar.gz
	
	# update SOAPfuse config file according to http://soap.genomics.org.cn/soapfuse.html
	# add cytoBand file from ucsc and update SOAPfuse config
	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz hg19-GRCh37.59/
	gunzip cytoBand.txt.gz

STAR-Fusion
===========

.. code-block:: bash

	# for STAR-Fusion, python 2 is required so create a separate environment
	https://github.research.chop.edu/rathik/star-fusion-detection-pipeline

deFUSE
======

.. code-block:: bash

	# for deFUSE, python 2 is required so create a separate environment
	# change perl in defuse_create_ref.pl to /usr/bin/env perl
	# download deFUSE reference database
	defuse_create_ref.pl -d /mnt/isilon/cbmi/variome/reference/defuse_db/hg19/

Prepare reference fasta and gtf:
================================

.. code-block:: bash

	# Code to prepare reference fasta and gtf (this might be inaccurate because I got the reference files from BCM):
	bash scripts/generate_ref.sh

	# make sure all reference fasta files are indexed: 
	samtools faidx <file.fasta|file.fa>

	# make sure the fasta reference used by bwa is indexed:
	bwa index protein_coding_canonical.T_chr.fa

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

