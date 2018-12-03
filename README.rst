.. |date| date::

******************************
PDX Mouse Subtraction Pipeline
******************************

:authors: Oliver Hampton, Chase Miller, Liu Xi, Maria Cardenas
:contact: Komal Rathi (rathik@email.chop.edu)
:organization: DBHi, CHOP
:status: Completed
:date: |date|

.. meta::
   :keywords: pdx, mouse, 2016
   :description: pdx mouse subtraction pipeline.

Introduction
============

The goal of this repo is to make the Mouse subtraction pipeline from BCM (Wheeler Lab) reproducible.

Installation
============

1. Create python3 environment:

.. code-block:: bash

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
	conda install -c bioconda bamutil

2. Create python2 environment (STAR-Fusion v1.0.1 is python 2.7 compatible):

.. code-block:: bash

	conda create --name star-fusion-env python=2.7
	source activate star-fusion-env
	conda install -c bioconda star-fusion
	conda install -c bioconda trinity
	conda install -c conda-forge -c bioconda samtools bzip2
	conda install -c conda-forge configparser

	# install some non-standard perl modules:
	perl -MCPAN -e shell
	install DB_File
	install URI::Escape
	install Set::IntervalTree
	install Carp::Assert
	install JSON::XS

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

	# change PA_all_fq_postfix in config file to .fq

deFUSE
======

.. code-block:: bash

	# for deFUSE, python 2 is required so use the python2 environment created for STAR-Fusion

	# Install via source:
	wget https://bitbucket.org/dranew/defuse/get/0f198c242b82.zip
	unzip 0f198c242b82.zip
	
	# in the tools directory, download boost
	cd tools && wget https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.tar.gz
	tar -zxvf boost_1_68_0.tar.gz
	export CPLUS_INCLUDE_PATH=/mnt/isilon/maris_lab/target_nbl_ngs/PPTC-PDX-genomics/mouse_subtraction_pipeline/scripts/dranew-defuse-0f198c242b82/tools/boost_1_68_0
	cd tools && make

	# download deFUSE reference database
	# change perl in defuse_create_ref.pl to /usr/bin/env perl
	defuse_create_ref.pl -d /mnt/isilon/cbmi/variome/reference/defuse_db/hg19/

GATK
====

.. code-block:: bash
	
	# get reference files
	wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
	wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
	wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz

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

    1. pindel_0.2.5b5_tdonly
    2. ERCCPlot.jar
    3. RnaSeqLimsData.pl

Steps to run the RNA-pipeline:
==============================

The RNA pipeline is divided into four steps:

1. Snakefile_Phase1: Align PDX RNA-seq data to hybrid genome, split into human and mouse bams and create human specific fastq files.
2. Snakefile_Phase2: Realign to human reference, do QC, run htseq and pindel. 
3. Snakefile_fusions_py2: Run python2 dependent fusion callers like STAR-Fusion and deFUSE 
4. Snakefile_soapfuse: Run python3 dependent fusion caller like SOAPfuse

Each snakefile has a corresponding bash script to run the pipeline:

.. code-block:: bash
	
	# Run phase 1
	cd rna-hybrid && bash run_phase1.sh

	# Run phase 2
	cd rna-hybrid && bash run_phase2.sh

	# Run python2 based fusion callers
	cd rna-hybrid && bash run_fusions_py2.sh

	# Run python3 based fusion callers
	cd rna-hybrid && bash run_soapfuse.sh


Steps to run the DNA-pipeline:
==============================

.. code-block:: bash
	
	cd dna-pipeline && bash run_dna.sh



