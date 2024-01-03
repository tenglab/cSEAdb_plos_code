#!/bin/bash

#PBS -N fc
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=8,mem=20gb
#PBS -m bea
#PBS -M Xiang.Liu@moffitt.org


# run job parallelly on HPC using GNU parallel tool

export PROJECT_DIR="/share/lab_teng/xiangliu/super_enhancer/se_dataportal/NCI_60/"
export PARALLEL_WLIST="$PROJECT_DIR/01_bowtie2_wlist"
#export ROSE_DIR="/home/4472271/rose"

#module load samtools
#module load R/3.6.3
#module load macs/2.1.0

#source ~/.bashrc
cd $ROSE_DIR

parallel -j 4 -k < $PARALLEL_WLIST
