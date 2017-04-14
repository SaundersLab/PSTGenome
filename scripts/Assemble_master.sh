#!/bin/sh
#SBATCH --mem 4000
#SBTACH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p tgac-medium


cd /usr/users/ga004/buntingd/PSTGenome
source /usr/users/ga004/buntingd/PSTGenome/genome/bin/activate
export LUIGI_CONFIG_PATH=$SL/FP_pipeline/luigi.cfg

srun python Assemble.py pe_libs.txt lmp_libs.txt \
            --base-dir $SL/PSTGenome/data \
            --workers 100\
