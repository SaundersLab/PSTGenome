#!/bin/sh
#SBATCH --mem 12000
#SBTACH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p tgac-medium


cd /usr/users/ga004/buntingd/PSTGenome
source /usr/users/ga004/buntingd/PSTGenome/genome/bin/activate
export LUIGI_CONFIG_PATH=$SL/FP_pipeline/luigi.cfg

srun python Assemble.py pe_libs.txt lmp_libs.txt \
            --K-list '[80, 100, 160, 200, 240, 300, 400, 500, 640]' \
            --base-dir $SL/PSTGenome/data \
            --workers 100\
            --no-lock
