#!/bin/sh
#SBATCH --mem 8000
#SBTACH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p nbi-medium


cd /usr/users/ga004/buntingd/PSTGenome
source /usr/users/ga004/buntingd/PSTGenome/genome/bin/activate
export LUIGI_CONFIG_PATH=$SL/FP_project/FP_pipeline/luigi.cfg

srun python Assemble.py pe_libs.txt lmp_libs.txt \
            --K-list '[200, 216, 232, 240, 280, 300, 400, 500, 640]' \
            --base-dir $SL/PSTGenome/data \
            --workers 100
