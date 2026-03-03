#!/bin/bash

#SBATCH --job-name=ctwas_test
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=24:00:00
#SBATCH --output= ######## CHANGE THIS ########.out
#SBATCH --error=  ######## CHANGE THIS ########.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jottensemeier@lji.org

#CHANGE THIS
WORKDIR="######## CHANGE THIS ########"
cd ${WORKDIR};

#DO NOT CHANGE THIS
profile="../profiles/slurm"

#CHANGE THIS
# config="../../tests/configs/Autoimmune.yaml"
logs="../runs"

# Remove necessary to make sure the run info is correct
# rm -r ${logs}
mkdir -p ${logs}

configfile="../../config/config.yaml"

snakemake --profile ${profile} \
          --configfile ${configfile} \
          --stats   ${WORKDIR}/${logs}/Test.stats >& ${logs}/Test.log