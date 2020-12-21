#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --mem=8g
#SBATCH --time=2-00:00:00
#SBATCH --output=./test_output
#env
./run_main.sh /net/nfs/opt/matlab2020a/ 
