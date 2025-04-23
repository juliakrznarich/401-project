#!/bin/bash
#SBATCH --job-name=openmm_md
#SBATCH --partition=your_partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:30:00
#SBATCH --mem=4GB
#SBATCH --output=md.out
#SBATCH --error=md.err

# Use the Python from your Conda environment directly
/mnt/home/krznari4/.conda/envs/openmm_env/bin/python -c "import openmm; print('✅ OpenMM successfully imported!')"

# Now run your script using that Python
/mnt/home/krznari4/.conda/envs/openmm_env/bin/python /mnt/ffs24/home/krznari4/401_project/openmm_example.py
