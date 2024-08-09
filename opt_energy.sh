#!/bin/bash
#PBS -P gaussian
#PBS -m abe
#PBS -q parallel12
#PBS -l select=1:ncpus=12:mem=16GB
#PBS -l walltime=720:00:00

dir_opt="opt" # Directory for .gjf files for optimizations
dir_energy="energy" # Directory for .gjf files for energy calculation


# initialize HPC  
cd $PBS_O_WORKDIR;
source /etc/profile.d/rec_modules.sh
module load miniconda
bash
source ~/.bashrc
conda activate chemistry # your conda environment with required libraries.

run_gaussian_in_dir() {
    local dir=$1
    cd "$PBS_O_WORKDIR/$dir"
    for gjf_file in *.gjf; 
    do
        base_name=$(basename "$gjf_file" .gjf)
        g16 <"$gjf_file" > "$base_name.log"
        formchk -3 "$base_name.chk" "$base_name.fchk"
    done   
    cd $PBS_O_WORKDIR
}

# PART 1: Generate .gjf files from SMILES
python ./smiles_for_optimize.py # edit your molecule SMILES, properties, and calculation levels here

# PART 2: Run Gaussian for .gjf files in $dir_opt
run_gaussian_in_dir $dir_opt

# PART 3: Generate .gjf from previous optimization. Imaginary frequency is not allowed by default
python ./opt_to_energy.py # edit calculation levels here

# PART 4: Run Gaussian for .gjf files in $dir_energy
run_gaussian_in_dir $dir_energy

