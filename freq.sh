#!/bin/bash
#PBS -P gaussian
#PBS -m abe
#PBS -q parallel12
#PBS -l select=1:ncpus=12:mem=16GB
#PBS -l walltime=720:00:00

dir_opt="opt" # Directory for .gjf files for optimizations
dir_freq="freq" # Directory for .gjf files for frequency calculation
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

run_gaussian_in_dir $dir_freq

