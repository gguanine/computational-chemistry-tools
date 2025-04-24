#!/bin/bash
#PBS -P gaussian
#PBS -m abe
#PBS -q parallel
#PBS -l select=1:ncpus=24:mem=48GB
#PBS -l walltime=720:00:00

nprocs=24
mem="48GB"
dir_opt="opt" # Directory for .gjf files for optimizations
dir_energy="energy" # Directory for .gjf files for energy calculation
command="# b3lyp def2svpp scrf=(SMD, solvent=acetonitrile) em=gd3"


# initialize HPC  
cd $PBS_O_WORKDIR;
source /etc/profile.d/rec_modules.sh
module load miniconda
bash
source ~/.bashrc
conda activate chemistry # your conda environment with required libraries.

run_gaussian_in_dir() {
    local dir="$1"
    cd "$PBS_O_WORKDIR/$dir" || { echo "Directory $dir not found."; return 1; }

    shopt -s nullglob
    gjf_files=(*.gjf)
    if [ ${#gjf_files[@]} -eq 0 ]; then
        echo "No .gjf files found in $dir"
        return 0
    fi

    for gjf_file in "${gjf_files[@]}"; do
        base_name=$(basename "$gjf_file" .gjf)
        echo "Running Gaussian for $base_name.gjf"
        g16 -p=${nprocs} -m=${mem} < "$gjf_file" > "$base_name.log"
        formchk -3 "$base_name.chk" "$base_name.fchk"
    done

    cd "$PBS_O_WORKDIR"
}


# Run Gaussian for .gjf files in $dir_opt
run_gaussian_in_dir $dir_opt

# Generate .gjf from previous optimization. Imaginary frequency is not allowed by default
python ./logs2gjfs.py -i "$dir_opt" -o "$dir_energy" --cmd "$command" # edit calculation levels here

# Run Gaussian for .gjf files in $dir_energy
run_gaussian_in_dir $dir_energy

