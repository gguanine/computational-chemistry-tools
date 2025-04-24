#!/bin/bash

nprocs=24
mem="48GB"

# Remove existing input.sh if it exists
if [ -f ./input.sh ]; then
    rm ./input.sh
fi

# Create the initial input.sh PBS script
cat > input.sh <<-EOM
#!/bin/bash
#PBS -P gaussian
#PBS -m abe
#PBS -q parallel
#PBS -l select=1:ncpus=${nprocs}:mem=${mem}
#PBS -l walltime=720:00:00

cd \$PBS_O_WORKDIR

EOM

# Loop through all .gjf files and add g16 + formchk commands
for f in ./*.gjf; do
    # Skip if no .gjf files are found
    [ -e "$f" ] || continue

    filename=$(basename "$f" .gjf)
    echo "writing $filename"
    
    printf "g16 -p=${nprocs} -m=${mem} < ${filename}.gjf > ${filename}.log\n" >> input.sh
    printf "formchk -3 ${filename}.chk ${filename}.fchk\n\n" >> input.sh
done

# qsub ./input.sh
