#!/bin/bash

if [ -f ./input.sh ]; then
    rm ./input.sh
fi

cat > input.sh <<-EOM
#!/bin/bash
#PBS -P gaussian
#PBS -m abe
#PBS -q parallel12
#PBS -l select=1:ncpus=12:mem=16GB
#PBS -l walltime=720:00:00

cd \$PBS_O_WORKDIR

EOM

for f in ./*.gjf
do
    filename=$(basename $f .gjf)
    echo "wtriting $filename"
    printf "g16 < %s.gjf > %s.log \nformchk -3 %s.chk > %s.fchk\n\n" $filename $filename $filename $filename >> input.sh
done

qsub ./input.sh
