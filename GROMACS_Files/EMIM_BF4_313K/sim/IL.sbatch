#!/bin/bash

#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
 
module load gromacs

NP=$SLURM_NPROCS

export inputfile='minimize'

cd ../

mkdir $inputfile

cd  $inputfile

gmx grompp -f ../sim/$inputfile.mdp -c ../sim/init_config.pdb -p ../sim/data.top -o $inputfile.tpr -maxwarn 2

mdrun_mpi -s $inputfile.tpr  -deffnm $inputfile.tpr > $inputfile.out

cd ../

export inputfile='annealing'

mkdir $inputfile

cd $inputfile

gmx grompp -f ../sim/$inputfile.mdp -po $inputfile.out.mdp -p ../sim/data.top -pp $inputfile.processed.top -c ../minimize/minimize.tpr.gro -o $inputfile.tpr -maxwarn 2

mpirun -np $NP mdrun_mpi -s $inputfile.tpr  -deffnm $inputfile.tpr > $inputfile.out



cd ../


export inputfile='nvt_eq'

mkdir $inputfile

cd $inputfile

gmx grompp -f ../sim/$inputfile.mdp -po $inputfile.out.mdp -p ../sim/data.top -pp $inputfile.processed.top -c ../annealing/annealing.tpr.gro -o $inputfile.tpr -maxwarn 2

mpirun -np $NP mdrun_mpi -s $inputfile.tpr  -deffnm $inputfile.tpr > $inputfile.out


cd ../

export inputfile='npt_eq'

mkdir $inputfile

cd $inputfile

gmx grompp -f ../sim/$inputfile.mdp -po $inputfile.out.mdp -p ../sim/data.top -pp $inputfile.processed.top -c ../nvt_eq/nvt_eq.tpr.gro -o $inputfile.tpr -maxwarn 2

mpirun -np $NP mdrun_mpi -s $inputfile.tpr  -deffnm $inputfile.tpr > $inputfile.out




cd ../

export inputfile='npt_prod1'

mkdir $inputfile

cd $inputfile

gmx grompp -f ../sim/$inputfile.mdp -po $inputfile.out.mdp -p ../sim/data.top -pp $inputfile.processed.top -c ../npt_eq/npt_eq.tpr.gro -o $inputfile.tpr -maxwarn 2

mpirun -np $NP mdrun_mpi -s $inputfile.tpr  -deffnm $inputfile.tpr > $inputfile.out



cd ../

export inputfile='npt_prod2'

mkdir $inputfile

cd $inputfile

gmx grompp -f ../sim/$inputfile.mdp -po $inputfile.out.mdp -p ../sim/data.top -pp $inputfile.processed.top -c ../npt_prod1/npt_prod1.tpr.gro -o $inputfile.tpr -maxwarn 2

mpirun -np $NP mdrun_mpi -s $inputfile.tpr  -deffnm $inputfile.tpr > $inputfile.out



