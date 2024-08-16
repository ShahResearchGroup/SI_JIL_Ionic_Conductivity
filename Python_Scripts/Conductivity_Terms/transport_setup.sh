#!/bin/bash

#Script was developed on 03/14/2023 by Amey Thorat and Ashutosh Verma

#this script converts trajectory files for center of mass of each molecule from gromacs tpr file. It implements bash multiprocessing
# step 1 creates individual index files for every molecule
# step 2 converts atomwise trajectory file to center of mass co-ordinates per molecule at each time step

# user needs to specify total number of species and total of molecules (can be estimated automatically as well)

export n_species="2"	#enter number of species present in the system

export n_total=""	#enter a number molecules here if you want to manually specify number of molecules. Leave "" if to be read from .gro file

export n_p_threads="32"

#read the total number of molecules from -2 (last but one) line of the *.gro file and store it into n_total

if [[ $n_total = "" ]] ; then

	b=( $(tail -2 *.gro) )	#takes last two lines into variable b

	n=$(echo ${b[0]} | grep -o .)	#n is b[0] split into characters

	for i in $n
	do
 		re='^[0-9]+$'			# numbers between 0 to 9. 
 		if [[ $i =~ $re ]] ; then	#if character is a number	
    		n_total+=$i			#concatanate
 		else				#else exit loop
    		break
 		fi		
	done
fi

echo "Total number of molecules: $n_total"

#creating directories to store individual index files and center of mass files for each molecule
mkdir -p index_files
mkdir -p com_files
mkdir -p data
mkdir -p summary


#creating directries inside data to store setup and analysis files
mkdir -p data/delta
mkdir -p data/species

#removing any previous files
rm -f index_files/*.ndx
rm -f com_files/*.pdb
rm -f data/species/*.csv
rm -f data/delta/*.pkl
rm -f summary/*.csv


#n_choice
n_choice=$((n_species+2))
echo "Choice for calculation of COM: $n_choice"

#defining task
task(){

       for (( j=1; j<=(($n_p_threads)); j++ ))
       do 
         k=$(( $i*$n_p_threads+$j ))
         if [ $k -le $n_total ]
       then
         echo $k &
	# using tpr file is essential otherwise .gro gives error for mass of virtual particles 
	 { echo r $k; echo q; } | gmx make_ndx -f *.tpr -o index_files/r_$k.ndx						#step 1

	#-com gives center of mass. you may use other options like -pbc, -nojump etc. as required
	 { echo $n_choice; echo "$!"; } | gmx traj -f *.xtc -s *.tpr -n index_files/r_$k.ndx -dt 20 -pbc -com -nojump -oxt com_files/com_$k.pdb	#step 2
       else
         break 1
       fi
       done

}

for (( i=0; i<=(($n_total/$n_p_threads+1)); i++ ))
do
  task "$i" &        
done

wait

echo "Center of mass calculation complete"


