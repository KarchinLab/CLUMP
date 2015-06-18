#!/bin/bash

#USAGE bash runlargescaleclump inputfile controlfile numberofpermutations

inputfile=$1
controlfile=$2
numberofpermutations=$3


awk '{print $2}' $inputfile | sort | uniq > inputfile.genes.txt

mkdir -p input_genes
mkdir -p control_genes
mkdir -p clump_results

for gene in `cat inputfile.genes.txt`
do
grep -w $gene $inputfile > input_genes/$gene
grep -w $gene $controlfile > control_genes/$gene
python combined.clump.py -f input_genes/$gene -p protein.2.length.txt -c control_genes/$gene -z $numberofpermutations > clump_results/$gene
##GENES without enough mutations will appear as a blank text file in the output directory 
done

cat clump_results/* > clump.permutations.testing.txt
