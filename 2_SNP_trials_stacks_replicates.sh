#!/bin/bash

# This script was used to run stacks for M values 1-4 and corresponding n values to assess SNP error

# set working directory & popmap directory
work=~/Gentianella/whole_data/SRA_Sequence/trimmed_100BP/
popmap=~/Gentianella/whole_data/Stacks_Trials/popmap_reps.txt

# For Loop M values
for M in 1 2 3 4; do
    # Calculate the value of n based on the value of M
    if [ ${M} -eq 1 ]; then
        n_high=$((${M} + 1))
        n_mid=$((${M}))
        echo "running stacks for $M and $n_high and $n_mid"
# Create new directory
    out=~/Gentianella/whole_data/Stacks_Trials/denovo_reps.M${M}.n${n_high}.m3
    mkdir $out
    # Move into the new directory
    cd $out
        # Stacks command
        denovo_map.pl --samples $work --popmap $popmap \
            --out-path $out \
            --paired \
            -M $M \
            -n $n_high \
            -T 6 \
            --min-samples-per-pop 0.8 \
            -X "ustacks: --force-diff-len" \
            -X "ustacks: -m 3"
    # Create new directory
    out=~/Gentianella/whole_data/Stacks_Trials/denovo_reps.M${M}.n${n_mid}.m3
    mkdir $out
    # Move into the new directory
    cd $out
        # Stacks command
        denovo_map.pl --samples $work --popmap $popmap \
            --out-path $out \
            --paired \
            -M $M \
            -n $n_mid \
            -T 6 \
            --min-samples-per-pop 0.8 \
            -X "ustacks: --force-diff-len" \
            -X "ustacks: -m 3"
    else
        n_low=$((${M} - 1))
        echo "running stack for $M and $n_low"
        # Create new directory
    out=~/Gentianella/whole_data/Stacks_Trials/denovo_reps.M${M}.n${n_low}.m3
    mkdir $out
    # Move into the new directory
    cd $out
        # Stacks command
        denovo_map.pl --samples $work --popmap $popmap \
            --out-path $out \
            --paired \
            -M $M \
            -n $n_low \
            -T 6 \
            --min-samples-per-pop 0.8 \
            -X "ustacks: --force-diff-len" \
            -X "ustacks: -m 3"
	n_high=$((${M} + 1))
    echo "running stacks for $M and $n_high"
	# create new directory
    out=~/Gentianella/whole_data/Stacks_Trials/denovo_reps.M${M}.n${n_high}.m3
    mkdir $out
    # Move into the new directory
    cd $out
        # Stacks command
        denovo_map.pl --samples $work --popmap $popmap \
            --out-path $out \
            --paired \
            -M $M \
            -n $n_high \
            -T 6 \
            --min-samples-per-pop 0.8 \
            -X "ustacks: --force-diff-len" \
            -X "ustacks: -m 3"
    n_mid=$((${M}))
    echo "running stacks for $M and $n_mid"
	# create new directory
    out=~/Gentianella/whole_data/Stacks_Trials/denovo_reps.M${M}.n${n_mid}.m3
    mkdir $out
    # Move into the new directory
    cd $out
        # Stacks command
        denovo_map.pl --samples $work --popmap $popmap \
            --out-path $out \
            --paired \
            -M $M \
            -n $n_mid \
            -T 6 \
            --min-samples-per-pop 0.8 \
            -X "ustacks: --force-diff-len" \
            -X "ustacks: -m 3"
    fi
done

#Calculate the r80 loci
data=~/Gentianella/whole_data/Stacks_Trials/ ##This is the location of the folders

for d in ${data}denovo_reps*
    do
echo "moving into $d"
cd $d #move into directory
echo "now in $d"
populations -P $d -M $data/popmap_reps.txt -r 80 --vcf
loc=$(cat populations.sumstats.tsv | #go into the output file
grep -v '^#' | #remove comment line
cut -f 1 | # Cut first column (Locus ID)
sort -n -u | wc -l) # Sort numerically and unique values
echo "There are ${loc} r80 loci" #print the number of r80 loci
cd $data
echo "now in $data"
echo "${loc}" >>r80_loci_reps.txt
done
