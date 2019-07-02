#!/bin/bash

# module load ucsc
# module load igvtools

export PATH=$PATH:~/IGV_2.5.3

# Process S2 profiles
ls ../../data/S2_exp1/*.bw | while read x
do
	pathName=`dirname "$x"`
	filename=`basename "$x"`
	exp_name="${filename:0:${#filename}-3}"

	bigWigToBedGraph $x ${pathName}/${exp_name}.bedGraph
	igvtools toTDF ${pathName}/${exp_name}.bedGraph ${pathName}/${exp_name}.tdf ../../annotations/dm6.chrom.sizes
	rm ${pathName}/${exp_name}.bedGraph
done


ls ../../data/S2_exp2/*.bw | while read x
do
	pathName=`dirname "$x"`
	filename=`basename "$x"`
	exp_name="${filename:0:${#filename}-3}"

	bigWigToBedGraph $x ${pathName}/${exp_name}.bedGraph
	igvtools toTDF ${pathName}/${exp_name}.bedGraph ${pathName}/${exp_name}.tdf ../../annotations/dm6.chrom.sizes
	rm ${pathName}/${exp_name}.bedGraph
done


# Process Kc167 profiles
ls ../../data/Kc167_exp1/*.bw | while read x
do
	pathName=`dirname "$x"`
	filename=`basename "$x"`
	exp_name="${filename:0:${#filename}-3}"

	bigWigToBedGraph $x ${pathName}/${exp_name}.bedGraph
	igvtools toTDF ${pathName}/${exp_name}.bedGraph ${pathName}/${exp_name}.tdf ../../annotations/dm6.chrom.sizes
	rm ${pathName}/${exp_name}.bedGraph
done

