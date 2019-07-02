#!/bin/bash

# module load ucsc

# Process S2 profiles
ls ../../data/S2_exp1/*.wig | while read x
do
	pathName=`dirname "$x"`
	filename=`basename "$x"`
	exp_name="${filename:0:${#filename}-4}"

	wigToBigWig $x ../../annotations/dm6.chrom.sizes ${pathName}/${exp_name}.bw
done


ls ../../data/S2_exp2/*.wig | while read x
do
	pathName=`dirname "$x"`
	filename=`basename "$x"`
	exp_name="${filename:0:${#filename}-4}"

	wigToBigWig $x ../../annotations/dm6.chrom.sizes ${pathName}/${exp_name}.bw
done


# Process Kc167 profiles
ls ../../data/Kc167_exp1/*.wig | while read x
do
	pathName=`dirname "$x"`
	filename=`basename "$x"`
	exp_name="${filename:0:${#filename}-4}"

	wigToBigWig $x ../../annotations/dm6.chrom.sizes ${pathName}/${exp_name}.bw
done
