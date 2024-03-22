#!/bin/bash

ucsc_path="/home/jordan990301/Projects/HiCPAP/examples/data/ucsc"
store_path="data/reference_gc"

mkdir -p "${store_path}/hg18"
mkdir -p "${store_path}/hg19"
mkdir -p "${store_path}/mm9"

chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

for ref in "hg18" "hg19"
do
    for resolution in 1000000 100000
    do 
        for chrom in "${chroms[@]}"
        do
            $ucsc_path/hgGcPercent -wigOut -file="${store_path}/${ref}/${ref}_gc${resolution}_chr${chrom}.txt" -chr="chr${chrom}" -win="${resolution}" $ref "${ucsc_path}/${ref}.2bit"
        done
    done
done

chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y)
ref="mm9"

for resolution in 1000000 100000
do 
    for chrom in "${chroms[@]}"
    do
        $ucsc_path/hgGcPercent -wigOut -file="${store_path}/${ref}/${ref}_gc${resolution}_chr${chrom}.txt" -chr="chr${chrom}" -win="${resolution}" $ref "${ucsc_path}/${ref}.2bit"
    done
done