#!/bin/bash

# We created GC-content information files through the tool and datasets provided by UCSC Genome Browser.
# Download the tool (hgGcPercent): https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
# Command used to create the hg18 GC-content information: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/gc5Base/
# Command used to create the hg19 GC-content information: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/gc5Base/

ucsc_path="${DATA_STORE}/data/ucsc"
store_path="./reference_gc"

### From http://hgdownload.cse.ucsc.edu/downloads.html 
mkdir -p "${ucsc_path}"
wget "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgGcPercent" -P "${ucsc_path}"
wget "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit" -P "${ucsc_path}"
wget "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit" -P "${ucsc_path}"
wget "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit" -P "${ucsc_path}"

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
            "${ucsc_path}/hgGcPercent" -wigOut -file="${store_path}/${ref}/${ref}_gc${resolution}_chr${chrom}.txt" -chr="chr${chrom}" -win="${resolution}" $ref "${ucsc_path}/${ref}.2bit"
        done
    done
done

chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y)
ref="mm9"

for resolution in 1000000 100000
do 
    for chrom in "${chroms[@]}"
    do
        "${ucsc_path}/hgGcPercent" -wigOut -file="${store_path}/${ref}/${ref}_gc${resolution}_chr${chrom}.txt" -chr="chr${chrom}" -win="${resolution}" $ref "${ucsc_path}/${ref}.2bit"
    done
done