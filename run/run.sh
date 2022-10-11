
## Get the MCIC scripts
git clone https://github.com/mcic-osu/mcic-scripts.git

## Put all FASTQ files in one directory
mkdir data/fastq
find data/original_fastqs -name "*fastq.gz" -exec cp -v {} data/fastq \;
#! There are 761 files in total, but we only copied 759
#! There are more R2 files than R1

find data/original_fastqs -name "*fastq.gz" | wc -l # Total number of files in new dir
ls data/fastq/*fastq.gz | wc -l  # Total number of files in new dir

for R2 in data/fastq/*_R2_*fastq.gz; do
    R1=${R2/_R2_/_R1_}
    [[ ! -f $R1 ]] && echo $R1
done
# data/fastq/HGB30_S64_L001_R1_001.fastq.gz
# data/fastq/HGB37_S53_L001_R1_001.fastq.gz
# data/fastq/HGB44_S42_L001_R1_001.fastq.gz
# data/fastq/HGB56_S91_L001_R1_001.fastq.gz
# data/fastq/HGB58_S20_L001_R1_001.fastq.gz
# data/fastq/HGB66_S21_L001_R1_001.fastq.gz
# data/fastq/HGB83_S47_L001_R1_001.fastq.gz
# data/fastq/HGB9_S2_L001_R1_001.fastq.gz


#! THESE ARE THE FILES FOR WHICH THE R1 IS MISSING
shopt -s globstar
for R2 in data/original_fastqs/**/*_R2_*fastq.gz; do
    R1=${R2/_R2_/_R1_}
    [[ ! -f $R1 ]] && echo $R1
done
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0009_L001-ds.c81acd975b674cd581cf33e0539a97ed/HGB9_S2_L001_R1_001.fastq.gz
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0030_L001-ds.8acdb5b033ea408b83eac68a99012b63/HGB30_S64_L001_R1_001.fastq.gz
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0037_L001-ds.eb72f74ac6e9434397eb755d75209cb3/HGB37_S53_L001_R1_001.fastq.gz
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0044_L001-ds.2afaec1a19864594982aa062a2b1e806/HGB44_S42_L001_R1_001.fastq.gz
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0056_L001-ds.a8e35d098e384f8a962fe1df27b57009/HGB56_S91_L001_R1_001.fastq.gz
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0058_L001-ds.c8c1f7207d9d47cda91517ecf771f9d7/HGB58_S20_L001_R1_001.fastq.gz
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0066_L001-ds.59cab63de53545c4bde717dce6aea2f0/HGB66_S21_L001_R1_001.fastq.gz
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0084_L001-ds.fc2225add5a548abbdc70d07db40ae6f/HGB83_S47_L001_R1_001.fastq.gz

#! THESE IS THE FILE FOR WHICH THE R2 IS MISSING
for R1 in data/original_fastqs/**/*_R1_*fastq.gz; do
    R2=${R1/_R1_/_R2_}
    [[ ! -f $R2 ]] && echo $R2
done
# data/original_fastqs/21-076489-293353067/FASTQ_Generation_2021-09-16_12_35_48Z-460880431/21-076489-0062_L001-ds.2d8dcf157ea04fbca2d11ea01dadb122/HGB62_S68_L001_R2_001.fastq.gz

find data/original_fastqs -name "*HGB62_*"
find data/original_fastqs -name "*HGB30_*"

for R1 in data/fastq/*_R1_*fastq.gz; do
    R2=${R1/_R1_/_R2_}
    [[ ! -f $R1 ]] && echo $R1
done

find data/original_fastqs -name "*_R2_*" | wc -l
find data/original_fastqs -name "*_R1_*" | wc -l

## Run cutadapt
for R1 in data/fastq/*_R1_*fastq.gz; do
    echo sbatch mcic-scripts/metabarcoding/cutadapt.sh \
        -i "$R1" -o results/cutadapt \
        -f GAGTGYCAGCMGCCGCGGTAA -r ACGGACTACNVGGGTWTCTAAT
done

## Infer the tree
module load R/4.0.2-gnu9.1
sbatch mcic-scripts/metabarcoding/tree_build.R -i results/dada/seqtab_V4.rds -o results/dada/tree.rds

