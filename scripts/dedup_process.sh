#!/bin/bash
#SBATCH --job-name=aptaseq_analysis
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=100G

source activate python3
unset PYTHONPATH

base_dir="/home/FCAM/slawson/miseq/210610_M01212_0292_000000000-JPPMK/Data/Intensities/BaseCalls"

samples="10k-5_S1 10k-50_S2 10k-100_S3 50k-5_S4 50k-50_S5 50k-100_S6 100k-5_S7 100k-50_S8"

for sample in $samples; do
    zcat $base_dir/$sample"_L001_R1_001.fastq.gz" > $base_dir/$sample"_L001_R1_001.fastq"
    python /home/FCAM/slawson/UMI-Reducer/AmpUMI/AmpUMI.py Process --fastq $base_dir/$sample"_L001_R1_001.fastq" --fastq_out $base_dir/AptaSeq/deduped/$sample"_R2_processed.fastq" --umi_regex "^IIIIIIIIII" > out.txt
    python $base_dir/AptaSeq/group_analysis.py $sample
done

python $base_dir/AptaSeq/final_stats.py
rm out.txt
rm groups.txt
rm skews.txt
rm controls.txt
rm unique.txt
