#!/bin/bash

#SBATCH --job-name=star_index # Job name
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=128000M
#SBATCH -t 05:00:00
#SBATCH --output=slurmout/star-index_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmout/star-index_%A.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}

cd ${outpath}
FASTA="/project/do2_projects/AFRI/practice/References/combined_genome.fa"

GTF="/project/do2_projects/AFRI/practice/References/combined_annotation.gtf"

mkdir star.overlap100.gencode.M35
cd star.overlap100.gencode.M35

module load star

call="STAR
    --runThreadN ${SLURM_NTASKS} \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles ${FASTA} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang 100"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
