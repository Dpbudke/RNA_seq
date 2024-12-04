#!/bin/bash

#SBATCH --job-name="htstream" # Job name
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 9
#SBATCH -t 04:00:00
#SBATCH --array=1-246
#SBATCH -o "slurmout/htstream_%A_%a.out" # File to which STDOUT will be written
#SBATCH -e "slurmout/htstream_%A_%a.err" # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dpbudke@ucdavis.edu

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

inpath="concatenated_reads"
outpath="01-HTS_Preproc"
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load htstream/1.3.3

call="hts_Stats -L ${outpath}/${sample}/${sample}.json -N 'initial stats' \
          -1 ${inpath}/${sample}/*R1_concatenated.fastq.gz \
          -2 ${inpath}/${sample}/*R2_concatenated.fastq.gz | \
      hts_SeqScreener -A ${outpath}/${sample}/${sample}.json -N 'screen phix' | \
      hts_SeqScreener -A ${outpath}/${sample}/${sample}.json -N 'count the number of rRNA reads'\
          -r -s ./References/mouse_rrna.fasta | \
      hts_SuperDeduper -A ${outpath}/${sample}/${sample}.json -N 'remove PCR duplicates' | \
      hts_AdapterTrimmer -A ${outpath}/${sample}/${sample}.json -N 'trim adapters' | \
      hts_PolyATTrim -A ${outpath}/${sample}/${sample}.json -N 'remove polyAT' | \
      hts_NTrimmer -A ${outpath}/${sample}/${sample}.json -N 'remove any remaining N characters' | \
      hts_QWindowTrim -A ${outpath}/${sample}/${sample}.json -N 'quality trim the ends of reads' | \
      hts_LengthFilter -A ${outpath}/${sample}/${sample}.json -N 'remove reads < 50bp' \
          -n -m 50 | \
      hts_Stats -A ${outpath}/${sample}/${sample}.json -N 'final stats' \
          -f ${outpath}/${sample}/${sample}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
