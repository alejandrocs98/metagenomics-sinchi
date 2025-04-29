#!/bin/bash 

#SBATCH -N 1 
#SBATCH --mem=40000
#SBATCH -p short 
#SBATCH -n 12
#SBATCH --time=48:00:00
#SBATCH --output=Test_coverage_AllSamples_%A_%a.out
#SBATCH --array=1-50%20   # This will create 50 tasks numbered 1-50 and allow 20 concurrent jobs to run

ID=$( sed -n ${SLURM_ARRAY_TASK_ID}p ~/06_PhageAttack_Viruses/03_Cdhit_BlastLeoScript_Approach/03_CovCalculations/lookup_bbmap.txt )

module load bbtools/38.96 

WORKINGDIR=/hpcfs/home/ciencias_biologicas/la.chica10/06_PhageAttack_Viruses/Test_Coverage
INPUTCONTIGS=${WORKINGDIR}/Representative_Contigs_95pident_JustPhagePrediction.fasta
FASTQDIR=~/06_PhageAttack_Viruses/03_Cdhit_BlastLeoScript_Approach/03_CovCalculations/p07

# SAMPLE=Isol.F.Day.10

SAMS=${WORKINGDIR}/SAMs
STATS=${WORKINGDIR}/STATs

for STRATEGY in all toss; do
bbmap.sh ref=${INPUTCONTIGS} \
in=${FASTQDIR}/${ID}_R1.all.fastq \
in2=${FASTQDIR}/${ID}_R2.all.fastq \
nodisk \
out=${SAMS}/${STRATEGY}_${ID}.aln.sam.gz \
outu=${SAMS}/${STRATEGY}_${ID}_unmapped.aln.sam.gz \
ambiguous=${STRATEGY} \
slow=t \
physcov=t \
maxindel=100 minid=90 \
ow=t \
threads=12

pileup.sh in=${SAMS}/${STRATEGY}_${ID}.aln.sam.gz \
out=${STATS}/${STRATEGY}_${ID}.covstats \
rpkm=${STATS}/${STRATEGY}_${ID}.rpkm \
secondary=t \
ref=${INPUTCONTIGS} \
threads=12 \
32bit=t
done 