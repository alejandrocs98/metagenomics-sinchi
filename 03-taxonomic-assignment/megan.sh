#!/bin/bash

# Fijar el directorio 06-reads-taxonomic-assignment/megan
cd ~/curso-uniandes/usuario-sinchi/06-reads-taxonomic-assignment/megan/

# Activar el ambiente de Conda de MEGAN
source ~/anaconda3/bin/activate
conda activate megan6

# path=~/curso-uniandes/usuario-sinchi/02-reads-clean/
path=../
#ext=_clean_paired_001_subsample.fastq
ext=_clean_paired_001_subsample_100k.fastq
diamond_db=~/databases/diamond/nr
megan_deb=~/databases/megan/megan-map-Feb2022-ue.db
for sample in GM1776_1_L001_R GM1776_2_L001_R; do
  out=$(basename $sample _L001_R)_R
  forward=${path}${sample}1${ext}
  reverse=${path}${sample}2${ext}
  diamond blastx -d ${diamond_db} -q $forward -o ${out}1.daa --outfmt 100
  diamond blastx -d ${diamond_db} -q $reverse -o ${out}2.daa --outfmt 100

  daa-meganizer -i ${out}1.daa -mdb $megan_deb
  daa-meganizer -i ${out}1.daa -mdb $megan_deb
done
