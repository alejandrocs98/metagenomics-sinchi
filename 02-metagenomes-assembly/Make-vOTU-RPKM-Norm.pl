#!/usr/bin/perl -w

# Written by Alejandro Reyes
# El programa solo necesita llamarlo desde la carpeta donde están todos los cov_stats y el único parametro que lleva es el nombre del archivo de salida.
# Todos los archivos deben seguir el siguiente esquema: toss_SAMPLEID.covstats.txt y all_SAMPLEID.covstats.txt
# Donde SAMPLEID es el nombre de cada muestra.
# Por último asume que son PE 2x150, si eso cambia toca cambiar el código. Puedo poner una variable que el usuario dé, porque no hay forma de detectar eso desde los covstats.

use strict;
use List::Util qw(max min);

if (@ARGV != 1) {
	die "\nUsage: Make-vOTU-RPKM-Norm.pl <outfile>\n\n";
}

my $outfile = shift;


my @coverage_files = `ls toss_*.covstats.txt`;

open(OUT, ">$outfile") or die ("Couldn't open outfile $outfile\n");
print OUT "Contig\tSample\tRPKM\tCovRatio\tNeededRatio\tExpectedRatio\n";

foreach my $f (@coverage_files){
	print "Processing File $f";
	&process_cov($f);
}close OUT;





sub process_cov{
  my $file = shift(@_);

my $all=$file;
$all=~ s/toss/all/;
my $sample= $file;
chop($sample);
$sample=~ s/toss_//;
$sample=~ s/.covstats.txt//;



my $CovTossFile = $file;
open (IN, "<$CovTossFile") or die ("Couldn't open file: $CovTossFile\n");

my $CovAllFile = $all;
open (COVALL, "<$CovAllFile") or die ("Couldn't open file: $CovAllFile\n");

my $TotalReadsMap = 0;
my %contig_len=();
my %extra_reads=();
my %needed_reads=();
my %toss_cov=();
my %deltaMapBase=();
my %contigs_all=();
my %contigs_toss=();
my %extra=();
my %CovRatio=();
my %NeededRatio=();
my %ExpRatio=();

while (my $line = <IN>){
  chomp $line;
	next unless ($line=~/contig/);
  my @temp = split /\s+/, $line;
	my $name=shift(@temp);
  @{$contigs_toss{$name}}=@temp;
	$TotalReadsMap+=$contigs_toss{$name}[6];
	$TotalReadsMap+=$contigs_toss{$name}[5];
}
close IN;

while (my $line = <COVALL>){
  chomp $line;
	next unless ($line=~/contig/);
  my @temp = split /\s+/, $line;
	my $name=shift(@temp);
	die ("Name $name does not exist in Toss\n") unless $contigs_toss{$name};
  @{$contigs_all{$name}}=@temp;
	$contig_len{$name}=$contigs_all{$name}[1];
	#print "Estas son las reads procesando: $contigs_all{$name}[6] + $contigs_all{$name}[5] - $contigs_toss{$name}[6] + $contigs_toss{$name}[5])\n";
	$extra_reads{$name}=abs(($contigs_all{$name}[6]+$contigs_all{$name}[5])-($contigs_toss{$name}[6]+$contigs_toss{$name}[5]));
	$toss_cov{$name}= ($contigs_toss{$name}[4] > 0) ? ($contigs_toss{$name}[5]+$contigs_toss{$name}[6])/$contigs_toss{$name}[4] : 0;
	$deltaMapBase{$name}=abs($contigs_all{$name}[4]-$contigs_toss{$name}[4]);
	#print "Para Needed reads se necesita la variacion en el mapeo que es $contigs_all{$name}[4] - $contigs_toss{$name}[4]\n";
	$needed_reads{$name}= sprintf "%.0f", ($toss_cov{$name}*$deltaMapBase{$name}*0.9);
	#print "Extra es el min entre $extra_reads{$name} y $needed_reads{$name}\n";
	$extra{$name}=min(abs($extra_reads{$name}),abs($needed_reads{$name}));
	$TotalReadsMap+=$extra{$name};
	$CovRatio{$name}=log(($contigs_all{$name}[3]+1)/($contigs_toss{$name}[3]+1))/log(10);
	#print "Needed reads es $needed_reads{$name}, Extra reads es: $extra_reads{$name}\n";
	$NeededRatio{$name}=log(($needed_reads{$name}+1)/($extra_reads{$name}+1))/log(10);
	$ExpRatio{$name}=($contigs_all{$name}[3]>0 && ($contigs_toss{$name}[6]+$contigs_toss{$name}[5])>0) ? log((100*(1-exp(-1*(($contigs_toss{$name}[6]+$contigs_toss{$name}[5]+$needed_reads{$name})*150)/$contig_len{$name})))/$contigs_all{$name}[3])/log(10) : 0;
}
close COVALL;


# Ahora si va a hacer la nueva tabla con los resultados por contig.
# RPKM es Reads por 1Kb de contig por 1M de reads mapeadas.
# Los reads son los reads mapeados en Toss + el min (reads extras mapeada en All; reads necesitadas para igual coverage de Toss)


my $RPKM=0;

#print "Las reads totales son $TotalReadsMap\n";


foreach my $k (keys %contigs_toss){
	$RPKM=($contigs_toss{$k}[5]+$contigs_toss{$k}[6]+$extra{$k})*1E9/($contig_len{$k}*$TotalReadsMap);
	print OUT "$k\t$sample\t$RPKM\t$CovRatio{$k}\t$NeededRatio{$k}\t$ExpRatio{$k}\n";
}
}

