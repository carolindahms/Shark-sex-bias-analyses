#!/bin/bash -l
SAF=/media/A/SHARK/MAPGREY/RATIOS/SAF # path for allele frequencies
SFS=/media/A/SHARK/MAPGREY/RATIOS/SFS # path for SFS
ANGSD="/home/marevo/Softwares/angsd/misc"
FILTERS=" -uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 30 -minQ 20" # filters as per Laturgie paper

a="$1" # pop1
b="$2" # pop2

$ANGSD/realSFS $SAF/${a}_GL1_autos.saf.idx $SAF/${b}_GL1_autos.saf.idx -P 4 -fold 1 -bootstrap 20 -resample_chr 1  > $SFS/${a}-${b}_autos_boots.sfs
