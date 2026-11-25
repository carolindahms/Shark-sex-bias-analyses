#!/bin/bash -l
#calculate fst from bootstraps

####	path to Folder or SAF
SAF="/media/A/SHARK/MAPGREY/RATIOS/SAF"
####	Path to folder of SFS
SFS="/media/A/SHARK/MAPGREY/RATIOS/singleSFS"
ANGSD="/home/marevo/Softwares/angsd/misc"


while read a b n1 n2; do
 $ANGSD/realSFS $SAF/$a"_females_GL1_X.saf.idx" $SAF/$b"_females_GL1_X.saf.idx" -r scaffold_36:4700000-17341187 -P 8 -fold 1 > $SFS/$a-$b"_X.sfs"
done < pop_pairs
