#!/bin/bash -l

####	path to Folder or SAF
SAF="/media/A/SHARK/MAPGREY/RATIOS/SAF"
####	Path to folder of SFS
SFS="/media/A/SHARK/MAPGREY/RATIOS/singleSFS"
ANGSD="/home/marevo/Softwares/angsd/misc"


while read a b n1 n2; do
  $ANGSD/realSFS fst index $SAF/$a"_GL1_autos.saf.idx" $SAF/$b"_GL1_autos.saf.idx" -sfs $SFS/$a-$b"_autos.sfs" -fstout $a-$b"_autos"
  $ANGSD/realSFS fst stats $a-$b"_autos.fst.idx"
done < pop_pairs
