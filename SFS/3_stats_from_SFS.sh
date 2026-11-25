#!/bin/bash -l
#calculate fst from bootstraps

#### Bamlist prefix , assume your Bamlists are called "Bamlist_POP"
BAMLIST="/media/A/SHARK/MAPGREY/BAMLISTS_FILT"
####	path to Genome
REF="/media/A/SHARK/REF/sCarAmb.v1.hapA.fa"
####	path to Folder or SAF
SAF="/media/A/SHARK/MAPGREY/RATIOS/SAF"
####	Path to folder of SFS
SFS="/media/A/SHARK/MAPGREY/RATIOS/singleSFS"
ANGSD="/home/marevo/Softwares/angsd/misc"
FILTERS=" -uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 30 -minQ 20 -doHWE 1 -hetbias_pval 1e-4 -hwe_pval 1e-4" # filters as per Laturgie paper

pop_pairs="$1" #tab separated list of pop1name pop2name popsize1 popsize2

#FOR AUTOS

while read a b n1 n2; do
    # Modify the awk line to add 1 to n1 and n2 directly in the print statement
    awk -v n1="$n1" -v n2="$n2" '{print (2*n1) + 1, (2*n2) + 1, "unfolded"; print}' $SFS/$a-$b"_autos.sfs" > $SFS/$a-$b"_autos_header.sfs"
    python allStats_moments_1Dsfs_boots.py  $SFS/$a-$b"_autos_header.sfs" "${a}" "${b}_autos" "${n1}" "${n2}"
done < $pop_pairs


#FOR X
while read a b n1 n2; do
    # Modify the awk line to add 1 to n1 and n2 directly in the print statement
    awk -v n1="$n1" -v n2="$n2" '{print (2*n1) + 1, (2*n2) + 1, "unfolded"; print}' $SFS/$a-$b"_X.sfs" > $SFS/$a-$b"_X_header.sfs"
    python allStats_moments_1Dsfs_boots.py $SFS/$a-$b"_X_header.sfs" "${a}" "${b}_X" "${n1}" "${n2}"

#done < $pop_pairs_f

