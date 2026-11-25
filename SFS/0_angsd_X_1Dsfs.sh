#calculate fst from bootstraps

#### Bamlist prefix , assume your Bamlists are called "Bamlis_POP"
BAMLIST="/media/A/SHARK/MAPGREY/BAM/blacktips"
####	path to Genome
REF="/media/B/SHARK_GENOME/HiC/mappedPT_yahs_scaffolds_final.fa"
####	path to Folder or SAF
SAF="/media/A/SHARK/MAPGREY/RATIOS/SAF/WGS/2mb"
####	Path to folder of SFS
SFS="/media/A/SHARK/MAPGREY/RATIOS/SFS/WGS"
ANGSD="/home/marevo/Softwares/angsd/misc"
FILTERS=" -uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 30 -minQ 20" # filters as per Laturgie paper

#SAF generated from filtered bams (without TE and cds) /media/A/SHARK/MAPGREY/BAMLISTS_FILT/

for i in  $BAMLIST/blacktip_0-50kb_X $BAMLIST/blacktip_50-100kb_X $BAMLIST/blacktip_above100kb_X; do
   a=$(cat $i | wc -l) # get number of individuals in each pop
   b=$(echo $a*0.25 | bc | awk '{print int($1+0.5)}') # equal to 50% of individuals (to be used as the minimum individuals seuqenced at depth >1 as filter to retain site)
	echo "retain sites sequenced in $b individuals"
	angsd -P 4 -bam $i -doSaf 1 -doMajorMinor 1 -doCounts 1 -doMaf 1 $FILTERS -minInd $b -r scaffold_36 -anc $REF -GL 1 -ref $REF -out $SAF/$(basename $i)_GL1_X 
done 

#SFS without PAR on X
for i in blacktip_0-50kb_X blacktip_50-100kb_X blacktip_above100kb_X; do
 $ANGSD/realSFS $SAF/${i}_GL1_X.saf.idx -P 4 -fold 1 -r scaffold_36:4700000-17341187 -bootstrap 100 > $BOOT/${i}_female_X_noPAR_boots.sfs
done 

