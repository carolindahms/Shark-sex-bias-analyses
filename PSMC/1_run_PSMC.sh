### Specify bam file path and reference
reference="/media/B/SHARK_GENOME/HiC/mappedPT_yahs_scaffolds_final.fa"
PREFIX=$1 #pop name
bam=$2 #bam path

### Obtain consensus sequence from bam file
mkdir ./$PREFIX
cd ./$PREFIX
bcftools mpileup -C 50 -Ou -f  $reference $bam | bcftools call -c - | /home/marevo/Softwares/bcftools-1.16/misc/vcfutils.pl vcf2fq -d 10 -D 60 | gzip > $PREFIX".fq.gz"
echo "Consensus sequence conversion complete: $(date)"

/home/marevo/Softwares/psmc/utils/fq2psmcfa -q20 $PREFIX".fq.gz" > $PREFIX".psmcfa"
/home/marevo/Softwares/psmc/utils/splitfa $PREFIX".psmcfa" > $PREFIX".split.psmcfa"

seq 100 | parallel -j 10 -I% --max-args 1 /home/marevo/Softwares/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-%.psmc $PREFIX".split.psmcfa"  &> boot.out 

#test if better when splitting first time window to "2+2" instead of 4
#seq 10 | parallel -j 10 -I% --max-args 1 /home/marevo/Softwares/psmc/psmc -N25 -t15 -r5 -b -p "2+2+25*2+4+6" -o round-%_2+2.psmc $PREFIX".split.psmcfa"  &> boot2.out 
#seq 1 | parallel -j 1 -I% --max-args 1 /home/marevo/Softwares/psmc/psmc -N25 -t15 -r5 -b -p "1+1+1+1+25*2+4+6" -o round-%_+1.psmc $PREFIX".split.psmcfa"  &> boot1.out 

cat $PREFIX".psmc" round-*.psmc > $PREFIX".combined.psmc"
#cat $PREFIX".psmc" round-*_2+2.psmc > $PREFIX"_2+2.combined.psmc"

#rm round-*.psmc

#generation time robbins 2006, mut rate from Lesturgie
/home/marevo/Softwares/psmc/utils/psmc_plot.pl -u 1.9e-8 -g 16.4 -R $PREFIX".100bootplot_mu=1.9e-8.psmc" $PREFIX".combined.psmc"
