#!/bin/bash -l
#calculate fst from bootstraps

####	path to Folder or SAF
SAF="/media/A/SHARK/MAPGREY/RATIOS/SAF/WGS/2mb"
####	Path to folder of SFS
SFS="/media/A/SHARK/MAPGREY/RATIOS/BOOTS/WGS/2mb/2D"
ANGSD="/home/marevo/Softwares/angsd/misc"


nohup $ANGSD/realSFS $SAF/blacktip_0-50kb_GL1_autos.saf.idx $SAF/mal_0-50kb_GL1_autos.saf.idx -P 8 -fold 1 -bootstrap 10 -resample_chr 1 > $SFS/blacktip-mal_0-50kb_boots_autos.sfs &
nohup $ANGSD/realSFS $SAF/blacktip_50-100kb_GL1_autos.saf.idx $SAF/mal_50-100kb_GL1_autos.saf.idx -P 8 -fold 1 -bootstrap 10 -resample_chr 1 > $SFS/blacktip-mal_50-100kb_boots_autos.sfs &
nohup $ANGSD/realSFS $SAF/blacktip_above100kb_GL1_autos.saf.idx $SAF/mal_above100kb_GL1_autos.saf.idx -P 8 -fold 1 -bootstrap 10 -resample_chr 1 > $SFS/blacktip-mal_above100kb_boots_autos.sfs &


nohup $ANGSD/realSFS $SAF/blacktip_0-50kb_X_GL1_X.saf.idx $SAF/mal_0-50kb_females_GL1_X.saf.idx -P 2 -fold 1 -bootstrap 10 -resample_chr 1 > $SFS/blacktip-mal_0-50kb_boots_X.sfs &
nohup $ANGSD/realSFS $SAF/blacktip_50-100kb_X_GL1_X.saf.idx $SAF/mal_50-100kb_females_GL1_X.saf.idx -P 2 -fold 1 -bootstrap 10 -resample_chr 1 > $SFS/blacktip-mal_50-100kb_boots_X.sfs &
nohup $ANGSD/realSFS $SAF/blacktip_above100kb_X_GL1_X.saf.idx $SAF/mal_above100kb_females_GL1_X.saf.idx -P 2 -fold 1 -bootstrap 10 -resample_chr 1 > $SFS/blacktip-mal_above100kb_boots_X.sfs &

