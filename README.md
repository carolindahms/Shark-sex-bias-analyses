# Shark-sex-bias-analyses
Scripts and pipelines for the shark sex bias paper

1. In angsd get SAF and generate SFS
2. Get weighted Fst with realSFS (5_angsd_fst.sh)
4. Calculate nucleotide diversity with pi.R from 1D SFS
5. Calculated dA from 2D SFS (with allStats_moments_1Dsfs_params.py)
6. Adjust Qpi for mutational bias and plot (1_pi_Qpi_plots.R)
7. Get Qpi for outgroup species (2_Qpi_adj_OUTGROUPS.R)
8. Test Qpi across distance classes from coding regions (3_Qpi_genes.R)
9. Plot Fsts and bootstrap QFst with pop replacement
   
