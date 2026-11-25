library(stringr)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(ggpmisc)
library(reshape2)

#DArTseq
angsd<-read.table("AngsdWFst_nomal16.txt",h=T)
#RADseq
angsdL<-read.table("AngsdWFST_Lest.txt",h=T)

QFst <- angsd %>%
  group_by(pop1, pop2) %>%
  summarise(
    QFst = log(1-2*Fst_autos)/log(1-2*Fst_X)  
  ) 

# Use Reynolds Distance to include pops with Fst>0.5
QFst_RD <- angsd %>%
  group_by(pop1, pop2) %>%
  summarise(
    QFst = -log(1-Fst_autos)/-log(1-Fst_X)   # Calculate mean Fst, adjust if necessary
  ) 
#Histogram
ggplot(QFst_RD,aes(QFst))+
  geom_histogram(fill="#5f72c9")+
  geom_vline(xintercept = median(QFst_RD$QFst, na.rm = T),linetype="dashed",color="#d9ee8a",size=2)+ #median 0.62
  geom_vline(xintercept = 0.75,linetype="dashed",color="#d8514b")+
  theme_classic()+
  xlab(expression(paste("QF"["ST"] )))

angsd$Qfst<-log(1-angsd$Fst_autos)/-log(1-angsd$Fst_X)

#Make Fst matrix
QFst<-QFst_RD
QFst$pop1<-as.factor(QFst$pop1)
QFst$pop2<-as.factor(QFst$pop2)
QFst<-as.data.frame(QFst)

QFst<-QFst %>% add_row(pop1 = "mal", pop2 = "mal", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "chagos", pop2 = "chagos", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "chesterfield", pop2 = "chesterfield", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "entrecasteaux", pop2 = "entrecasteaux", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "gln", pop2 = "gln", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "grand_astrolabe", pop2 = "grand_astrolabe", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "indo", pop2 = "indo", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "matthew", pop2 = "matthew", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "ningaloo", pop2 = "ningaloo", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "north_gbr", pop2 = "north_gbr", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "scott", pop2 = "scott", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "south_gbr", pop2 = "south_gbr", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "south_nc", pop2 = "south_nc", QFst= NA)
QFst<-QFst %>% add_row(pop1 = "walpole", pop2 = "walpole", QFst= NA)

matPm1<-QFst %>% mutate_if(is.character, factor, levels=c("mal","chagos", 
                                                        "ningaloo", "rowley", "scott",
                                                        "indo","north_gbr",
                                                      "south_gbr",
                                                        "chesterfield", "entrecasteaux","gln", "grand_astrolabe",
                                                        "south_nc","walpole", "matthew")) %>% xtabs(QFst ~ pop1 + pop2, ., drop.unused.levels = F) %>% '+'(., t(.))

diag(matPm1)<-NA # diagonal NAs
melted_lespops <- reshape2::melt(matPm1) 

Our_Qfst_plot_lespops<-ggplot(data = melted_lespops, aes(pop1, pop2)) +
  
  geom_tile(aes(fill =value )) +
  scale_fill_gradientn(colours = c("#06519C" ,"#5f72c9" , "#d9ee8a","#d8514b", "#e65c4b"),
                      rescaler = ~ scales::rescale_mid(.x, mid = 0.75), limits=c(0,1.5)) + 

  ggtitle(expression(paste("QF"["ST"] )))  +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_fixed()+
  scale_y_discrete(labels=c("mal"="Maldives","chagos"="Chagos","ningaloo"="Ningaloo",
                             "rowley"="Rowley","scott"="Scott","indo"="Indo","north_gbr"="North GBR","south_gbr" ="South GBR",
                            "chesterfield"="Chesterfield", 
                            "entrecasteaux"="Enctrecasteaux","gln"="GLN",
                            "grand_astrolabe"="Grand Astrolabe",
                            "south_nc"="South NC","walpole"="Walpole","matthew"="Matthew")) +
  scale_x_discrete(labels=c("mal"="Maldives","chagos"="Chagos","ningaloo"="Ningaloo",
                            "rowley"="Rowley","scott"="Scott","indo"="Indo","north_gbr"="North GBR","south_gbr" ="South GBR",
                            "chesterfield"="Chesterfield", 
                            "entrecasteaux"="Enctrecasteaux","gln"="GLN",
                            "grand_astrolabe"="Grand Astrolabe",
                            "south_nc"="South NC","walpole"="Walpole","matthew"="Matthew"))


# Autosomal FST
merged_data<-angsd

merged_data<-merged_data %>% add_row(pop1 = "mal", pop2 = "mal", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "chagos", pop2 = "chagos", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "chesterfield", pop2 = "chesterfield", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "entrecasteaux", pop2 = "entrecasteaux", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "gln", pop2 = "gln", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "grand_astrolabe", pop2 = "grand_astrolabe", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "indo", pop2 = "indo", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "matthew", pop2 = "matthew", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "ningaloo", pop2 = "ningaloo", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "north_gbr", pop2 = "north_gbr", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "scott", pop2 = "scott", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "south_gbr", pop2 = "south_gbr", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "south_nc", pop2 = "south_nc", Fst_autos= NA, Fst_X= NA)
merged_data<-merged_data %>% add_row(pop1 = "walpole", pop2 = "walpole", Fst_autos= NA, Fst_X= NA)


matAUTO<-merged_data %>% mutate_if(is.character, factor, levels=c( "mal","chagos",
                                                           "ningaloo", "rowley","scott",
                                                          "indo",
                                                          "south_gbr","north_gbr",
                                                          "chesterfield", "entrecasteaux","gln","grand_astrolabe",
                                                          "south_nc","walpole", "matthew")) %>% xtabs(Fst_autos ~ pop1 + pop2, ., drop.unused.levels = F) %>% '+'(., t(.))





diag(matAUTO)<-NA # diagonal is missing (no data)
melted_AUTO <- reshape2::melt(matAUTO) # this "melts" the table into a edge list
######
# X

matX<-merged_data %>% mutate_if(is.character, factor, levels=c( "mal","chagos",
                                                                 "ningaloo", "rowley","scott",
                                                                 "indo",
                                                                 "south_gbr","north_gbr",
                                                                 "chesterfield", "entrecasteaux","gln","grand_astrolabe",
                                                                 "south_nc","walpole", "matthew")) %>% xtabs(Fst_X ~ pop1 + pop2, ., drop.unused.levels = F) %>% '+'(., t(.))


diag(matX)<-NA # 
melted_X <- reshape2::melt(matX)

# Combine X and AUTO
combined_matrix <- matrix(NA, nrow = 15, ncol = 15)

# Assign the upper triangle
combined_matrix[upper.tri(combined_matrix)] <- matAUTO[upper.tri(matAUTO)]

# Assign the lower triangle
combined_matrix[lower.tri(combined_matrix, diag = TRUE)] <- matX[lower.tri(matX, diag = TRUE)]

row.names(combined_matrix)<-rownames(matAUTO)
colnames(combined_matrix)<-colnames(matAUTO)
diag(combined_matrix)<-NA 

melted_combo <- reshape2::melt(combined_matrix)

pdf("Autos_X_Fst.pdf", width=15,height = 12)
our_fst<-ggplot(data = melted_combo, aes(Var1, Var2)) +
  
  geom_tile(aes(fill =value )) +
  scale_fill_gradientn(colours = c("#06519C" ,"#5f72c9" , "#d9ee8a","#d8514b", "#e65c4b"),
                       limits=c(0,0.7), na.value = "grey",
                       rescaler = ~ scales::rescale_mid(.x, mid = 0.2)) + 
  
  
  ggtitle(expression(paste("F"["ST"] )))  +
  
  xlab("X Chromosome") + ylab("Autosomes") +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  
  coord_fixed()+
  scale_y_discrete(labels=c("chagos"="Chagos","chesterfield"="Chesterfield", "mal"="Maldives",
                            "entrecasteaux"="Enctrecasteaux","gln"="GLN",
                            "grand_astrolabe"="Grand Astrolabe","indo"="Indo",
                            "matthew"="Matthew","ningaloo"="Ningaloo",
                            "south_gbr" ="South GBR",
                            "rowley"="Rowley","scott"="Scott","souht-gbr"="South GBR","north_gbr"="North GBR",
                            "south_nc"="South NC","walpole"="Walpole")) +
  scale_x_discrete(labels=c("chagos"="Chagos","chesterfield"="Chesterfield", "mal"="Maldives",
                            "entrecasteaux"="Enctrecasteaux","gln"="GLN",
                            "grand_astrolabe"="Grand Astrolabe","indo"="Indo",
                            "matthew"="Matthew","ningaloo"="Ningaloo",
                            "south_gbr" ="South GBR","north_gbr"="North GBR",
                            "rowley"="Rowley","scott"="Scott","souht-gbr"="South GBR",
                            "south_nc"="South NC","walpole"="Walpole"))

dev.off

