library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
#---------------------FUNCTIONS------------------------
process_files <- function(file) {
  # Read the data
  df <- read.table(file, header = TRUE) 
  
  # Extract the type from the filename
  type <- ifelse(grepl("auto", file), "auto", "x")
  df <- df %>% mutate(CHR = type)
  
  return(df)
}

compute_Q_CI <- function(df, n_mc = 10000) {
  
  A_vals <- df %>% filter(CHR == "auto") %>% pull(da)
  X_vals <- df %>% filter(CHR == "x") %>% pull(da)
  
  # Monte Carlo sampling
  A_star <- sample(A_vals, n_mc, replace = TRUE)
  X_star <- sample(X_vals, n_mc, replace = TRUE)
  
  Q_star <- X_star / A_star
  
  #for pi of pop1
  
  A_vals_p1 <- df %>% filter(CHR == "auto") %>% pull(pi)
  X_vals_p1 <- df %>% filter(CHR == "x") %>% pull(pi)
  
  # Monte Carlo sampling
  A_star_p1 <- sample(A_vals_p1, n_mc, replace = TRUE)
  X_star_p1 <- sample(X_vals_p1, n_mc, replace = TRUE)
  
  Q_pi <- X_star_p1 / A_star_p1
  
  
  #corrected pi/da
  pop1_qpida=(X_star_p1/X_star)/(A_star_p1/A_star)
  
  
  tibble(
    Qda_mean = mean(Q_star),
    Qda_median = median(Q_star),
    Qda_CI_lower = quantile(Q_star, 0.025),
    Qda_CI_upper = quantile(Q_star, 0.975),
    X_da_median=median(X_star),
    X_da_CI_lower=quantile(X_star, 0.025),
    X_da_CI_upper=quantile(X_star, 0.975),
    A_da_median=median(A_star),
    A_da_CI_lower=quantile(A_star, 0.025),
    A_da_CI_upper=quantile(A_star, 0.975),
    Qpi_mean = mean(Q_pi),
    Qpi_median = median(Q_pi),
    Qpi_CI_lower = quantile(Q_pi, 0.025),
    Qpi_CI_upper = quantile(Q_pi, 0.975),
    Qpi.da_median = median(pop1_qpida),
    Qpi.da_CI_lower = quantile(pop1_qpida, 0.025),
    Qpi.da_CI_upper = quantile(pop1_qpida, 0.975),
   
    
    
    
  )
}

da_long <- function(df) {
  df %>%
    #select(-matches("^(X_da|A_da)")) %>%
    
    pivot_longer(
      cols = everything(),
      names_to = "raw_metric",
      values_to = "value"
    ) %>%
    
    mutate(
      stat = case_when(
        str_detect(raw_metric, "median")   ~ "median",
        str_detect(raw_metric, "CI_lower") ~ "CI_lower",
        str_detect(raw_metric, "CI_upper") ~ "CI_upper",
        TRUE                               ~ "mean"
      ),
      
      main_stat = case_when(
        str_detect(raw_metric, "Qpi\\.da|Qpi\\.da") ~ "Qpi.da",
        str_detect(raw_metric, "Qpi|Qpi")          ~ "Qpi",
        str_detect(raw_metric, "Qda")                ~ "Qda",
        str_detect(raw_metric, "X_da")                ~ "X_da",
        str_detect(raw_metric, "A_da")                ~ "A_da",
        TRUE                                         ~ "other" 
      ),
      
      pop = case_when(
        str_detect(raw_metric, "Qpi|pop") ~ "pop1",
        str_detect(raw_metric, "Qpi|pop") ~ "pop2",
        str_detect(raw_metric, "Qda|X_da|A_da") ~ "all", # label for Qda population blanket
        TRUE                                ~ NA_character_
      )
    ) %>%
    
    filter(
      stat %in% c("median", "CI_lower", "CI_upper"), 
      main_stat %in% c("Qpi", "Qpi.da", "Qda", "X_da","A_da")       
    ) %>%
    
    select(main_stat, pop, stat, value) %>%
    pivot_wider(names_from = stat, values_from = value)
}
#------------------------------------------------------

file_list_pi <- list.files(pattern = ".*pi\\.txt$", full.names = TRUE)
file_list_da <- list.files(pattern = ".*1D\\.txt$", full.names = TRUE)

#for limbatus, cuvier and acronotus
acro<-file_list_pi[c(3,4)]
acroDa<-file_list_da[c(1,2)]

cuv<-file_list_pi[c(5,6)]
cuvDa<-file_list_da[c(5,6)]

lim<-file_list_pi[c(1,2)]
limDa<-file_list_da[c(7,8)]


data_list <- lapply(acro, process_files)
data_listda <- lapply(acroDa, process_files)

data_list <- lapply(cuv, process_files)
data_listda <- lapply(cuvDa, process_files)

data_list <- lapply(lim, process_files)
data_listda <- lapply(limDa, process_files)

#~~~~~~~~
pi_data <- bind_rows(data_list)
pi20<-pi_data[c(1:20,101:120),]

da_data<-bind_rows(data_listda)

# optional: set order inside CHR 
setorder(pi20, CHR)    # or setorder(DTpi, CHR, POS)
setorder(da_data, CHR)    # or setorder(DTda, CHR, POS)
com<-as.data.frame(cbind(pi20$pi, da_data$da, pi20$CHR))
colnames(com)<-c("pi","da","CHR")
com$pi<-as.numeric(com$pi)
com$da<-as.numeric(com$da)

da_results <- com %>%
  do(compute_Q_CI(.))


#-----outgroups----
acro_bl_f<-da_long(da_results)
write.table(acro_bl_f, file ="acronotus.key-blacktips_Qpi.da.txt", row.names = F, quote = F)

cuv_bl_f<-da_long(da_results)
write.table(cuv_bl_f, file ="cuvier-blacktips_adj_Qpi.da.txt", row.names = F, quote = F)

lim_bl_f<-da_long(da_results)
write.table(lim_bl_f, file ="limbatus-blacktips_adj_Qpi.da.txt", row.names = F, quote = F)



