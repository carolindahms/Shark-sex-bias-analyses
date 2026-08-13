setwd("C:/Users/karol/Documents/Sharks/dA")
library(dplyr)
library(tidyr)
#---------------------FUNCTIONS------------------------
process_files <- function(file) {
  # Read the data
  df <- read.table(file, header = TRUE)  # Adjust if the header is absent
  
  # Extract the type from the filename
  type <- ifelse(grepl("auto", file), "auto", "x")
  # Add the new column
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
  
  A_vals_p1 <- df %>% filter(CHR == "auto") %>% pull(pi1)
  X_vals_p1 <- df %>% filter(CHR == "x") %>% pull(pi1)
  
  # Monte Carlo sampling
  A_star_p1 <- sample(A_vals_p1, n_mc, replace = TRUE)
  X_star_p1 <- sample(X_vals_p1, n_mc, replace = TRUE)
  
  Q_pi1 <- X_star_p1 / A_star_p1
  
  #for pi of pop2
  
  A_vals_p2 <- df %>% filter(CHR == "auto") %>% pull(pi2)
  X_vals_p2 <- df %>% filter(CHR == "x") %>% pull(pi2)
  
  # Monte Carlo sampling
  A_star_p2 <- sample(A_vals_p2, n_mc, replace = TRUE)
  X_star_p2 <- sample(X_vals_p2, n_mc, replace = TRUE)
  
  Q_pi2 <- X_star_p2 / A_star_p2
  
  #corrected pi/da
  pop1_qpida=(X_star_p1/X_star)/(A_star_p1/A_star)
  pop2_qpida=(X_star_p2/X_star)/(A_star_p2/A_star)
  
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
    Qpi1_mean = mean(Q_pi1),
    Qpi1_median = median(Q_pi1),
    Qpi1_CI_lower = quantile(Q_pi1, 0.025),
    Qpi1_CI_upper = quantile(Q_pi1, 0.975),
    Qpi2_mean = mean(Q_pi2),
    Qpi2_median = median(Q_pi2),
    Qpi2_CI_lower = quantile(Q_pi2, 0.025),
    Qpi2_CI_upper = quantile(Q_pi2, 0.975),
    Qpi1.da_median = median(pop1_qpida),
    Qpi1.da_CI_lower = quantile(pop1_qpida, 0.025),
    Qpi1.da_CI_upper = quantile(pop1_qpida, 0.975),
    Qpi2.da_median = median(pop2_qpida),
    Qpi2.da_CI_lower = quantile(pop2_qpida, 0.025),
    Qpi2.da_CI_upper = quantile(pop2_qpida, 0.975)
    
    
    
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
        str_detect(raw_metric, "Qpi1\\.da|Qpi2\\.da") ~ "Qpi.da",
        str_detect(raw_metric, "Qpi1|Qpi2")          ~ "Qpi",
        str_detect(raw_metric, "Qda")                ~ "Qda",
        str_detect(raw_metric, "X_da")                ~ "X_da",
        str_detect(raw_metric, "A_da")                ~ "A_da",
        TRUE                                         ~ "other" 
      ),
      
      pop = case_when(
        str_detect(raw_metric, "Qpi1|pop1") ~ "pop1",
        str_detect(raw_metric, "Qpi2|pop2") ~ "pop2",
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

file_list_o <- list.files(pattern = ".*1D\\.txt$", full.names = TRUE)
bl<-file_list_o[1:2]
ind<-file_list_o[3:4]
sv<-file_list_o[5:6]

data_list <- lapply(sv, process_files)
data_list <- lapply(bl, process_files)
data_list <- lapply(ind, process_files)

da_data <- bind_rows(data_list)

da_results <- da_data %>%
  do(compute_Q_CI(.))


chagos_indo<-da_long(da_results)
chagos_indo$pop<-gsub("pop1","Chagos(Indo)", chagos_indo$pop)
chagos_indo$pop<-gsub("pop2","Indo", chagos_indo$pop)
chagos_indo$pop<-gsub("all","Chagos-Indo", chagos_indo$pop)

write.table(chagos_indo, file ="chagos-indo_Qpi.da.txt", row.names = F, quote = F)
write.table(chagos_bl, file ="chagos-blacktips_adj_Qpi.da.txt", row.names = F, quote = F)
write.table(chagos_sv, file ="chagos-silvertips_adj_Qpi.da.txt", row.names = F, quote = F)

chagos_bl<-da_long(da_results)
chagos_bl$pop<-gsub("pop1","Chagos(Blacktips)", chagos_bl$pop)
chagos_bl$pop<-gsub("pop2","Blacktips", chagos_bl$pop)
chagos_bl$pop<-gsub("all","Chagos-Blacktips", chagos_bl$pop)

chagos_sv<-da_long(da_results)
chagos_sv$pop<-gsub("pop1","Chagos(Silvertips)", chagos_sv$pop)
chagos_sv$pop<-gsub("pop2","Silvertips", chagos_sv$pop)
chagos_sv$pop<-gsub("all","Chagos-Silvertips", chagos_sv$pop)

bind<-rbind(chagos_bl,chagos_indo, chagos_sv)
pldf<-bind
p<-ggplot(pldf[pldf$main_stat=="Qpi"|pldf$main_stat=="Qpi.da",], aes(y = median, x = pop,color=main_stat)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, position = position_dodge(width=0.4)) +
  geom_point(size = 3, position=position_dodge(width=0.4)) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey60") +
  theme(legend.position = "inside", legend.position.inside = c(0.86,0.82),legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))+
  theme_bw() 

p+ylim(0,0.0015)


