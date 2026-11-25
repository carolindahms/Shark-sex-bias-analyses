setwd("C:/Users/karol/Documents/Sharks/sexbias paper/Data/qpi/genes/2mb")

process_files <- function(file) {
  # Read the data
  df <- read.table(file, header = TRUE)
  type <- ifelse(grepl("autos", file), "autos", "X")
  pop_name <- sub("^([^_]+_[^_]+)_.*", "\\1", basename(file))
  df <- df %>% mutate(CHR = type, POP=pop_name)
  return(df)
}
compute_Qpi_CI <- function(df, n_mc = 10000) {
  
  A_vals <- df %>% filter(CHR == "autos") %>% pull(pi)
  X_vals <- df %>% filter(CHR == "X") %>% pull(pi)
  
  # Monte Carlo sampling
  A_star <- sample(A_vals, n_mc, replace = TRUE)
  X_star <- sample(X_vals, n_mc, replace = TRUE)
  
  Q_star <- X_star / A_star
  
  tibble(
    Api_median = median(A_star),
    Xpi_median = median(X_star),
    Qpi_median = median(Q_star),
    Qpi_mean = mean(Q_star),
    Qpi_CI_lower = quantile(Q_star, 0.025),
    Qpi_CI_upper = quantile(Q_star, 0.975),
    Api_CI_lower = quantile(A_star, 0.025),
    Api_CI_upper = quantile(A_star, 0.975),
    Xpi_CI_lower = quantile(X_star, 0.025),
    Xpi_CI_upper = quantile(X_star, 0.975)
  )
}

#-----------------------
file_list <- list.files(pattern = ".*pi\\.txt$", full.names = TRUE)

data_list <- lapply(file_list, process_files)
pi_data <- bind_rows(data_list)

pi_results <- pi_data %>%
  group_by(POP) %>%
  do(compute_Qpi_CI(.))

pi_results[c("Population", "Distance")] <- str_split_fixed(pi_results$POP, "_", 2)

pi_results$Qpi_median_CarMel<-pi_results$Qpi_median/0.938636171 # mean of da (Chagos-C.melanopterus) from 3 distance categories
pi_results$Qpi_CI_lower_CarMel<-pi_results$Qpi_CI_lower/0.938636171
pi_results$Qpi_CI_upper_CarMel<-pi_results$Qpi_CI_upper/0.938636171

#Raw Qpi
ggplot(pi_results, aes(y = Qpi_median, x = Distance, color=Population, group=Population)) +
  geom_errorbar(aes(ymin = Qpi_CI_lower, ymax = Qpi_CI_upper), width = 0,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3,alpha=0.7, position = position_dodge(width = 0.5)) +

  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey60",linewidth = 0.8) +
  geom_hline(yintercept = mean(pi_results$Qpi_median[pi_results$Population=="mal"]), 
             linetype = "dashed", color = "#ffa500ff",linewidth = 0.8) +
  
  geom_hline(yintercept = mean(pi_results$Qpi_median[pi_results$Population == "rowley"], 
                               na.rm = TRUE), 
             linetype = "dashed", 
             color = "#92b145ff", linewidth = 0.8) +
  labs(
    y = expression(Q[pi]) ,
    x = "Distance to coding sequences",
    title="Unadjusted Qpi"
  ) +
  theme_bw()+
  theme(#legend.position = "none", 
        
        panel.grid.minor  = element_blank(),
        panel.grid.major = element_blank()
  ) +
  scale_color_manual(name = "Region",values = c("#ffa500ff","#92b145ff")) 

#Mutational bias adjusted Qpi
ggplot(pi_results, aes(y = Qpi_median_CarMel, x = Distance, color=Population, group=Population)) +
  geom_errorbar(aes(ymin = Qpi_CI_lower_CarMel, ymax = Qpi_CI_upper_CarMel), width = 0,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3,alpha=0.7, position = position_dodge(width = 0.5)) +
  
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey60",linewidth = 0.8) +
  geom_hline(yintercept = mean(pi_results$Qpi_median_CarMel[pi_results$Population=="mal"]), 
             linetype = "dashed", color = "#ffa500ff",linewidth = 0.8) +
  
  geom_hline(yintercept = mean(pi_results$Qpi_median_CarMel[pi_results$Population == "rowley"], 
                               na.rm = TRUE), 
             linetype = "dashed", 
             color = "#92b145ff", linewidth = 0.8) +
  labs(
    y = expression(Q[pi]) ,
    x = "Distance to coding sequences",
    title="C.melanopterus adjusted Qpi"
  ) +
  theme_bw()+
  theme(#legend.position = "none", 
    
    panel.grid.minor  = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_color_manual(name = "Region",values = c("#ffa500ff","#92b145ff")) 




