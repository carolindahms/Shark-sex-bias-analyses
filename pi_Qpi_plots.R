library(dplyr)
library(tidyr)

##-------functions------

process_files <- function(file) {
  # Read the data
  df <- read.table(file, header = TRUE)  # Adjust if the header is absent
  
  # Extract the type from the filename
  type <- ifelse(grepl("autos", file), "autos", "X")
  pop_name <- sub("_(.*)", "", basename(file))
  # Add the new column
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

new_levels <- c("C.acronotus.GA","C.limbatus","G.cuvier","Juan", "Zelee", "mal", "chagos", "ningaloo", "rowley", "scott", 
                "indo", "north-gbr", "south-gbr", 
                "Bampton","chesterfield", "gln", "entrecasteaux", 
                "grand-astrolabe", "petit-astrolabe", 
                "south-nc", "walpole", "matthew","Niku","Palmyra","Fakarava")

pi_results$POP <- as.factor(pi_results$POP)
pi_results$POP <- factor(pi_results$POP, levels = new_levels)

ocean_o2<-c("WPA","Other species", "Other species","CPA","Other species", 
            "WIO", "WPA","CPA","WIO",'CIO',
            "WPA","WPA","WPA","WPA", "WPA",
            "CIO","WPA", "EIO","WPA","WPA",
            "EIO","EIO", "WPA","WPA","WPA")

pi_results$Region<-ocean_o2

pi_results$lon<-c(158.14727, 32, 34, 214.4499, 36,
                  42.58667 ,185.481, 197.9219,46.11179,72.33611 ,
                  158.61110 ,162.92411,163.9241082, 165.588336, 129.783225,
                  73.18841, 171.361578, 113.62166, 145.459672, 166.55546,
                  119.349853, 121.942425, 151.933803, 166.452488, 168.9946959)

#################################################################################
## mutational bias adjusted Qpi based on proxy (da)
## between Chagos vs Indonesia and outgroups blacktips and silvertips
#chagos-blacktip Xda/Ada ratio = 0.7744253 , multiply with raw Qpi for adjustment
#chagos-silvertip Xda/Ada = 0.7028057
#chagos-Indo Xda/Ada = 0.8264699

pi_results$Qpi_median_CarMel<-pi_results$Qpi_median/0.7744253
pi_results$Qpi_median_CarAlb<-pi_results$Qpi_median/0.7028057
pi_results$Qpi_median_Indo<-pi_results$Qpi_median/0.8264699

pi_results$Qpi_CI_lower_CarMel<-pi_results$Qpi_CI_lower/0.7744253 #Xci_lower/Aci_lower
pi_results$Qpi_CI_lower_CarAlb<-pi_results$Qpi_CI_lower/0.7028057
pi_results$Qpi_CI_lower_Indo<-pi_results$Qpi_CI_lower/0.8264699

pi_results$Qpi_CI_upper_CarMel<-pi_results$Qpi_CI_upper/0.7744253 #Xci_up/Aci_up
pi_results$Qpi_CI_upper_CarAlb<-pi_results$Qpi_CI_upper/0.7028057
pi_results$Qpi_CI_upper_Indo<-pi_results$Qpi_CI_upper/0.8264699


#for outgroups use only C. melanopterus adj values
#Acronotus GA 
pi_results[2,14]<-0.62600606 #median
pi_results[2,17]<-0.39767682 #lower
pi_results[2,20]<-0.81366433 #upper

#limbatus  
pi_results[3,14]<-0.37867063
pi_results[3,17]<-0.31929675
pi_results[3,20]<-0.42679149

#cuvier  
pi_results[5,14]<-0.24010937 
pi_results[5,17]<-0.17524379  
pi_results[5,20]<-0.31596299

pi_ordered <- pi_results[order(as.integer(pi_results$POP)), ]
write.table(pi_ordered, "allQpi_daAdjusted_lon.txt", quote=F, row.names = F)

#########################
## QPi across longitude

qpi<-ggplot(pi_results, aes(y = Qpi_median_CarMel, x = as.numeric(lon),color=Region, label=POP)) +
  geom_errorbar(aes(ymin = Qpi_CI_lower_CarMel, ymax = Qpi_CI_upper_CarMel), width = 0) +
  geom_point(size = 3,alpha=0.7) +
  ylim(0.1,1 )+
  # geom_text_repel(angle=90,nudge_y = -0.2, vjust = 1,  box.padding = 0.5,   hjust=1)+
  
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey60") +
  
  labs(
   # y = expression(Q[pi] == pi[X] / pi[A]),
    y = expression(Q[pi] /dA) ,
    x = ""
  ) +
  theme_bw()+
  scale_x_continuous(breaks=c(60, 90, 120, 150, 180, 210))+
  theme(legend.position = "none", 
        #legend.direction = "horizontal",
        #legend.position.inside = c(0.73,0.85),
        #legend.text = element_text(size = 8),
        #legend.title = element_text(size = 9),
        #legend.background = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.grid.major = element_blank()
        ) +
  scale_color_manual(name = "Region",values = c("#ffa500ff","#3a2b59ff","#92b145ff",
                                                "grey","#B21807","#8a85edff")) 







#compare with uncorrected Qpis:
library(data.table)
DT <- as.data.table(pi_results)

no <- DT[, .(POP, lon, Region, Qpi_median, Qpi_CI_lower, Qpi_CI_upper)]
no[, correction := "No"]

yes <- DT[, .(POP, lon,Region,
              Qpi_median = Qpi_median_CarMel,
              Qpi_CI_lower = Qpi_CI_lower_CarMel,
              Qpi_CI_upper = Qpi_CI_upper_CarMel)]
yes[, correction := "Yes"]

outDT <- rbindlist(list(no, yes), use.names = TRUE)
setcolorder(outDT, c("POP","lon","Region","Qpi_median","Qpi_CI_lower","Qpi_CI_upper","correction"))

nudge_amt <- 1

ggplot(outDT, aes(x = as.numeric(lon), y = Qpi_median, color = Region, shape=correction)) +
  geom_errorbar(aes(ymin = Qpi_CI_lower, ymax = Qpi_CI_upper, group = correction),
                width = 0,
                position = position_nudge(x = ifelse(outDT$correction == "Yes", nudge_amt, -nudge_amt))) +
  geom_point(aes(group = correction),
             size = 3,
             position = position_nudge(x = ifelse(outDT$correction == "Yes", nudge_amt, -nudge_amt))) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey60") +
  labs(y = expression(Q[pi] == pi[X] / pi[A]), x = "") +
  theme_bw() +
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180, 210)) +
  theme(legend.direction = "horizontal",
        legend.position = c(0.73, 0.85),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_color_manual(name = "Region",
                     values = c("#ffa500ff","#3a2b59ff","#92b145ff",
                                "grey","#B21807","#8a85edff"))


  
