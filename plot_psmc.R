library(ggplot2)
library(scales)
library(dplyr)
library(purrr)
library(stringr)
setwd("C:/Users/karol/Documents/Sharks/psmc")

# Define population prefixes
pop_prefixes <- c("RS20.100bootplot_mu=1.9e-8.psmc", "Nin84.100bootplot_mu=1.9e-8.psmc", "MAL19.100bootplot_mu=1.9e-8.psmc")

pop_colors <- c("RS20.100bootplot_mu=1.9e-8.psmc" = "#94b6efff", "Nin84.100bootplot_mu=1.9e-8.psmc" = "#d7c743ff", "MAL19.100bootplot_mu=1.9e-8.psmc" = "#f96e19ff")


main_data <- pop_prefixes %>%
  map_df(~ {
    file_path <- file.path(data_dir, paste0(.x, ".0.txt"))
    
    if (file.exists(file_path)) {
      df <- read.table(file_path, header = FALSE)
      df <- df[, 1:2]
      colnames(df) <- c("Time", "Ne")
      df$Population <- .x
      return(df)
    } else {
      warning(paste("Main file not found for:", .x))
      return(NULL)
    }
  })


boot_data <- pop_prefixes %>%
  map_df(~ {
    current_pop <- .x
    # Find all bootstrap files for this specific prefix
    files <- list.files(path = data_dir, 
                        pattern = paste0("^", current_pop, "\\.[1-9][0-9]*\\.txt$"), 
                        full.names = TRUE)
    
    files %>%
      map_df(~ {
        df <- read.table(.x, header = FALSE)
        df <- df[, 1:2]
        colnames(df) <- c("Time", "Ne")
        
        # Extract the unique bootstrap number from the file name
        file_name <- basename(.x)
        boot_num <- str_extract(file_name, "(?<=\\.)[0-9]+(?=\\.txt$)")
        
        df$Population <- current_pop
        df$Group_ID <- paste0(current_pop, "_", boot_num) 
        
        return(df)
      })
  })


ggplot() +
  geom_step(
    data = boot_data, 
    aes(x = Time, y = Ne, group = Group_ID, color = Population), 
    linewidth = 0.3, 
    alpha = 0.15
  ) +
  
  geom_step(
    data = main_data, 
    aes(x = Time, y = Ne, color = Population), 
    linewidth = 1.4
  ) +
  
  scale_color_manual(values = pop_colors) +
  scale_x_log10(
    labels = comma,
    limits = c(1000, 10000000) 
  ) + 
  scale_y_continuous(
    labels = comma
  ) +
  
  labs(
    title =  expression(paste("PSMC of historical ", N[e], " changes")) ,
    #subtitle = "Thick lines indicate median estimates; faint lines indicate 100 bootstrap replicates",
    x =  expression(paste("Years Before Present (g = 16.4, mu = 1.9 x ", 10^-8, ")")),
    y =  expression(paste("Effective Population Size ( ", N[e], " x ", 10^4, ")")),
    color = "Population"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey98")
  )


