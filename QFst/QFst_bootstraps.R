# Population bootstrap for median QFST (95% percentile CIs)
# Modified from Carolin to rely mostly on base R and 
# resample pops within basins.
# Apart from within basin, also get a mean overall
set.seed(123)
nboot <- 10000

# 1. Read data 
DartSeq <- read.table("AngsdWFst_dartseq.txt", header = TRUE)
RADseq <- read.table("AngsdWFST_radseq.txt", header = TRUE)

# 1. Identify basins to compare QFST

# DArTseq
basins_DartSeq <- list(
  CIO = c("chagos", "mal"),
  EIO = c("rowley", "ningaloo", "scott"),
  WPA = c("indo", "south_gbr", "north_gbr", "chesterfield", "walpole",
          "south_nc", "matthew", "grand_astrolabe", "entrecasteaux", "gln")
)
# RADseq
basins_RADseq <- list(
  WIO = c("Juan", "Zelee"),
  PA = c("Bampton", "Niku", "Palmyra", "Fakarava")
)


prepare_data <- function(df, basin_list, dataset) {
  basin_lookup <- setNames(rep(names(basin_list), lengths(basin_list)),
                           unlist(basin_list, use.names = FALSE))
  df$basin1 <- unname(basin_lookup[df$pop1])
  df$basin2 <- unname(basin_lookup[df$pop2])
  excluded <- is.na(df$basin1) | is.na(df$basin2)
  if (any(excluded)) {
    message(dataset, ": excluding ", sum(excluded), " pairs with unmapped populations")
  }
  df <- df[!excluded, ]
  # -log1p(-x),  numerically stable for FST close to zero.
  df$Qfst <- -log1p(-df$Fst_autos) / -log1p(-df$Fst_X)
  df <- df[is.finite(df$Qfst) & is.finite(df$Fst_autos) &
             is.finite(df$Fst_X), ]
  if (dataset == "DartSeq") {
    #avoids  dropping reversed comparisons.
    df$comparison <- ifelse(df$basin1 == df$basin2, "within",
      ifelse((df$basin1 == "CIO" & df$basin2 == "EIO") |
             (df$basin1 == "EIO" & df$basin2 == "CIO"), "CIO-EIO",
      ifelse((df$basin1 == "CIO" & df$basin2 == "WPA") |
             (df$basin1 == "WPA" & df$basin2 == "CIO"), "CIO-WPA", "EIO-WPA")))
  } else {
    # Preserve the original script's combined within-basin category.
    df$comparison <- ifelse(df$basin1 == df$basin2, "within", "WIO-PA")
  }
  df
}

# 2. Weighted median: equivalent to repeating each pair according to its weight
weighted_median <- function(x, w) {
  keep <- is.finite(x) & w > 0
  x <- x[keep]; w <- w[keep]
  if (!length(x)) return(NA_real_)
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cum <- cumsum(w)
  half <- sum(w) / 2
  # Match median(rep(x, w)), including averaging two central values.
  (x[which(cum >= half)[1]] + x[which(cum > half)[1]]) / 2
}

# 3. Bootstrap populations within basins. In addition to the regional
#    comparisons, calculate an overall median using ALL available pairs of course excluding within basin comparison.

bootstrap_dataset <- function(df, dataset, nboot) {
  pop_info <- unique(rbind(
    data.frame(pop = df$pop1, basin = df$basin1),
    data.frame(pop = df$pop2, basin = df$basin2)
  ))
  stopifnot(!anyDuplicated(pop_info$pop))
  basin_pops <- split(pop_info$pop, pop_info$basin)
  measures <- c("Qfst", "Fst_X", "Fst_autos")
  # Include within-basin pairs ONLY in the overall result for our dataset.
  comparisons <- c(setdiff(unique(df$comparison), if (dataset == "DartSeq") "within" else character(0)), "ALL")
  boot <- array(NA_real_, c(nboot, length(comparisons), length(measures)),
                dimnames = list(NULL, comparisons, measures))

  # Each replicate is sampled ONCE and reused for all statistics/comparisons.
  for (b in seq_len(nboot)) {
    sampled <- unlist(lapply(basin_pops, function(p) sample(p, length(p), replace = TRUE)),
                      use.names = FALSE)
    counts <- table(factor(sampled, levels = pop_info$pop))
    weights <- as.numeric(counts[df$pop1]) * as.numeric(counts[df$pop2])
    for (comp in comparisons) {
      rows <- if (comp == "ALL") df$basin1 != df$basin2 else df$comparison == comp
      for (measure in measures) {
        boot[b, comp, measure] <- weighted_median(df[[measure]][rows], weights[rows])
      }
    }
  }

  out <- list()
  for (comp in comparisons) {
    rows <- if (comp == "ALL") rep(TRUE, nrow(df)) else df$comparison == comp
    for (measure in measures) {
      valid <- boot[, comp, measure]
      valid <- valid[is.finite(valid)]
      if (length(valid) < 100) stop("Too few valid replicates: ", dataset, " / ", comp)
      out[[length(out) + 1]] <- data.frame(
        dataset = dataset, comparison = comp, statistic = measure,
        observed = median(df[[measure]][rows]),
        CI_lower = unname(quantile(valid, 0.025)),
        CI_upper = unname(quantile(valid, 0.975)),
        n_pairs = sum(rows), n_valid = length(valid), nboot = nboot
      )
    }
  }
  do.call(rbind, out)
}

DartSeq <- prepare_data(DartSeq, basins_DartSeq, "DartSeq")
RADseq <- prepare_data(RADseq, basins_RADseq, "RADseq")
results <- rbind(
  bootstrap_dataset(DartSeq, "DartSeq", nboot),
  bootstrap_dataset(RADseq, "RADseq", nboot)
)

# 4. Save  medians and CIs.
print(results, row.names = FALSE)
write.csv(results, "Qfst_bootstrap_95CI.csv", row.names = FALSE)

# 5. QFST plot 

custom_colors <- c(
  "EIO-WPA" = "#8DB343",
  "WIO-PA" = "#B21807",
  "CIO-EIO"=  "#E1A12B",
  "CIO-WPA" = "#AE85AD"
)

datapoints<-rbind(RADseq,DartSeq)

ggplot(
 datapoints %>% filter(!is.na(Qfst) & comparison != "within"),
  aes(x = comparison, y = Qfst, fill = comparison, color = comparison)
) +
    geom_jitter(
    aes(x = comparison),
    width = 0.05,
    height = 0,
    size = 3,
    alpha = 0.8
  ) + 
  # Means and CI
  geom_errorbar(
    data = q[q$comparison!="ALL",],
    aes(x = comparison, y=observed, ymin = CI_lower, ymax = CI_upper),
    width = 0.15, linewidth = 0.8,  position = position_nudge(x = -0.18)
  ) +
  geom_point(
    data = q[q$comparison != "ALL", ],
    aes(x = comparison, y = observed, shape = dataset),
    size = 4.5,  position = position_nudge(x = -0.18)
  ) +
  scale_shape_manual(values = c("DartSeq" = 16, "RADseq" = 15)) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
  labs(
    title = expression("Regional QF"["ST"] ~ "comparisons"),
    x = "Comparison",
    y = expression(paste("QF"["ST"]))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12)
  )+
  geom_hline(yintercept = 0.75, color="grey", linetype="dashed")


