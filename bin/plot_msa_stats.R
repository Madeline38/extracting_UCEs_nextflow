#!/usr/bin/env Rscript
# where : Process PLOT_MSA_STATS : plot_msa_stats.R ${summary_stats} ${n_samples}
# usage : This script generates plots summarizing multiple sequence alignment (MSA) statistics. It creates boxplots for various alignment metrics, distribution of the number of taxa per locus, and performs PCA analysis on the data.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args      <- commandArgs(trailingOnly = TRUE)
stats     <- args[1]
n_samples <- as.numeric(args[2])

theme_set(theme_bw())

vars_comp <- c("No_of_taxa", "Alignment_length", "Missing_percent",
               "Proportion_variable_sites", "No_variable_sites",
               "Parsimony_informative_sites",
               "Proportion_parsimony_informative", "GC_content")

d <- read.delim(stats, check.names = FALSE) %>%
  mutate(Prop_taxa  = No_of_taxa / n_samples,
         taxa_class = cut(Prop_taxa, breaks = c(0, .5, .75, 1.0001),
                          labels = c("<50 %", "50-75 %", ">75 %"),
                          include.lowest = TRUE))

stopifnot(all(vars_comp %in% names(d)))

## 1. Boxplots de toutes les variables
d_long <- d %>% select(all_of(vars_comp)) %>% pivot_longer(everything())

p_box <- ggplot(d_long, aes(x = "", y = value)) +
  geom_jitter(width = .25, alpha = .15, size = .6, colour = "grey40") +
  geom_boxplot(fill = "#4FB0E8", alpha = .7, outlier.shape = NA) +
  facet_wrap(~name, scales = "free_y") +
  labs(x = NULL, y = NULL,
       title = paste0("Résumé des alignements (", nrow(d), " loci)"))
ggsave("msa_stats_boxplots.png", p_box, width = 11, height = 8, dpi = 200)

## 2. Distribution du nombre de taxons par locus
p_taxa <- ggplot(d, aes(No_of_taxa)) +
  geom_bar(fill = "#4FB0E8") +
  labs(x = "Nombre de taxons par locus", y = "Nombre de loci")
ggsave("msa_stats_taxa.png", p_taxa, width = 8, height = 4.5, dpi = 200)

## 3. Variables clés par classe de complétude (% d'échantillons présents)
p_class <- d %>%
  select(taxa_class, Alignment_length, Missing_percent,
         Proportion_variable_sites, Proportion_parsimony_informative,
         Parsimony_informative_sites, GC_content) %>%
  pivot_longer(-taxa_class) %>%
  ggplot(aes(taxa_class, value)) +
  geom_boxplot(fill = "#4FB0E8", alpha = .7, outlier.size = .4) +
  facet_wrap(~name, scales = "free_y") +
  labs(x = "Complétude du locus (% d'échantillons)", y = NULL)
ggsave("msa_stats_par_completude.png", p_class, width = 11, height = 7, dpi = 200)

## 4. ACP
tryCatch({
d_pca <- d %>%
  mutate(log_length  = log10(Alignment_length),
         log_missing = log1p(Missing_percent)) %>%
  select(No_of_taxa, log_length, log_missing, Proportion_variable_sites,
         Parsimony_informative_sites, No_variable_sites,
         Proportion_parsimony_informative, GC_content)

d_pca <- d_pca[complete.cases(d_pca), ]
col_taxa <- d_pca$No_of_taxa                                  # sauvegarde avant filtrage
d_pca <- d_pca[, apply(d_pca, 2, sd) > 0, drop = FALSE]       # retire les colonnes constantes

if (nrow(d_pca) >= 5 && ncol(d_pca) >= 3) {
  pc  <- prcomp(d_pca, scale. = TRUE)
  ve  <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  lab <- function(i) sprintf("PC%d (%.1f %%)", i, ve[i])
  sc  <- as.data.frame(pc$x[, 1:2])
  sc$No_of_taxa <- col_taxa

  p_pca <- ggplot(sc, aes(PC1, PC2, colour = No_of_taxa)) +
    geom_point(alpha = .6, size = 1.2) +
    scale_colour_viridis_c() +
    labs(x = lab(1), y = lab(2), title = "ACP des loci")
  ggsave("msa_stats_pca.png", p_pca, width = 7, height = 6, dpi = 200)

  cor_var <- sweep(pc$rotation, 2, pc$sdev, "*")
  arr   <- data.frame(var = rownames(cor_var), PC1 = cor_var[, 1], PC2 = cor_var[, 2])
  theta <- seq(0, 2 * pi, length.out = 200)
  p_circle <- ggplot() +
    geom_path(data = data.frame(x = cos(theta), y = sin(theta)), aes(x, y), colour = "grey70") +
    geom_segment(data = arr, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(.2, "cm"))) +
    geom_text(data = arr, aes(PC1, PC2, label = var), vjust = -.6, size = 3) +
    coord_fixed() +
    labs(x = lab(1), y = lab(2), title = "Cercle des corrélations")
  ggsave("msa_stats_pca_cercle.png", p_circle, width = 6.5, height = 6.5, dpi = 200)
} else {
  message("Trop peu de loci ou de variables pour l'ACP : étape ignorée")
}
}, error = function(e) message("ACP ignorée : ", conditionMessage(e)))
