library(cardx)
library(dplyr)
library(DT)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(gtsummary)
library(grid)
library(gridExtra)
library(mia)
library(miaViz)
library(miaTime)
library(multtest)
library(parameters)
library(patchwork)
library(quarto)
library(scater)
library(sechm)
library(tidyr)
library(TreeSummarizedExperiment)
library(tidyverse)
library(reshape2)
library(vegan)
library(lmerTest)
library(kableExtra)

outdir ="./output/"

run_diversity_tests <- function(tse, comparisons, 
                                group.var, index, 
                                adjust.method, paired = FALSE) {
  # Get all values and groups
  values <- colData(tse)[[index]]
  groups <- colData(tse)[[group.var]]
  
  # Create an empty results dataframe
  results <- data.frame(
    group1 = character(),
    group2 = character(),
    mean1 = numeric(),
    mean2 = numeric(),
    logFC = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through requested comparisons
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    
    # Subset data for just these two groups
    subset_indices <- which(groups %in% c(group1, group2))
    subset_values <- values[subset_indices]
    subset_groups <- groups[subset_indices]
    
    # Ensure paired data is properly aligned (if paired is TRUE)
    if (paired) {
      # Assuming 'subject' column exists for pairing
      subset_subjects <- tse$id[subset_indices]
      paired_data <- data.frame(
        subject = subset_subjects,
        group = subset_groups,
        value = subset_values
      ) %>%
        tidyr::pivot_wider(names_from = group, values_from = value) %>%
        drop_na()  # Drop pairs with missing data
      
      # Extract paired values
      paired_values <- paired_data %>% select(all_of(c(group1, group2)))
    }
    
    # Perform Wilcoxon test
    if (paired) {
      test_result <- wilcox.test(
        paired_values[[1]], paired_values[[2]],
        paired = TRUE, exact = FALSE
      )
    } else {
      test_result <- wilcox.test(
        subset_values ~ subset_groups,
        exact = FALSE
      )
    }
    
    # Calculate means
    mean1 <- mean(values[groups == group1])
    mean2 <- mean(values[groups == group2])
    
    if (paired) {
      # Note: when the values are paired, we can calculate FC separately
      # for every individual and summarize afterwards. 
      FC <- mean(values[groups == group2] / values[groups == group1])
    } else {
      FC <- mean2/mean1
    }
    
    # Add results
    results <- rbind(results, data.frame(
      group1 = group1,
      group2 = group2,
      mean1 = mean1,
      mean2 = mean2,
      logFC = log2(FC),
      p.value = test_result$p.value,
      stringsAsFactors = FALSE
    ))
    
  }
  
  # Add adjusted p-values using specified method
  # results$p.adj <- p.adjust(results$p.value, method = adjust.method)
  
  return(results)
}

run_lmer <- function(tse, target, formula_used, id.var, group.var, time.var){
  
  df <- colData(tse) %>% as.data.frame
  df$y <- df[[target]]
  df <- df[, c("y", id.var, group.var, time.var)]
  
  # Filter out missing data cases
  df_no_miss <- df[df %>% complete.cases,]
  
  m <- lmerTest::lmer(as.formula(formula_used), data = df_no_miss)
  return(m)
}