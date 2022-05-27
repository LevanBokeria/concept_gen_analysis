# Playing around with ANOVA tests
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)


# Clear the environment
rm(list = ls())

# Load the excel file
data_all <-
        read_excel(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed.xlsx'
        ) %>%
        as_tibble()


data_all <- data_all %>%
        mutate(congruency=as.factor(congruency),
               concept_phase_1=as.factor(concept_phase_1),
               concept_phase_2=as.factor(concept_phase_2),
               arr_phase_1_1=as.factor(arr_phase_1_1))

levels(data_all$arr_phase_1_1)

# Summary statistics
data_all %>%
        group_by(congruency,arr_phase_1_1) %>%
        get_summary_stats(phase_2_min_phase_1_ses_1_2_perf, type = 'mean_sd')

# Boxplot
ggboxplot(data_all, x = "congruency", y = "phase_2_min_phase_1_ses_1_2_perf")

# Check for outliers
data_all %>%
        group_by(congruency) %>%
        identify_outliers(phase_2_min_phase_1_ses_1_2_perf) %>%
        View()

ggqqplot(data_all, "phase_2_min_phase_1_ses_1_2_perf", facet.by = "congruency")

# Do the 3-way anova
res.aov_3 <- data_all %>% anova_test(phase_2_min_phase_1_ses_1_2_perf ~ 
                                           congruency*concept_phase_1*arr_phase_1_1)
res.aov_3

res.aov_2 <- data_all %>% anova_test(phase_2_min_phase_1_ses_1_2_perf ~ 
                                             congruency*arr_phase_1_1)
res.aov_2

bxp <- ggboxplot(
        data_all, x = "arr_phase_1_1", y = "phase_2_min_phase_1_ses_1_2_perf",
        color = "congruency", palette = "jco"
)
bxp
