# Libraries
library(BayesFactor)
library(assortedRFunctions)
library(tidyverse)
library(readxl)
library(viridis)
# library(MASS)

# Clear the environment
rm(list = ls())

# Load the unbalanced excel file
data_all <-
    read_excel(
        'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed.xlsx'
    ) %>%
    as_tibble()

# Load the balanced csv file
# data_all <-
#     read_csv(
#         'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed_balanced.csv'
#     ) %>%
#     as_tibble()

data_all <- data_all %>%
    mutate(congruency=as.factor(congruency),
           concept_phase_1=as.factor(concept_phase_1),
           concept_phase_2=as.factor(concept_phase_2),
           arr_phase_1_1=as.factor(arr_phase_1_1))

# Exclude the possible outlier
# data_all <- subset(data_all,!(ptp == '5b27a2477f78160001464118'))

# Get data from the last two arrangements
# data_all <- subset(data_all, arr_phase_1_1 == 3 | arr_phase_1_1 == 15)
data_all <- subset(data_all, arr_phase_1_1 == 1 | arr_phase_1_1 == 14)

# Analyze learning rate

learning_rate_congruent_1 = data_all$phase_2_min_phase_1_learning_rate_exp[data_all$congruency == 1]
learning_rate_congruent_0 = data_all$phase_2_min_phase_1_learning_rate_exp[data_all$congruency == 0]

bf_learning_rate <-
    reportBF(
        ttestBF(
            learning_rate_congruent_1,
            learning_rate_congruent_0,
            paired = FALSE
        ),
        4
    )

# Do non-parametric test on learning rates 
wilcox.test(learning_rate_congruent_1,learning_rate_congruent_0)




# Analyze session differences
ses_1_2_perf_congruent_1 = data_all$phase_2_min_phase_1_ses_1_2_perf[data_all$congruency == 1]
ses_1_2_perf_congruent_0 = data_all$phase_2_min_phase_1_ses_1_2_perf[data_all$congruency == 0]

bf_ses_1_2_perf <-
    reportBF(ttestBF(
        ses_1_2_perf_congruent_1,
        ses_1_2_perf_congruent_0,
        paired = FALSE
    ),
    4)

t.test(ses_1_2_perf_congruent_0,ses_1_2_perf_congruent_1,paired=FALSE)

# Effect size for the session differences
sd_c1 <- sd(ses_1_2_perf_congruent_1)
sd_c0 <- sd(ses_1_2_perf_congruent_0)

n_c1 <- length(ses_1_2_perf_congruent_1)
n_c0 <- length(ses_1_2_perf_congruent_0)

pooled_sd <- sqrt( 
    (
    (n_c1-1)*(sd_c1^2) + 
    (n_c0-1)*(sd_c0^2)
    ) / (n_c1 + n_c0 - 2)
    )

effect_size <- (mean(ses_1_2_perf_congruent_1) - mean(ses_1_2_perf_congruent_0)) /
    pooled_sd

# # Some plots

ggplot(data=data_all, mapping=aes(x=congruency, 
                                  y=phase_2_min_phase_1_ses_1_2_perf)) +
    geom_boxplot(notch=TRUE) +
    geom_jitter(width=0.1) + 
    ggtitle("A boxplot with jitter")

# Transform to long form
data_all %>%
    select(phase_1_ses_1_2_perf,phase_2_ses_1_2_perf) %>% 
    gather(key = 'Description', value='Accuracy') %>%
    ggplot(aes(x=Description, y=Accuracy)) + 
    geom_boxplot(notch = TRUE) + 
    geom_jitter(width=0.1)

