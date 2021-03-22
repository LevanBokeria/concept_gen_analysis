# Libraries
library(BayesFactor)
library(assortedRFunctions)
library(tidyverse)
library(readxl)
# library(MASS)

# Clear the environment
rm(list = ls())

# Load the excel file
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

# data_all <- subset(data_all,!(ptp == '5b27a2477f78160001464118'))

###############################################################################
###############################################################################
# Concept order

# Split into two concepts
data_nl_c1 <- data_all %>%
    subset(concept_phase_1 == 'neck_legs_space' & congruency == 1)
data_nl_c0 <- data_all %>%
    subset(concept_phase_1 == 'neck_legs_space' & congruency == 0)
data_bt_c1 <- data_all %>%
    subset(concept_phase_1 == 'beak_tail_space' & congruency == 1)
data_bt_c0 <- data_all %>%
    subset(concept_phase_1 == 'beak_tail_space' & congruency == 0)

# Run the bayesian analysis for order of concepts

BF_nl <- ttestBF(
    data_nl_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_nl_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)
BF_bt <- ttestBF(
    data_bt_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_bt_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)

t_test_nl <- t.test(
    data_nl_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_nl_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)
t_test_bt <- t.test(
    data_bt_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_bt_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)

# Plot the differences for order of concepts
ggplot(data=data_all, mapping=aes(x=congruency, 
                                    y=phase_2_min_phase_1_ses_1_2_perf,
                                  color=concept_phase_1)) +
    geom_boxplot()
    # geom_jitter(width=0.1)

###############################################################################
###############################################################################
# Arrangement order 

# Split into all the arrangements
data_arr_1_c1 <- data_all %>%
    subset(arr_phase_1_1 == 1 & congruency == 1)
data_arr_1_c0 <- data_all %>%
    subset(arr_phase_1_1 == 1 & congruency == 0)
data_arr_2_c1 <- data_all %>%
    subset(arr_phase_1_1 == 14 & congruency == 1)
data_arr_2_c0 <- data_all %>%
    subset(arr_phase_1_1 == 14 & congruency == 0)
data_arr_3_c1 <- data_all %>%
    subset(arr_phase_1_1 == 3 & congruency == 1)
data_arr_3_c0 <- data_all %>%
    subset(arr_phase_1_1 == 3 & congruency == 0)
data_arr_4_c1 <- data_all %>%
    subset(arr_phase_1_1 == 15 & congruency == 1)
data_arr_4_c0 <- data_all %>%
    subset(arr_phase_1_1 == 15 & congruency == 0)

# Run the bayesian analysis for order of concepts

BF_arr_1_first <- ttestBF(
    data_arr_1_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_1_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)
BF_arr_2_first <- ttestBF(
    data_arr_2_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_2_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)
BF_arr_3_first <- ttestBF(
    data_arr_3_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_3_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)
BF_arr_4_first <- ttestBF(
    data_arr_4_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_4_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)

t_test_arr_1_first <- t.test(
    data_arr_1_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_1_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)
t_test_arr_2_first <- t.test(
    data_arr_2_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_2_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)
t_test_arr_3_first <- t.test(
    data_arr_3_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_3_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)
t_test_arr_4_first <- t.test(
    data_arr_4_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_4_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)

# Plot the differences for order of concepts
# ggplot(data=data_all, mapping=aes(x=arr_phase_1_1, 
#                                   y=phase_2_min_phase_1_ses_1_2_perf,
#                                   color=congruency)) +
    # geom_violin()
# geom_jitter(width=0.1)



# Split such that incongruent still includes the other arrangements
data_arr_1_2_c0 <- data_all %>%
    subset((arr_phase_1_1 == 1 | arr_phase_1_1 == 14) & congruency == 0)

data_arr_3_4_c0 <- data_all %>%
    subset((arr_phase_1_1 == 15 | arr_phase_1_1 == 3) & congruency == 0)

# Run the bayesian analysis for order of concepts

BF_arr_1_2 <- ttestBF(
    data_arr_1_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_1_2_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)
BF_arr_2_1 <- ttestBF(
    data_arr_2_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_1_2_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)
BF_arr_3_4 <- ttestBF(
    data_arr_3_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_3_4_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)
BF_arr_4_3 <- ttestBF(
    data_arr_4_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_3_4_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
) %>% reportBF(4)

t_test_arr_1_2 <- t.test(
    data_arr_1_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_1_2_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)
t_test_arr_2_1 <- t.test(
    data_arr_2_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_1_2_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)
t_test_arr_3_4 <- t.test(
    data_arr_3_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_3_4_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)
t_test_arr_4_3 <- t.test(
    data_arr_4_c1 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],
    data_arr_3_4_c0 %>% select(phase_2_min_phase_1_ses_1_2_perf) %>% .[[1]],    
    paired = FALSE
)

############################################################################