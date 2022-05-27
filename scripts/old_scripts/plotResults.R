# Load libraries
library(tidyverse)
library(readxl)
library(ggpubr)
library(rstatix)
library(emmeans)


# Clear the environment
rm(list = ls())

# Load the results_table_qc_pass_ptp_analyzed
# results_table_qc_pass_ptp_analyzed <-
#         read_excel(
#                 'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed.xlsx'
#         ) %>%
#         as_tibble()

# Load the balanced csv file
results_table_qc_pass_ptp_analyzed <-
        read_csv(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed_balanced.csv'
        ) %>%
        as_tibble()


results_table_qc_pass_ptp_analyzed <- results_table_qc_pass_ptp_analyzed %>%
        mutate(congruency=as.factor(congruency),
               concept_phase_1=as.factor(concept_phase_1),
               concept_phase_2=as.factor(concept_phase_2),
               arr_phase_1_1=as.factor(arr_phase_1_1))

# Load the long_form_data_all_ptp
long_form_data_all_ptp <-
        read_excel(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/long_form_data_all_ptp.xlsx'
                , sheet = 1
        ) %>%
        as_tibble()

long_form_data_all_ptp <- long_form_data_all_ptp %>%
        mutate(congruency=as.factor(congruency),
               concept_phase_1=as.factor(concept_phase_1),
               concept_phase_2=as.factor(concept_phase_2),
               arr_phase_1_1=as.factor(arr_phase_1_1))

long_form_data_qc_pass_ptp <- subset(long_form_data_all_ptp, 
                                     data_submitted == 1 & 
                                     global_pass_incl_phase_2_fails == 1 & 
                                     debrief_qc_pass == 1 & 
                                     fb_int_qc_pass == 1 & 
                                     phase_1_rt_qc_pass == 1 & 
                                     phase_2_rt_qc_pass == 1)




# ###################### FILTER THE DATA ###################################
long_form_data_qc_pass_ptp <-
        subset(long_form_data_qc_pass_ptp,
                                     arr_phase_1_1 == 3 | arr_phase_1_1 == 15)
results_table_qc_pass_ptp_analyzed <-
        subset(results_table_qc_pass_ptp_analyzed,
                                             arr_phase_1_1 == 3 | arr_phase_1_1 == 15)


###################### GATHER ##############################################
# Gather the long form data for results_table_qc_pass_ptp_analyzed
results_table_qc_pass_ptp_analyzed_gathered <- 
        results_table_qc_pass_ptp_analyzed %>%
        select(ptp,congruency,phase_1_ses_1_2_perf,phase_2_ses_1_2_perf) %>%
        gather(key='phase',value='ses_1_2_perf',
               phase_1_ses_1_2_perf,
               phase_2_ses_1_2_perf)


##################### start plots ###########################################



# Plot the basic congruency effect, collapsing everything.
results_table_qc_pass_ptp_analyzed %>%
        ggplot(aes(x=congruency,y=phase_2_min_phase_1_ses_1_2_perf)) + 
        geom_boxplot(width=0.5,notch = TRUE) + 
        geom_jitter(width=0.1) + 
        ylab("P2-P1, ses 1&2 acc") + 
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14))

# Show phases separately
results_table_qc_pass_ptp_analyzed_gathered %>%
        ggplot(aes(x=congruency,y=ses_1_2_perf,fill=phase)) +
        geom_boxplot(width=0.5,notch=TRUE, outlier.shape = " ") +
        geom_point(pch=21,
                   position = position_jitterdodge(dodge.width=0.5,
                                                   jitter.width=0.1)) +
        ylab("Ses 1&2 acc") + 
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),
              legend.position = 'top') + 
        scale_fill_manual(values=c('#2c7fb8','orange'),
                          guide=guide_legend(title='',
                                             nrow = 1))


# Show the congruency effect for arrangements 
results_table_qc_pass_ptp_analyzed %>%
        ggplot(aes(x=arr_phase_1_1,
                   y=phase_2_min_phase_1_ses_1_2_perf,
                   fill=congruency)) + 
        geom_boxplot(width=0.5,notch = FALSE) + 
        geom_point(pch=21,
                   position = position_jitterdodge(dodge.width=0.5,
                                                   jitter.width=0.1)) +
        ylab("P2-P1, ses 1&2 acc") +
        xlab("Arrangements") + 
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),
              legend.position = 'top') + 
        # scale_x_discrete(labels=c('Arr1','Arr2','Arr3','Arr4')) + 
        scale_fill_manual(values=c('#e5f5f9','#2ca25f'),
                          guide=guide_legend(title='Congruency',
                                             nrow = 1))


# Show arrangement by congruency insteraction
results_table_qc_pass_ptp_analyzed %>%
        select(ptp,congruency,arr_phase_1_1,phase_2_min_phase_1_ses_1_2_perf) %>%
        reorder_levels(arr_phase_1_1, order = c("1", "14", "3", "15")) %>%
        ggplot(aes(x=arr_phase_1_1,
                   y=phase_2_min_phase_1_ses_1_2_perf,
                   fill=congruency)) +
        geom_boxplot(width=0.5,notch=FALSE,outlier.shape = " ") +
        geom_point(pch=21,
                   position = position_jitterdodge(dodge.width=0.5,
                                                   jitter.width=0.2)) +
        ylab("P2-P1 ses 1&2 perf") + 
        xlab("Arrangements") + 
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),
              legend.position = 'top') + 
        scale_x_discrete(labels=c('Arr1','Arr2','Arr3','Arr4')) +
        scale_fill_manual(values=c('#e5f5f9','#2ca25f'),
                          guide=guide_legend(title='Congruency',
                                             nrow = 1))



levels(results_table_qc_pass_ptp_analyzed$arr_phase_1_1)

# Summary statistics
results_table_qc_pass_ptp_analyzed %>%
        group_by(congruency,arr_phase_1_1) %>%
        get_summary_stats(phase_2_min_phase_1_ses_1_2_perf, type = 'mean_sd')

# Boxplot
ggboxplot(results_table_qc_pass_ptp_analyzed, x = "congruency", y = "phase_2_min_phase_1_ses_1_2_perf")

# Check for outliers
results_table_qc_pass_ptp_analyzed %>%
        group_by(congruency) %>%
        identify_outliers(phase_2_min_phase_1_ses_1_2_perf) %>%
        View()

ggqqplot(results_table_qc_pass_ptp_analyzed, "phase_2_min_phase_1_ses_1_2_perf", facet.by = "congruency")

# Do the 3-way anova
res.aov_3 <- results_table_qc_pass_ptp_analyzed %>% anova_test(phase_2_min_phase_1_ses_1_2_perf ~ 
                                           congruency*concept_phase_1*arr_phase_1_1)
res.aov_3

res.aov_2 <- results_table_qc_pass_ptp_analyzed %>% anova_test(phase_2_min_phase_1_ses_1_2_perf ~ 
                                             congruency*arr_phase_1_1)
res.aov_2

bxp <- ggboxplot(
        results_table_qc_pass_ptp_analyzed, x = "arr_phase_1_1", y = "phase_2_min_phase_1_ses_1_2_perf",
        color = "congruency", palette = "jco"
)
bxp
