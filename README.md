# UK Biobank Project

This repository contains R code used in articles using the UK Biobank Study.

## LPA and MVPA analysis
Code for all analysis in the article examining the associations of LPA and MVPA with mortality (Luo, Clare et al 2025).

| Description | R Code |
| --- | --- |
| 1 - Data extraction and cleaning. | [Data cleaning](Code/2025a/1_Data_cleaning.R) |
| 2 - Finalise data in preparation for imputation. | [Finalise data](Code/2025a/2_Data_finalise.R) |
| 3 - Multiple imputation - impute intermittent missing data | [Imputation](Code/2025a/3_Multiple_imputation.R) |
| 4 - Final data creation - finalise imputed data and structure for analysis | [Finalise data](Code/2025a/4_Finalise_imputed_data.R) |
| 5 - Evaluate number of degrees of freedom for splines. | [Model fit](Code/2025a/5_Model_fit.R) |
| 6 - Fit best performing model based on model fit tests. | [Final models](Code/2025a/6_Final_model_estimation.R) |
| 7 - Create primary analysis plots using ggplot2. | [Primary figures](Code/2025a/7_Primary_figures.R) |
| 8 - Create sensitivity analysis plots using ggplot2. | [Sensitivity figures](Code/2025a/8_Sensitivity_figures.R) |
| 9 - Estimate relative effects (RRs and RDs) from primary analysis. | [Relative results](Code/2025a/9_Relative_results_HPC.R) |
| 10 - Create relative effect plots. | [Relative results figures](Code/2025a/10_Relative_results_figures.R) |
| 11 - Stratified analysis by sex and age group. | [Age and sex](Code/2025a/11_By_age_and_sex.R) |
| 12 - Relative effects for specific counterfactual scenarios. | [Scenarios](Code/2025a/12_Scenario_comparisons.R) |

