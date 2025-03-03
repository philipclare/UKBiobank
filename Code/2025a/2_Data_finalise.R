######################################################################################
##   
## Analysis of light physical activity and survival
## Finalise data, ready for imputation
## By: Philip Clare & Susan Luo
## Date: 29/7/2024
## Licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
## OSF Registration: 10.17605/OSF.IO/5U7CH
##
######################################################################################
# 1. Setup Environment
#-------------------------------------------------------------------------------------

# 1.1. Define location of data
if (Sys.getenv("NCPUS")!="") {
  .libPaths("/home/z3312911/RPackages/")
  workdir <- "/home/z3312911/LPA/"
} else if (Sys.info()[['sysname']]=="Windows") {
  workdir <- "D:/The University of Sydney (Staff)/Susan Luo - data and R/"
} else if (Sys.info()[['sysname']]=="Darwin") {
  workdir <- "/Users/pjclare/Library/CloudStorage/OneDrive-SharedLibraries-TheUniversityofSydney(Staff)/Susan Luo - data and R/" # MAC
}

# 1.2. Define packages to be used, check if installed, and load
libs <- c("dplyr","fastDummies","tidyr")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

rm(libs)
rm(missing)

######################################################################################
# 2. Define function for cleaning data
#-------------------------------------------------------------------------------------

data_clean <- function (data) {
  data <- data %>%
    select(n_eid, center_0, age_acc, sex, ethnicity, edu_0_cat, HHincome_0_cat, TDI_bsl, HealthRating_0_cat, employment_0_cat, 
           smk_0_cat, alc_0_cat, alc_0_freq_cat, diet_score, bmi_0_cat, seasonality, valid_weartime_indays, f_death, death, MVPA_MINS_06_22_ML, 
           LPA_MINS_06_22_ML, MVPA_cat, LPA_cat) %>% 
    filter(f_death > 0) # n=71,875
  
  data$center_0 <- droplevels(data$center_0) 
  
  data <- data %>%
    rename(
      b_educ = edu_0_cat,
      b_income = HHincome_0_cat,
      b_health = HealthRating_0_cat,
      b_employ = employment_0_cat,
      b_bmi = bmi_0_cat,
      b_smk = smk_0_cat,
      center = center_0,
      MVPA = MVPA_MINS_06_22_ML,
      LPA = LPA_MINS_06_22_ML
    )
  
  levels(data$b_income)[levels(data$b_income)=='9'] <- NA
  levels(data$b_health)[levels(data$b_health)=='9'] <- NA
  
  data$b_alc <- data$alc_0_cat
  data$b_alc <- ifelse(data$b_alc==2 & data$alc_0_freq_cat<=3 & data$alc_0_freq_cat>0,3,data$b_alc)
  
  data <- data %>%
    mutate(
      b_alc = case_match(b_alc, -3 ~ NA, .default = b_alc),
      b_smk = case_match(b_smk, -3 ~ NA, .default = b_smk)
    )
  
  
  data$b_income <- droplevels(data$b_income) 
  data$b_health <- droplevels(data$b_health) 
  
  data$sex <- factor(data$sex,
                     levels=c(0,1),
                     labels=c("female","male"))
  data$ethnicity <- factor(data$ethnicity,
                           levels=c(1,2),
                           labels=c("white","nonwhite"))
  data$b_employ <- factor(data$b_employ,
                          levels=c(1,2),
                          labels=c("employed","other"))
  data$b_bmi <- factor(data$b_bmi,
                       levels=c(2,1,3,4),
                       labels=c("healthy","underweight","overweight","obese"))
  data$b_educ <- factor(data$b_educ,
                        levels=c(1,2,3,4,5),
                        labels=c("university","vocational","olevels","alevels","other"))
  data$b_income <- factor(data$b_income,
                          levels=c(1,2,3,4,5),
                          labels=c("lessthan18k","18to31k","31to52k","52to100k","100kplus"))
  data$b_health <- factor(data$b_health,
                          levels=c(1,2,3,4),
                          labels=c("excellent","good","fair","poor"))
  data$b_alc <- factor(data$b_alc,
                       levels=c(0,1,2,3),
                       labels=c("never","ex","current_infrequent","current_weekly"))
  data$b_smk <- factor(data$b_smk,
                       levels=c(0,1,2),
                       labels=c("never","ex","current"))
  
  data <- data %>%
    select(-c(MVPA_cat,LPA_cat,alc_0_cat,alc_0_freq_cat))

}

######################################################################################
# 3. Load and finalise data
#-------------------------------------------------------------------------------------

data <- readRDS(paste0(workdir,"Data/cleaned_data.rds")) # n=71,878
data <- data_clean(data)

complete_data <- data[complete.cases(data),] 

complete_data$seasonality <- factor(complete_data$seasonality)

s1 <- readRDS(paste0(workdir,"Data/data_for_sensitivity_shiftworker.rds")) 
s1 <- data_clean(s1)

s2 <- readRDS(paste0(workdir,"Data/data_for_sensitivity_poorhealthonly.rds")) 
s2 <- data_clean(s2)

s3 <- readRDS(paste0(workdir,"Data/data_for_sensitivity_disability_poorhealth.rds")) 
s3 <- data_clean(s3)

######################################################################################
# 4. Save data
#-------------------------------------------------------------------------------------

saveRDS(data, paste0(workdir,"Data/data_for_imputation.rds")) # n=70,891
saveRDS(s1, paste0(workdir,"Data/data_for_imputation_sensitivity1.rds")) # n=64,845
saveRDS(s2, paste0(workdir,"Data/data_for_imputation_sensitivity2.rds")) # n=69,526
saveRDS(s3, paste0(workdir,"Data/data_for_imputation_sensitivity3.rds")) # n=53,371
saveRDS(complete_data, paste0(workdir,"Data/complete_case_data.rds")) # n=64,096
