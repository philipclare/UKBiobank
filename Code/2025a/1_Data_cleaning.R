######################################################################################
##   
## Analysis of light physical activity and survival
## Clean data
## By: Susan Luo & Philip Clare
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
libs <- c("lubridate","sjmisc","expss","tidyverse","haven","matrixStats","Hmisc","margins",
          "lattice","survival","arsenal","survminer","coxme","parallel","splines","prediction")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

rm(libs)
rm(missing)

######################################################################################
# 2. Load data
#-------------------------------------------------------------------------------------

# 2.1. Load the data
data_ACC_PA <- read.csv("~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/OneDrive/UK Biobank - Melody/data from Liangkai/UK Biobank - acc PA/PA & GRS & T2D/data and R file/PA_data_final.csv") %>%
  select(ID, seasonality, valid_weartime_indays, avg_acc_06_22, contains("_06_22_ML")) %>%
  rename(n_eid = ID) ## 92,184 people at acc measurement

# 2.2. Load the data - dignosistic data
data_diag <- read.csv(file = "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Post-PhD/UK Biobank/Couple research/Data and R/diagnosetime_acc.csv", header = TRUE)

# 2.3. Load the data - covariates data
data_cov <- read.csv(file = "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Post-PhD/UK Biobank/ukb670620tab&r/data_bsl_covariates.csv", header = TRUE) %>% 
  rename(n_eid = f.eid)

# 2.4. combine the data together
data <- data_cov %>% 
  inner_join(data_ACC_PA, by = "n_eid", keep = FALSE) %>%
  inner_join(data_diag, by = "n_eid", keep = FALSE)  # n=92,171

rm(data_cov, data_ACC_PA, data_diag)

######################################################################################
# 3. Clean data
#-------------------------------------------------------------------------------------

# 3.1. Participants exclusion
data <- data %>%
  filter(CVDbase==0, cancerbase==0)  #n=71,878

data$center_0 <- as.factor(data$center_0)
data$sex <- as.factor(data$sex)
data$ethnicity <- as.factor(data$ethnicity)
data$edu_0_cat <- as.factor(data$edu_0_cat)
data$HHincome_0_cat <- as.factor(data$HHincome_0_cat)
data$employment_0_cat <- as.factor(data$employment_0_cat)
data$HealthRating_0_cat <- as.factor(data$HealthRating_0_cat)

######################################################################################
data <- data %>%
  mutate(MVPA_MINS_06_22_ML = MVPA_hrs_06_22_ML*60,
         LPA_MINS_06_22_ML = LPA_hrs_06_22_ML*60)
######################################################################################
### results for the tables
data$MVPA_cat <- cut(
  data$MVPA_MINS_06_22_ML,
  breaks = quantile(data$MVPA_MINS_06_22_ML, c(0, 0.25, 0.5, 0.75, 1)),
  right = TRUE,
  include.lowest = TRUE)
data$LPA_cat <- cut(
  data$LPA_MINS_06_22_ML,
  breaks = quantile(data$LPA_MINS_06_22_ML, c(0, 0.25, 0.5, 0.75, 1)),
  right = TRUE,
  include.lowest = TRUE)

cleaned_data <- data %>% 
  select(-f.2188.0.0, -f.826.0.0)

######################################################################################
# 4. Save data
#-------------------------------------------------------------------------------------

saveRDS(cleaned_data, paste0(workdir,"cleaned_data.rds"))


######################################################################################
# 5. Generate the datasets for sensitivity analyses
data_for_sensitivity_poorhealthonly <- 
  data %>% 
  filter(HealthRating_0_cat!=4 | is.na(HealthRating_0_cat)) %>% 
  select(-f.2188.0.0, -f.826.0.0)

data_for_sensitivity_disability_poorhealth <- 
  data %>% 
  filter(HealthRating_0_cat!=4 | is.na(HealthRating_0_cat)) %>%
  filter(f.2188.0.0!=1 | is.na(f.2188.0.0)) %>% 
  select(-f.2188.0.0, -f.826.0.0)

data_for_sensitivity_shiftworker <- 
  data %>% 
  filter(f.826.0.0==-3 | f.826.0.0==-1 | f.826.0.0==1 | is.na(f.826.0.0)) %>% 
  select(-f.2188.0.0, -f.826.0.0)

saveRDS(data_for_sensitivity_poorhealthonly, "data_for_sensitivity_poorhealthonly.rds")
saveRDS(data_for_sensitivity_disability_poorhealth, "data_for_sensitivity_disability_poorhealth.rds")
saveRDS(data_for_sensitivity_shiftworker, "data_for_sensitivity_shiftworker.rds")

