######################################################################################
##   
## Analysis of light physical activity and survival
## Combine imputations and finalise data for analysis
## By: Philip Clare
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
libs <- c("dplyr","fastDummies")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

rm(libs)
rm(missing)

######################################################################################
# 2. Load imputed data into list and process each
#-------------------------------------------------------------------------------------

primary_data <- lapply(seq(1,20), function (x) {
  
  data <- readRDS(file=paste0(workdir,"Data/Imputed/Imputed_data_primary_",x,".rds"))
  
  data <- dummy_cols(data, select_columns = c("seasonality","sex","ethnicity","b_educ","b_income","b_health","b_employ","b_bmi","b_smk","b_alc"),
                     remove_first_dummy = TRUE,
                     remove_selected_columns = TRUE)
  
  data$logage <- log(data$age_acc)
  
  data
  
})

s1_data <- lapply(seq(1,20), function (x) {
  
  data <- readRDS(file=paste0(workdir,"Data/Imputed/Imputed_data_s1_",x,".rds"))
  
  data <- dummy_cols(data, select_columns = c("seasonality","sex","ethnicity","b_educ","b_income","b_health","b_employ","b_bmi","b_smk","b_alc"),
                     remove_first_dummy = TRUE,
                     remove_selected_columns = TRUE)
  
  data$logage <- log(data$age_acc)
  
  data
  
})

s2_data <- lapply(seq(1,20), function (x) {
  
  data <- readRDS(file=paste0(workdir,"Data/Imputed/Imputed_data_s2_",x,".rds"))
  
  data <- dummy_cols(data, select_columns = c("seasonality","sex","ethnicity","b_educ","b_income","b_health","b_employ","b_bmi","b_smk","b_alc"),
                     remove_first_dummy = TRUE,
                     remove_selected_columns = TRUE)
  
  data$logage <- log(data$age_acc)
  
  data
  
})

s3_data <- lapply(seq(1,20), function (x) {
  
  data <- readRDS(file=paste0(workdir,"Data/Imputed/Imputed_data_s3_",x,".rds"))
  
  data <- dummy_cols(data, select_columns = c("seasonality","sex","ethnicity","b_educ","b_income","b_health","b_employ","b_bmi","b_smk","b_alc"),
                     remove_first_dummy = TRUE,
                     remove_selected_columns = TRUE)
  
  data$logage <- log(data$age_acc)
  
  data
  
})

s4_data <- lapply(1, function (x) {
  
  data <- readRDS(file=paste0(workdir,"Data/complete_case_data.rds"))
  
  data <- dummy_cols(data, select_columns = c("seasonality","sex","ethnicity","b_educ","b_income","b_health","b_employ","b_bmi","b_smk","b_alc"),
                     remove_first_dummy = TRUE,
                     remove_selected_columns = TRUE)
  
  data$logage <- log(data$age_acc)
  
  data
  
})

######################################################################################
# 5. Save final imputed data for analysis
#-------------------------------------------------------------------------------------

# 5.1 Save data as imputation list ready for analysis
saveRDS(primary_data,file=paste0(workdir,"Data/primary_analysis_data.rds"))
saveRDS(s1_data,file=paste0(workdir,"Data/s1_analysis_data.rds"))
saveRDS(s2_data,file=paste0(workdir,"Data/s2_analysis_data.rds"))
saveRDS(s3_data,file=paste0(workdir,"Data/s3_analysis_data.rds"))
saveRDS(s4_data,file=paste0(workdir,"Data/s4_analysis_data.rds"))
