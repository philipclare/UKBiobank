######################################################################################
##   
## Analysis of light physical activity and survival
## Check distribution and spline df for best fitting model
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
libs <- c("survival","rstpm2","parallel","splines","XLConnect")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

rm(missing)

######################################################################################
# 2. Load and finalise data
#-------------------------------------------------------------------------------------

data <- readRDS(paste0(workdir,"Data/s4_analysis_data.rds"))[[1]]

######################################################################################
# 3. Create cluster for parallel processing
#-------------------------------------------------------------------------------------

cl <- makeCluster(4)
clusterExport(cl,c("libs","data"))
clusterEvalQ(cl, lapply(libs, library, character.only = TRUE))

######################################################################################
# 4. Check degrees of freedom of survival distribution
#-------------------------------------------------------------------------------------

# 4.1. Run all models and save AIC/BIC
df_dist <- seq(1,8)
fit_dist <- do.call(rbind,parLapply(cl, df_dist, function (x) {
  fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = 1) + ns(LPA, df = 1) + ns(MVPA*LPA, df = 1) + 
                 ns(logage, df = 1) + ns(logage*MVPA, df = 1) + ns(logage*LPA, df = 1) + 
                 sex_male + I(sex_male*MVPA) + I(sex_male*LPA) + ethnicity_nonwhite + 
                 b_educ_vocational + b_educ_olevels + b_educ_alevels + b_educ_other + 
                 b_income_18to31k + b_income_31to52k + b_income_52to100k + b_income_100kplus + b_employ_other +
                 b_health_good + b_health_fair + b_health_poor + 
                 b_bmi_underweight + I(b_bmi_underweight*MVPA) + I(b_bmi_underweight*LPA) + 
                 b_bmi_overweight + I(b_bmi_overweight*MVPA) + I(b_bmi_overweight*LPA) + 
                 b_bmi_obese + I(b_bmi_obese*MVPA) + I(b_bmi_obese*LPA) + 
                 b_alc_ex + b_alc_current_infrequent + b_alc_current_weekly + 
                 b_smk_ex + I(b_smk_ex*MVPA) + I(b_smk_ex*LPA) + b_smk_current + 
                 I(b_smk_current*MVPA) + I(b_smk_current*LPA) + diet_score + 
                 seasonality_Spring + seasonality_Summer + seasonality_Winter + 
                 TDI_bsl + valid_weartime_indays + ns(b_met, df = 1) + b_sleep_cat_enough + b_sleep_cat_toomuch +
                 b_diabetes + b_hypertension + b_highcholesterol + b_affective_disorder,
               data = data,
               df = x)
  c(x,AIC(fit),BIC(fit))
  

  
}))

best_dist_aic <- arrayInd(which.min(fit_dist[,2]), dim(fit_dist))
best_dist_bic <- arrayInd(which.min(fit_dist[,3]), dim(fit_dist))

# 4.2. Print the best fitting models based on AIC/BIC
paste0("The best fitting model based on AIC is df=",fit_dist[best_dist_aic[[1]],1])
paste0("The best fitting model based on BIC is df=",fit_dist[best_dist_bic[[1]],1])

df <- fit_dist[best_dist_aic[[1]],1]

######################################################################################
# 5. Check splines for MVPA and LPA
#-------------------------------------------------------------------------------------

# 5.1. Run all models and save AIC/BIC
df_pa <- seq(1,6)
fit_act <- do.call(rbind,parLapply(cl, df_pa, function (x,df_pa,df) {
  do.call(rbind,lapply(df_pa, function (y,x,df_pa,df) {
    do.call(rbind,lapply(df_pa, function (z,y,x,df) {
      fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = x) + 
                     ns(LPA, df = y) + 
                     ns(MVPA*LPA, df = z) +  
                     ns(logage, df = 1) + ns(logage*MVPA, df = 1) + ns(logage*LPA, df = 1) + 
                     sex_male + I(sex_male*MVPA) + I(sex_male*LPA) + ethnicity_nonwhite + 
                     b_educ_vocational + b_educ_olevels + b_educ_alevels + b_educ_other + 
                     b_income_18to31k + b_income_31to52k + b_income_52to100k + b_income_100kplus + b_employ_other +
                     b_health_good + b_health_fair + b_health_poor + 
                     b_bmi_underweight + I(b_bmi_underweight*MVPA) + I(b_bmi_underweight*LPA) + 
                     b_bmi_overweight + I(b_bmi_overweight*MVPA) + I(b_bmi_overweight*LPA) + 
                     b_bmi_obese + I(b_bmi_obese*MVPA) + I(b_bmi_obese*LPA) + 
                     b_alc_ex + b_alc_current_infrequent + b_alc_current_weekly + 
                     b_smk_ex + I(b_smk_ex*MVPA) + I(b_smk_ex*LPA) + b_smk_current + 
                     I(b_smk_current*MVPA) + I(b_smk_current*LPA) + diet_score + 
                     seasonality_Spring + seasonality_Summer + seasonality_Winter + 
                     TDI_bsl + valid_weartime_indays + ns(b_met, df = 1) + b_sleep_cat_enough + b_sleep_cat_toomuch +
                     b_diabetes + b_hypertension + b_highcholesterol + b_affective_disorder,
                   data = data,
                   df = df)
      c(x,y,z,AIC(fit),BIC(fit))
    },x=x,y=y,df=df))
  },x=x,df_pa=df_pa,df=df))
},df_pa=df_pa,df=df))

best_act_aic <- arrayInd(which.min(fit_act[,4]), dim(fit_act))
best_act_bic <- arrayInd(which.min(fit_act[,5]), dim(fit_act))

# 5.2. Print the best fitting models based on AIC/BIC
paste0("The best fitting model based on AIC has df=",fit_act[best_act_aic[[1]],1]," for MVPA, df=",fit_act[best_act_aic[[1]],2]," for LPA, and df=",fit_act[best_act_aic[[1]],3]," for MVPA*LPA")
paste0("The best fitting model based on BIC has df=",fit_act[best_act_bic[[1]],1]," for age, df=",fit_act[best_act_bic[[1]],2]," for age*mvpa, and df=",fit_act[best_act_bic[[1]],3]," for age*lpa")

mvpa_df <- fit_act[best_act_aic[[1]],1]
lpa_df <- fit_act[best_act_aic[[1]],2]
mvpalpa_df <- fit_act[best_act_aic[[1]],3]

######################################################################################
# 6. Check splines for age and age*activity interactions
#-------------------------------------------------------------------------------------

# 6.1. Run all models and save AIC/BIC
df_age <- seq(1,6)
fit_age <- do.call(rbind,parLapply(cl, df_age, function (x,df_age,df,mvpa_df,lpa_df,mvpalpa_df) {
  do.call(rbind,lapply(df_age, function (y,x,df_age,df,mvpa_df,lpa_df,mvpalpa_df) {
    do.call(rbind,lapply(df_age, function (z,y,x,df) {
      fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = mvpa_df) + ns(LPA, df = lpa_df) + ns(MVPA*LPA, df = mvpalpa_df) + 
                     ns(logage, df = x) + ns(logage*MVPA, df = y) + ns(logage*LPA, df = z) + 
                     sex_male + I(sex_male*MVPA) + I(sex_male*LPA) + ethnicity_nonwhite + 
                     b_educ_vocational + b_educ_olevels + b_educ_alevels + b_educ_other + 
                     b_income_18to31k + b_income_31to52k + b_income_52to100k + b_income_100kplus + b_employ_other +
                     b_health_good + b_health_fair + b_health_poor + 
                     b_bmi_underweight + I(b_bmi_underweight*MVPA) + I(b_bmi_underweight*LPA) + 
                     b_bmi_overweight + I(b_bmi_overweight*MVPA) + I(b_bmi_overweight*LPA) + 
                     b_bmi_obese + I(b_bmi_obese*MVPA) + I(b_bmi_obese*LPA) + 
                     b_alc_ex + b_alc_current_infrequent + b_alc_current_weekly + 
                     b_smk_ex + I(b_smk_ex*MVPA) + I(b_smk_ex*LPA) + b_smk_current + 
                     I(b_smk_current*MVPA) + I(b_smk_current*LPA) + diet_score + 
                     seasonality_Spring + seasonality_Summer + seasonality_Winter + 
                     TDI_bsl + valid_weartime_indays + ns(b_met, df = 1) + b_sleep_cat_enough + b_sleep_cat_toomuch +
                     b_diabetes + b_hypertension + b_highcholesterol + b_affective_disorder,
                   data = data,
                   df = df)
      c(x,y,z,AIC(fit),BIC(fit))
    },x=x,y=y,df=df))
  },x=x,df_age=df_age,df=df,mvpa_df=mvpa_df,lpa_df=lpa_df,mvpalpa_df=mvpalpa_df))
},df_age=df_age,df=df,mvpa_df=mvpa_df,lpa_df=lpa_df,mvpalpa_df=mvpalpa_df))

best_age_aic <- arrayInd(which.min(fit_age[,4]), dim(fit_age))
best_age_bic <- arrayInd(which.min(fit_age[,5]), dim(fit_age))

# 6.2. Print the best fitting models based on AIC/BIC
paste0("The best fitting model based on AIC has df=",fit_age[best_age_aic[[1]],1]," for age, df=",fit_age[best_age_aic[[1]],2]," for age*MVPA, and df=",fit_age[best_age_aic[[1]],3]," for age*LPA")
paste0("The best fitting model based on BIC has df=",fit_age[best_age_bic[[1]],1]," for age, df=",fit_age[best_age_bic[[1]],2]," for age*MVPA, and df=",fit_age[best_age_bic[[1]],3]," for age*LPA")

age_df <- fit_age[best_age_aic[[1]],1]
agemvpa_df <- fit_age[best_age_aic[[1]],2]
agelpa_df <- fit_age[best_age_aic[[1]],3]

######################################################################################
# 7. Check splines for baseline mvpa
#-------------------------------------------------------------------------------------

# 7.1. Run all models and save AIC/BIC
df_b_mvpa <- seq(1,6)
fit_b_mvpa <- do.call(rbind,parLapply(cl, df_b_mvpa, function (x,df,mvpa_df,lpa_df,mvpalpa_df,age_df,agemvpa_df,agelpa_df) {
  
  fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = mvpa_df) + ns(LPA, df = lpa_df) + ns(MVPA*LPA, df = mvpalpa_df) + 
                 ns(logage, df = age_df) + ns(logage*MVPA, df = agemvpa_df) + ns(logage*LPA, df = agelpa_df) + 
                 sex_male + I(sex_male*MVPA) + I(sex_male*LPA) + ethnicity_nonwhite + 
                 b_educ_vocational + b_educ_olevels + b_educ_alevels + b_educ_other + 
                 b_income_18to31k + b_income_31to52k + b_income_52to100k + b_income_100kplus + b_employ_other +
                 b_health_good + b_health_fair + b_health_poor + 
                 b_bmi_underweight + I(b_bmi_underweight*MVPA) + I(b_bmi_underweight*LPA) + 
                 b_bmi_overweight + I(b_bmi_overweight*MVPA) + I(b_bmi_overweight*LPA) + 
                 b_bmi_obese + I(b_bmi_obese*MVPA) + I(b_bmi_obese*LPA) + 
                 b_alc_ex + b_alc_current_infrequent + b_alc_current_weekly + 
                 b_smk_ex + I(b_smk_ex*MVPA) + I(b_smk_ex*LPA) + b_smk_current + 
                 I(b_smk_current*MVPA) + I(b_smk_current*LPA) + diet_score + 
                 seasonality_Spring + seasonality_Summer + seasonality_Winter + 
                 TDI_bsl + valid_weartime_indays + ns(b_met, df = x) + b_sleep_cat_enough + b_sleep_cat_toomuch +
                 b_diabetes + b_hypertension + b_highcholesterol + b_affective_disorder,
               data = data,
               df = df)
  c(x,AIC(fit),BIC(fit))
  
},df=df,mvpa_df=mvpa_df,lpa_df=lpa_df,mvpalpa_df=mvpalpa_df,age_df=age_df,agemvpa_df=agemvpa_df,agelpa_df=agelpa_df))

best_b_mvpa_aic <- arrayInd(which.min(fit_b_mvpa[,2]), dim(fit_b_mvpa))
best_b_mvpa_bic <- arrayInd(which.min(fit_b_mvpa[,3]), dim(fit_b_mvpa))

# 6.2. Print the best fitting models based on AIC/BIC
paste0("The best fitting model based on AIC has df=",fit_b_mvpa[best_b_mvpa_aic[[1]],1])
paste0("The best fitting model based on BIC has df=",fit_b_mvpa[best_b_mvpa_bic[[1]],1])

b_mvpa_df <- fit_b_mvpa[best_b_mvpa_aic[[1]],1]

######################################################################################
# 7. Fit final complete case model for shiny app
#-------------------------------------------------------------------------------------

fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = mvpa_df) + ns(LPA, df = lpa_df) + ns(MVPA*LPA, df = mvpalpa_df) + 
               ns(logage, df = age_df) + ns(logage*MVPA, df = agemvpa_df) + ns(logage*LPA, df = agelpa_df) + 
               sex_male + I(sex_male*MVPA) + I(sex_male*LPA) + ethnicity_nonwhite + 
               b_educ_vocational + b_educ_olevels + b_educ_alevels + b_educ_other + 
               b_income_18to31k + b_income_31to52k + b_income_52to100k + b_income_100kplus + b_employ_other +
               b_health_good + b_health_fair + b_health_poor + 
               b_bmi_underweight + I(b_bmi_underweight*MVPA) + I(b_bmi_underweight*LPA) + 
               b_bmi_overweight + I(b_bmi_overweight*MVPA) + I(b_bmi_overweight*LPA) + 
               b_bmi_obese + I(b_bmi_obese*MVPA) + I(b_bmi_obese*LPA) + 
               b_alc_ex + b_alc_current_infrequent + b_alc_current_weekly + 
               b_smk_ex + I(b_smk_ex*MVPA) + I(b_smk_ex*LPA) + b_smk_current + 
               I(b_smk_current*MVPA) + I(b_smk_current*LPA) + diet_score + 
               seasonality_Spring + seasonality_Summer + seasonality_Winter + 
               TDI_bsl + valid_weartime_indays + ns(b_met, df = b_mvpa_df) + b_sleep_cat_enough + b_sleep_cat_toomuch +
               b_diabetes + b_hypertension + b_highcholesterol + b_affective_disorder,
             data = data,
             df = df)

######################################################################################
# 8. Save final model and model fit for appendices
#-------------------------------------------------------------------------------------

saveRDS(fit,paste0(workdir,"Data/model_fit.rds"))
saveRDS(fit_dist,paste0(workdir,"Results/Model fit - distribution.rds"))
saveRDS(fit_act,paste0(workdir,"Results/Model fit - activity splines.rds"))
saveRDS(fit_age,paste0(workdir,"Results/Model fit - age splines.rds"))

wb1 <- loadWorkbook(paste0(workdir,"Results/model_fit.xlsx"), create = TRUE)
createSheet(wb1, name = "survival")
createSheet(wb1, name = "activity")
createSheet(wb1, name = "age")
writeWorksheet(wb1,fit_dist,"survival",startRow = 1, startCol = 1, header = FALSE)
writeWorksheet(wb1,fit_act,"activity",startRow = 1, startCol = 1, header = FALSE)
writeWorksheet(wb1,fit_age,"age",startRow = 1, startCol = 1, header = FALSE)
saveWorkbook(wb1)

