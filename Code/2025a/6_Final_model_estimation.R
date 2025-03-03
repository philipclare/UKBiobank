######################################################################################
##   
## Analysis of light physical activity and survival
## Run best fitting model and save model fit object
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
libs <- c("Amelia","arm","dplyr","ggplot2","ggpubr","msm","parallel","rstpm2","splines","survival","tidyr")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

######################################################################################
# 2. Define functions to run model and get predictions in logit scale
#-------------------------------------------------------------------------------------

get_fit <- function (data) {
  
  fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = 2) + ns(LPA, df = 2) + ns(MVPA*LPA, df = 4) + 
                 ns(logage, df = 1) + ns(logage*MVPA, df = 2) + ns(logage*LPA, df = 1) + 
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
                 TDI_bsl + valid_weartime_indays,
               data = data,
               df = 3)
}

get_fit2 <- function (data) {
  
  fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = 2) + ns(LPA, df = 2) + ns(MVPA*LPA, df = 4) + 
                 ns(logage, df = 1) + ns(logage*MVPA, df = 2) + ns(logage*LPA, df = 1) + 
                 sex_male + I(sex_male*MVPA) + I(sex_male*LPA) + ethnicity_nonwhite + 
                 b_educ_vocational + b_educ_olevels + b_educ_alevels + b_educ_other + 
                 b_income_18to31k + b_income_31to52k + b_income_52to100k + b_income_100kplus + b_employ_other +
                 b_health_good + b_health_fair +
                 b_bmi_underweight + I(b_bmi_underweight*MVPA) + I(b_bmi_underweight*LPA) +
                 b_bmi_overweight + I(b_bmi_overweight*MVPA) + I(b_bmi_overweight*LPA) +
                 b_bmi_obese + I(b_bmi_obese*MVPA) + I(b_bmi_obese*LPA) +
                 b_alc_ex + b_alc_current_infrequent + b_alc_current_weekly + 
                 b_smk_ex + I(b_smk_ex*MVPA) + I(b_smk_ex*LPA) + b_smk_current +
                 I(b_smk_current*MVPA) + I(b_smk_current*LPA) + diet_score +
                 seasonality_Spring + seasonality_Summer + seasonality_Winter +
                 TDI_bsl + valid_weartime_indays,
               data = data,
               df = 3)
}

get_results_lpa <- function (fit) {
  
  newdata <- data.frame(fit@data)[,-36]
  predict <- do.call(rbind,lapply(seq(0,600,15), function (y) {
    newdata$LPA <- y
    newdata$MVPA <- quantile(fit@data$MVPA,0.001)
    newdata$f_death <- 10
    p <- 1-predict(fit,type="meansurv",
                 newdata=newdata,
                 full=TRUE,
                 se.fit=TRUE)
    cbind(y,p)
  }))
  
  colnames(predict) <- c("LPA","time","est","lower","upper")
  predict$se <- (predict$upper-predict$est)/qnorm(0.975)
  
  p <- pivot_wider(predict[,c(1,3)],
                   names_from = c("LPA"),
                   values_from = c("est"))
  se <- pivot_wider(predict[,c(1,6)],
                    names_from = c("LPA"),
                    values_from = c("se"))
  
  list(p,se)
}

get_results_mvpa <- function (fit) {
  
  newdata <- data.frame(fit@data)[,-36]
  predict <- do.call(rbind,lapply(seq(0,100,2), function (x) {
    newdata$LPA <- quantile(fit@data$LPA,0.001)
    newdata$MVPA <- x
    newdata$f_death <- 10
    p <- 1-predict(fit,type="meansurv",
                 newdata=newdata,
                 full=TRUE,
                 se.fit=TRUE)
    cbind(x,p)
  }))
  
  colnames(predict) <- c("MVPA","time","est","lower","upper")
  predict$se <- (predict$upper-predict$est)/qnorm(0.975)
  
  p <- pivot_wider(predict[,c(1,3)],
                   names_from = c("MVPA"),
                   values_from = c("est"))
  se <- pivot_wider(predict[,c(1,6)],
                    names_from = c("MVPA"),
                    values_from = c("se"))
  
  list(p,se)
}

get_results_both <- function (fit) {
  
  newdata <- data.frame(fit@data)[,-36]
  quantiles <- quantile(newdata$MVPA,p=c(0.1,0.5,0.9))
  
  predict <- do.call(rbind,lapply(quantiles, function (x) {
    do.call(rbind,lapply(seq(0,600,15), function (y) {
      newdata$MVPA <- x
      newdata$LPA <- y
      newdata$f_death <- 10
      p <- 1-predict(fit,type="meansurv",
                   newdata=newdata,
                   full=TRUE,
                   se.fit=TRUE)
      cbind(x,y,p)
    }))
  }))
  colnames(predict) <- c("MVPA","LPA","time","est","lower","upper")
  predict$se <- (predict$upper-predict$est)/qnorm(0.975)
  predict$MVPA <- factor(predict$MVPA,
                         labels=c("Tenth","Median","Ninetieth"))
  
  p <- pivot_wider(predict[,c(1,2,4)],
                   names_from = c("MVPA","LPA"),
                   names_sep = "_",
                   values_from = c("est"))
  se <- pivot_wider(predict[,c(1,2,7)],
                    names_from = c("MVPA","LPA"),
                    names_sep = "_",
                    values_from = c("se"))
  
  list(p,se)
}

######################################################################################
# 3. Load imputed data
#-------------------------------------------------------------------------------------

primary_data <- readRDS(paste0(workdir,"Data/primary_analysis_data.rds"))
s1_data <- readRDS(paste0(workdir,"Data/s1_analysis_data.rds"))
s2_data <- readRDS(paste0(workdir,"Data/s2_analysis_data.rds"))
s3_data <- readRDS(paste0(workdir,"Data/s3_analysis_data.rds"))
s4_data <- readRDS(paste0(workdir,"Data/s4_analysis_data.rds"))

######################################################################################
# 4. Create cluster for parallel processing
#-------------------------------------------------------------------------------------

cl <- makeCluster(2)
clusterExport(cl,c("libs","data"))
clusterEvalQ(cl, lapply(libs, library, character.only = TRUE))

######################################################################################
# 5. Run model in each imputation and pool using Rubin's rules
#-------------------------------------------------------------------------------------

# 5.1 Primary analysis
start1 <- Sys.time()
pr_model_fit <- parLapply(cl,primary_data, get_fit)
pr_results_lpa <- parLapply(cl, pr_model_fit, get_results_lpa)
pr_results_mvpa <- parLapply(cl,pr_model_fit, get_results_mvpa)
pr_results_both <- parLapply(cl, pr_model_fit, get_results_both)
end1 <- Sys.time()
end1-start1

# 5.2 First sensitivity analysis
start2 <- Sys.time()
s1_model_fit <- parLapply(cl,s1_data, get_fit)
s1_results_both <- parLapply(cl, s1_model_fit, get_results_both)
end2 <- Sys.time()
end2-start2

# 5.3 Second sensitivity analysis
start3 <- Sys.time()
s2_model_fit <- parLapply(cl, s2_data, get_fit2)
s2_results_both <- parLapply(cl, s2_model_fit, get_results_both)
end3 <- Sys.time()
end3-start3

# 5.4 Third sensitivity analysis
start4 <- Sys.time()
s3_model_fit <- parLapply(cl, s3_data, get_fit2)
s3_results_both <- parLapply(cl, s3_model_fit, get_results_both)
end4 <- Sys.time()
end4-start4

# 5.5 Fourth sensitivity analysis
start5 <- Sys.time()
s4_model_fit <- parLapply(cl, s4_data, get_fit)
s4_results_both <- parLapply(cl, s4_model_fit, get_results_both)
end5 <- Sys.time()
end5-start5

######################################################################################
# 6. Save results
#-------------------------------------------------------------------------------------

saveRDS(pr_results_lpa,paste0(workdir,"Results/Model Results - Primary - LPA.rds"))
saveRDS(pr_results_mvpa,paste0(workdir,"Results/Model Results - Primary - MVPA.rds"))
saveRDS(pr_results_both,paste0(workdir,"Results/Model Results - Primary - Combined.rds"))

saveRDS(s1_results_both,paste0(workdir,"Results/Model Results - S1 - Combined.rds"))
saveRDS(s2_results_both,paste0(workdir,"Results/Model Results - S2 - Combined.rds"))
saveRDS(s3_results_both,paste0(workdir,"Results/Model Results - S3 - Combined.rds"))
saveRDS(s4_results_both,paste0(workdir,"Results/Model Results - S4 - Combined.rds"))