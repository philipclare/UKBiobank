######################################################################################
##   
## Analysis of light physical activity and survival
## Run best fitting model and generate relative risks and risk differences
## By: Philip Clare
## Date: 2/12/2024
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
  workdir <- "C:/Users/pcla5984/The University of Sydney (Staff)/Susan Luo - data and R/"
} else if (Sys.info()[['sysname']]=="Darwin") {
  workdir <- "/Users/pjclare/Library/CloudStorage/OneDrive-SharedLibraries-TheUniversityofSydney(Staff)/Susan Luo - data and R/" # MAC
}

# 1.2. Define packages to be used, check if installed, and load
libs <- c("Amelia","arm","dplyr","flexsurv","ggplot2","ggpubr","marginaleffects","msm","parallel","splines","survival","tidyr","rstpm2")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

# 1.3. Set argument either for testing or passed by scheduler on HPC
if (Sys.info()[['sysname']]=="Linux") {
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
} else if (Sys.info()[['sysname']]=="Windows" | Sys.info()[['sysname']]=="Darwin") {
  args <- 1 # for testing on local computer
}

######################################################################################
# 2. Define functions to run model and get predictions in logit scale
#-------------------------------------------------------------------------------------

get_fit <- function (data) {
  
  fit <- stpm2(Surv(f_death, death) ~ ns(MVPA, df = 2) + ns(LPA, df = 2) + ns(MVPA*LPA, df = 4) + 
                 ns(logage, df = 1) + ns(logage*MVPA, df = 2) + ns(logage*LPA, df = 3) + 
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
               data = data)
}

diff_lpa <- function (fit) {
  
  newdata <- data.frame(fit@data)
  newdata$f_death <- 10
  quantile <- quantile(newdata$MVPA,p=0.001)
  
  predict <- avg_predictions(fit,
                             newdata=newdata,
                             variables=list(LPA=seq(0,600,15),
                                            MVPA=quantile),
                             vcov=TRUE,
                             type="fail")
  
  ref_index <- which(predict$LPA==60)
  var <- vcov(predict)
  
  predict$rd <- predict$estimate-predict$estimate[ref_index]
  predict$rd.se <- do.call(rbind,lapply(seq(1,nrow(predict)), function (x,var,ref_index) {
    m <- predict$estimate[c(x,ref_index)]
    newvar <- var[c(x,ref_index),c(x,ref_index)]
    msm::deltamethod(~x1-x2,m,newvar)
  },var=var,ref_index=ref_index))
  
  predict$rr <- predict$estimate/predict[ref_index,]$estimate
  predict$rr.se <- c(do.call(rbind,lapply(seq(1,nrow(predict)), function (x,var,ref_index) {
    m <- predict$estimate[c(x,ref_index)]
    newvar <- var[c(x,ref_index),c(x,ref_index)]
    msm::deltamethod(~x1/x2,m,newvar)
  },var=var,ref_index=ref_index)))
  
  predict
  
}

diff_mvpa <- function (fit) {

  newdata <- data.frame(fit@data)
  newdata$f_death <- 10
  quantile <- quantile(newdata$LPA,p=0.001)
  
  predict <- avg_predictions(fit,
                             newdata=newdata,
                             variables=list(LPA=quantile,
                                            MVPA=seq(0,100,2)),
                             vcov=TRUE,
                             type="fail")
  
  ref_index <- which(predict$MVPA==0)
  var <- vcov(predict)
  
  predict$rd <- predict$estimate-predict$estimate[ref_index]
  predict$rd.se <- do.call(rbind,lapply(seq(1,nrow(predict)), function (x,var,ref_index) {
    m <- predict$estimate[c(x,ref_index)]
    newvar <- var[c(x,ref_index),c(x,ref_index)]
    msm::deltamethod(~x1-x2,m,newvar)
  },var=var,ref_index=ref_index))
  
  predict$rr <- predict$estimate/predict[ref_index,]$estimate
  predict$rr.se <- c(do.call(rbind,lapply(seq(1,nrow(predict)), function (x,var,ref_index) {
    m <- predict$estimate[c(x,ref_index)]
    newvar <- var[c(x,ref_index),c(x,ref_index)]
    msm::deltamethod(~x1/x2,m,newvar)
  },var=var,ref_index=ref_index)))
  
  predict
  
}

diff_both <- function (fit) {
  
  newdata <- data.frame(fit@data)
  newdata$f_death <- 10
  quantiles <- quantile(newdata$MVPA,p=c(0.1,0.5,0.9))
  
  p <- lapply(quantiles, function (x) {
    predict <- avg_predictions(fit,
                               newdata=newdata,
                               variables=list(LPA=seq(0,600,15),
                                              MVPA=x),
                               vcov=TRUE,
                               type="fail")
    
    ref_index <- which(predict$LPA==60)
    var <- vcov(predict)
    
    predict$rd <- predict$estimate-predict$estimate[ref_index]
    predict$rd.se <- do.call(rbind,lapply(seq(1,nrow(predict)), function (x,var,ref_index) {
      m <- predict$estimate[c(x,ref_index)]
      newvar <- var[c(x,ref_index),c(x,ref_index)]
      msm::deltamethod(~x1-x2,m,newvar)
    },var=var,ref_index=ref_index))
    
    predict$rr <- predict$estimate/predict[ref_index,]$estimate
    predict$rr.se <- c(do.call(rbind,lapply(seq(1,nrow(predict)), function (x,var,ref_index) {
      m <- predict$estimate[c(x,ref_index)]
      newvar <- var[c(x,ref_index),c(x,ref_index)]
      msm::deltamethod(~x1/x2,m,newvar)
    },var=var,ref_index=ref_index)))
    
    predict
    
  })
  
}

######################################################################################
# 3. Load imputed data
#-------------------------------------------------------------------------------------

data <- readRDS(paste0(workdir,"Data/primary_analysis_data.rds"))

data <- data[[args[1]]]

######################################################################################
# 4. Run model in each imputation and pool using Rubin's rules
#-------------------------------------------------------------------------------------

start <- Sys.time()
model_fit <- get_fit(data)
lpa <- diff_lpa(model_fit)
mvpa <- diff_mvpa(model_fit)
both <- diff_both(model_fit)
end <- Sys.time()
end-start

######################################################################################
# 5. Extract RRs and RDs
#-------------------------------------------------------------------------------------

b_lpa_pr <- pivot_wider(lpa[,c("LPA","estimate")],
                        names_from = c("LPA"),
                        names_sep = "_",
                        values_from = c("estimate"))
se_lpa_pr <- pivot_wider(lpa[,c("LPA","std.error")],
                         names_from = c("LPA"),
                         names_sep = "_",
                         values_from = c("std.error"))

b_lpa_rr <- pivot_wider(lpa[,c("LPA","rr")],
                        names_from = c("LPA"),
                        names_sep = "_",
                        values_from = c("rr"))
se_lpa_rr <- pivot_wider(lpa[,c("LPA","rr.se")],
                         names_from = c("LPA"),
                         names_sep = "_",
                         values_from = c("rr.se"))

b_lpa_rd <- pivot_wider(lpa[,c("LPA","rd")],
                        names_from = c("LPA"),
                        names_sep = "_",
                        values_from = c("rd"))
se_lpa_rd <- pivot_wider(lpa[,c("LPA","rd.se")],
                         names_from = c("LPA"),
                         names_sep = "_",
                         values_from = c("rd.se"))

b_mvpa_pr <- pivot_wider(mvpa[,c("MVPA","estimate")],
                        names_from = c("MVPA"),
                        names_sep = "_",
                        values_from = c("estimate"))
se_mvpa_pr <- pivot_wider(mvpa[,c("MVPA","std.error")],
                         names_from = c("MVPA"),
                         names_sep = "_",
                         values_from = c("std.error"))

b_mvpa_rr <- pivot_wider(mvpa[,c("MVPA","rr")],
                        names_from = c("MVPA"),
                        names_sep = "_",
                        values_from = c("rr"))
se_mvpa_rr <- pivot_wider(mvpa[,c("MVPA","rr.se")],
                         names_from = c("MVPA"),
                         names_sep = "_",
                         values_from = c("rr.se"))

b_mvpa_rd <- pivot_wider(mvpa[,c("MVPA","rd")],
                        names_from = c("MVPA"),
                        names_sep = "_",
                        values_from = c("rd"))
se_mvpa_rd <- pivot_wider(mvpa[,c("MVPA","rd.se")],
                         names_from = c("MVPA"),
                         names_sep = "_",
                         values_from = c("rd.se"))

b_both_pr <- lapply(seq(1,3), function (y) {
  dat <- both[[y]][,c("LPA","estimate")]
  pivot_wider(dat,
              names_from = c("LPA"),
              names_sep = "_",
              values_from = c("estimate"))
})
se_both_pr <- lapply(seq(1,3), function (y) {
  dat <- both[[y]][,c("LPA","std.error")]
  pivot_wider(dat,
              names_from = c("LPA"),
              names_sep = "_",
              values_from = c("std.error"))
})

b_both_rr <- lapply(seq(1,3), function (y) {
  dat <- both[[y]][,c("LPA","rr")]
  pivot_wider(dat,
              names_from = c("LPA"),
              names_sep = "_",
              values_from = c("rr"))
})
se_both_rr <- lapply(seq(1,3), function (y) {
  dat <- both[[y]][,c("LPA","rr.se")]
  pivot_wider(dat,
              names_from = c("LPA"),
              names_sep = "_",
              values_from = c("rr.se"))
})

b_both_rd <- lapply(seq(1,3), function (y) {
  dat <- both[[y]][,c("LPA","rd")]
  pivot_wider(dat,
              names_from = c("LPA"),
              names_sep = "_",
              values_from = c("rd"))
})
se_both_rd <- lapply(seq(1,3), function (y) {
  dat <- both[[y]][,c("LPA","rd.se")]
  pivot_wider(dat,
              names_from = c("LPA"),
              names_sep = "_",
              values_from = c("rd.se"))
})

res_lpa <- list(b_lpa_pr,se_lpa_pr,b_lpa_rr,se_lpa_rr,b_lpa_rd,se_lpa_rd)
res_mvpa <- list(b_mvpa_pr,se_mvpa_pr,b_mvpa_rr,se_mvpa_rr,b_mvpa_rd,se_mvpa_rd)
res_both <- list(b_both_pr,se_both_pr,b_both_rr,se_both_rr,b_both_rd,se_both_rd)

######################################################################################
# 5. Save results
#-------------------------------------------------------------------------------------

saveRDS(res_lpa,paste0(workdir,"Results/lpa relative ",args[[1]],".rds"))
saveRDS(res_mvpa,paste0(workdir,"Results/mvpa relative ",args[[1]],".rds"))
saveRDS(res_both,paste0(workdir,"Results/both relative ",args[[1]],".rds"))