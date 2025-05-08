######################################################################################
##   
## Analysis of light physical activity and survival
## Compare specific counterfactual scenarios
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
libs <- c("Amelia","arm","dplyr","flexsurv","ggplot2","marginaleffects","msm","parallel","splines","survival","tidyr","rstpm2","XLConnect")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

######################################################################################
# 2. Define functions to run model and get predictions in logit scale
#-------------------------------------------------------------------------------------

get_diffs <- function (data) {
  
  comp_table <- matrix(c(22,0,11,22,44,60,60,60,60,60,0,0,11,22,44,300,300,300,300,300),nrow=5)
  
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
               data = data,
               df = 3)
  
  predict <- do.call(rbind,lapply(seq(1,5), function (x){
    mvpa_ref <- comp_table[x,1]
    lpa_ref <- comp_table[x,2]
    mvpa_cmp <- comp_table[x,3]
    lpa_cmp <- comp_table[x,4]
    
    diff <- predict(fit,type="meansurvdiff",
                    newdata=data.frame(data[,!(names(data) %in% c("MVPA","LPA","f_death"))],MVPA=mvpa_ref, LPA=lpa_ref,f_death=10),
                    exposed=function(data) transform(data, MVPA=mvpa_cmp, LPA=lpa_cmp),
                    var=c("MVPA","LPA"),
                    full=TRUE,
                    se.fit=TRUE)
    
    newdata <- data.frame(fit@data)
    newdata$f_death <- 10
    
    predict <- avg_predictions(fit,
                               newdata=newdata,
                               variables=list(LPA=c(lpa_ref,lpa_cmp),
                                              MVPA=c(mvpa_ref,mvpa_cmp)),
                               vcov=TRUE,
                               type="fail")
    
    ref_index <- which(predict$LPA==lpa_ref & predict$MVPA==mvpa_ref)
    var <- vcov(predict)
    
    predict$rd <- predict$estimate-predict$estimate[ref_index]
    predict$rd.se <- do.call(rbind,lapply(seq(1,nrow(predict)), function (x,var,ref_index) {
      m <- predict$estimate[c(x,ref_index)]
      newvar <- var[c(x,ref_index),c(x,ref_index)]
      msm::deltamethod(~x1-x2,m,newvar)
    },var=var,ref_index=ref_index))
    
    predict$rr <- predict$estimate/predict[ref_index,]$estimate
    predict$rr.se <- c(0,do.call(rbind,lapply(seq(2,nrow(predict)), function (x,var,ref_index) {
      m <- predict$estimate[c(x,ref_index)]
      newvar <- var[c(x,ref_index),c(x,ref_index)]
      msm::deltamethod(~x1/x2,m,newvar)
    },var=var,ref_index=ref_index)))
    
    comp_index <- which(predict$LPA==lpa_cmp & predict$MVPA==mvpa_cmp)
    
    c(mvpa_ref,lpa_ref,mvpa_cmp,lpa_cmp,diff$Estimate,(diff$upper-diff$Estimate)/qnorm(0.975),predict$rd[comp_index],predict$rd.se[comp_index],predict$rr[comp_index],predict$rr.se[comp_index])
  }))

  colnames(predict) <- c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp","est","se","rd","rd.se","rr","rr.se")

  p <- pivot_wider(data.frame(predict[,c(1,2,3,4,5)]),
                   names_from = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                   names_sep = "_",
                   values_from = c("est"))
  se <- pivot_wider(data.frame(predict[,c(1,2,3,4,6)]),
                    names_from = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                    names_sep = "_",
                    values_from = c("se"))
  
  rd <- pivot_wider(data.frame(predict[,c(1,2,3,4,7)]),
                   names_from = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                   names_sep = "_",
                   values_from = c("rd"))
  rd.se <- pivot_wider(data.frame(predict[,c(1,2,3,4,8)]),
                    names_from = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                    names_sep = "_",
                    values_from = c("rd.se"))
  
  rr <- pivot_wider(data.frame(predict[,c(1,2,3,4,9)]),
                   names_from = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                   names_sep = "_",
                   values_from = c("rr"))
  rr.se <- pivot_wider(data.frame(predict[,c(1,2,3,4,10)]),
                    names_from = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                    names_sep = "_",
                    values_from = c("rr.se"))
  
  list(p,se,rd,rd.se,rr,rr.se)
}

######################################################################################
# 3. Load imputed data
#-------------------------------------------------------------------------------------

data <- readRDS(paste0(workdir,"Data/primary_analysis_data.rds"))

######################################################################################
# 4. Create cluster for parallel processing
#-------------------------------------------------------------------------------------

cl <- makeCluster(2)
clusterExport(cl,c("libs","data"))
clusterEvalQ(cl, lapply(libs, library, character.only = TRUE))

######################################################################################
# 4. Run model in each imputation and pool using Rubin's rules
#-------------------------------------------------------------------------------------

start <- Sys.time()
diffs <- parLapply(cl, data, get_diffs)
end <- Sys.time()
end-start

######################################################################################
# 5. Pool results over imputations
#-------------------------------------------------------------------------------------

b <- do.call(rbind,lapply(diffs, function (x) {
  x[[1]]
}))
se <- do.call(rbind,lapply(diffs, function (x) {
  x[[1]]
}))

rd <- do.call(rbind,lapply(diffs, function (x) {
  x[[3]]
}))
rd.se <- do.call(rbind,lapply(diffs, function (x) {
  x[[4]]
}))

rr <- do.call(rbind,lapply(diffs, function (x) {
  x[[5]]
}))
rr.se <- do.call(rbind,lapply(diffs, function (x) {
  x[[6]]
}))

b_meld_data <- mi.meld(q=as.matrix(b),se=as.matrix(se))
rd_meld_data <- mi.meld(q=as.matrix(rd),se=as.matrix(rd.se))
rr_meld_data <- mi.meld(q=as.matrix(rr),se=as.matrix(rr.se))


tab_data <- merge(merge(merge(pivot_longer(data.frame(rd_meld_data[[1]]),
                               cols = everything(),
                               names_to = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                               names_pattern = "X(.*)_(.*)_(.*)_(.*)",
                               values_to = "rd"),
                  pivot_longer(data.frame(rd_meld_data[[2]]),
                               cols = everything(),
                               names_to = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                               names_pattern = "X(.*)_(.*)_(.*)_(.*)",
                               values_to = "rd.se")),
                  pivot_longer(data.frame(rr_meld_data[[1]]),
                               cols = everything(),
                               names_to = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                               names_pattern = "X(.*)_(.*)_(.*)_(.*)",
                               values_to = "rr")),
                  pivot_longer(data.frame(rr_meld_data[[2]]),
                               cols = everything(),
                               names_to = c("MVPA_ref","LPA_ref","MVPA_cmp","LPA_cmp"),
                               names_pattern = "X(.*)_(.*)_(.*)_(.*)",
                               values_to = "rr.se"))

tab_data$rd.conf.low <- tab_data$rd - qnorm(0.975)*tab_data$rd.se
tab_data$rd.conf.high <- tab_data$rd + qnorm(0.975)*tab_data$rd.se
tab_data$rr.conf.low <- tab_data$rr - qnorm(0.975)*tab_data$rr.se
tab_data$rr.conf.high <- tab_data$rr + qnorm(0.975)*tab_data$rr.se

View(tab_data)

wb1 <- loadWorkbook(paste0(workdir,"Results/scenario_comparisons 20250501.xlsx"), create = TRUE)
createSheet(wb1, name = "Comp")
writeWorksheet(wb1,tab_data,"Comp",startRow = 1, startCol = 1, header = TRUE)
saveWorkbook(wb1)