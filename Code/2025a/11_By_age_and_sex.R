######################################################################################
##   
## Analysis of light physical activity and survival
## Estimate results by age group and sex
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
libs <- c("Amelia","arm","dplyr","ggplot2","ggpubr","msm","parallel","rstpm2","splines","survival","tidyr")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

######################################################################################
# 2. Define function to run model and get predictions in logit scale
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
               data = data)
}

get_results_sex <- function (fit) {
  
  newdata <- data.frame(fit@data)[,-36]
  
  mvpa <- quantile(newdata$MVPA,p=c(0.1,0.5,0.9))
  sex <- c(0,1)
  
  grid <- expand.grid(mvpa, sex, stringsAsFactors = FALSE)
  
  predict <- do.call(rbind,lapply(seq(1,nrow(grid)), function (x) {
    do.call(rbind,lapply(seq(0,600,15), function (y) {
      
      newdata$MVPA <- grid[x,1]
      newdata$LPA <- y
      newdata$sex_male <- grid[x,2]
      newdata$f_death <- 10
      
      p <- 1-predict(fit,type="meansurv",
                   newdata=newdata,
                   full=TRUE,
                   se.fit=TRUE)
      
      cbind(grid[x,],y,p)
      
    }))
  }))
  
  colnames(predict) <- c("MVPA","sex","LPA","time","est","lower","upper")
  predict$se <- (predict$upper-predict$est)/qnorm(0.975)
  predict$MVPA <- factor(predict$MVPA,
                         labels=c("Tenth","Median","Ninetieth"))
  
  p <- pivot_wider(predict[,c(1,2,3,5)],
                   names_from = c("sex","MVPA","LPA"),
                   names_sep = "_",
                   values_from = c("est"))
  se <- pivot_wider(predict[,c(1,2,3,8)],
                    names_from = c("sex","MVPA","LPA"),
                    names_sep = "_",
                    values_from = c("se"))
  
  list(p,se)
}

get_results_age <- function (fit) {
  
  data <- data.frame(fit@data)[,-36]
  
  mvpa <- quantile(data$MVPA,p=c(0.1,0.5,0.9))
  age <- c(1,2,3)
  
  grid <- expand.grid(mvpa, age, stringsAsFactors = FALSE)
  
  age_list <- list(log(seq(40,49,1)),
                   log(seq(50,59,1)),
                   log(seq(60,69,1)))
  
  predict <- do.call(rbind,lapply(seq(1,nrow(grid)), function (x) {
    do.call(rbind,lapply(seq(0,600,15), function (y) {
      
      newdata <- do.call(rbind,lapply(age_list[[grid[x,2]]],function (y) {
        z <- data
        z$logage <- y
        z
      }))
      newdata$MVPA <- grid[x,1]
      newdata$LPA <- y
      newdata$f_death <- 10
      
      p <- 1-predict(fit,type="meansurv",
                   newdata=newdata,
                   full=TRUE,
                   se.fit=TRUE)
      
      cbind(grid[x,],y,p)
      
    }))
  }))
  
  colnames(predict) <- c("MVPA","age","LPA","time","est","lower","upper")
  predict$se <- (predict$upper-predict$est)/qnorm(0.975)
  predict$MVPA <- factor(predict$MVPA,
                         labels=c("Tenth","Median","Ninetieth"))
  
  p <- pivot_wider(predict[,c(1,2,3,5)],
                   names_from = c("age","MVPA","LPA"),
                   names_sep = "_",
                   values_from = c("est"))
  se <- pivot_wider(predict[,c(1,2,3,8)],
                    names_from = c("age","MVPA","LPA"),
                    names_sep = "_",
                    values_from = c("se"))
  
  list(p,se)
}

meld_results <- function (res) {
  b <- do.call(rbind,lapply(res, function (x) {
    x[[1]]
  }))
  se <- do.call(rbind,lapply(res, function (x) {
    x[[2]]
  }))
  
  meld_data <- mi.meld(q=as.matrix(b),se=as.matrix(se))
}

construct_cis <- function (res) {
  res$conf.low <- res$est - qnorm(0.975)*res$se
  res$conf.high <- res$est + qnorm(0.975)*res$se
  
  res$logitest <- logit(res$est)
  res$logitse <- sapply(seq(1,nrow(res)), function (x) {
    se <- msm::deltamethod(g=~log(x1/(1-x1)),mean=res$est[x],cov=(res$se[x])^2,ses=TRUE)
  })
  res$logit.low <- invlogit(res$logitest - qnorm(0.975)*res$logitse)
  res$logit.high <- invlogit(res$logitest + qnorm(0.975)*res$logitse)
  
  res$logest <- log(res$est)
  res$logse <- sapply(seq(1,nrow(res)), function (x) {
    se <- msm::deltamethod(g=~log(x1),mean=res$est[x],cov=(res$se[x])^2,ses=TRUE)
  })
  res$log.low <- exp(res$logest - qnorm(0.975)*res$logse)
  res$log.high <- exp(res$logest + qnorm(0.975)*res$logse)
  res
}

######################################################################################
# 3. Load imputed data
#-------------------------------------------------------------------------------------

data <- readRDS(paste0(workdir,"Data/primary_analysis_data.rds"))

######################################################################################
# 4. Create cluster for parallel processing
#-------------------------------------------------------------------------------------

cl <- makeCluster(4)
clusterExport(cl,c("libs","data"))
clusterEvalQ(cl, lapply(libs, library, character.only = TRUE))

######################################################################################
# 5. Run model in each imputation and pool using Rubin's rules
#-------------------------------------------------------------------------------------

start <- Sys.time()
model_fit <- parLapply(cl, data, get_fit)
results_sex <- parLapply(cl, model_fit, get_results_sex)
results_age <- parLapply(cl, model_fit, get_results_age)
end <- Sys.time()
end-start

######################################################################################
# 6. Pool results over imputations
#-------------------------------------------------------------------------------------

# 6.1 Results by sex
sex_res <- meld_results(results_sex)

sex_res <- merge(pivot_longer(data.frame(sex_res[[1]]),
                               cols = everything(),
                               values_to = "est",
                               names_pattern = "X(.*)_(.*)_(.*)",
                               names_to = c("sex","MVPA","LPA")),
                  pivot_longer(data.frame(sex_res[[2]]),
                               cols = everything(),
                               values_to = "se",
                               names_pattern = "X(.*)_(.*)_(.*)",
                               names_to = c("sex","MVPA","LPA")),
                  by=c("sex","MVPA","LPA"))

sex_res$LPA <- as.numeric(sex_res$LPA)
sex_res$sex <- factor(sex_res$sex,
                      labels=c("Female","Male"))
sex_res$MVPA <- factor(sex_res$MVPA,
                       levels=c("Tenth","Median","Ninetieth"),
                       labels=c("10th percentile","Median","90th percentile"))

sex_res <- construct_cis(sex_res)

# 6.2 Results by age category
age_res <- meld_results(results_age)

age_res <- merge(pivot_longer(data.frame(age_res[[1]]),
                               cols = everything(),
                               values_to = "est",
                               names_pattern = "X(.*)_(.*)_(.*)",
                               names_to = c("age","MVPA","LPA")),
                  pivot_longer(data.frame(age_res[[2]]),
                               cols = everything(),
                               values_to = "se",
                               names_pattern = "X(.*)_(.*)_(.*)",
                               names_to = c("age","MVPA","LPA")),
                  by=c("age","MVPA","LPA"))

age_res$LPA <- as.numeric(age_res$LPA)
age_res$MVPA <- factor(age_res$MVPA,
                      levels=c("Tenth","Median","Ninetieth"),
                      labels=c("10th percentile","Median","90th percentile"))

age_res <- construct_cis(age_res)

######################################################################################
# 7. Define theme
#-------------------------------------------------------------------------------------

figure_theme <- theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
        axis.line = element_line(colour = 'grey80', linewidth = 0.3),
        axis.ticks = element_line(colour = "grey80", linewidth = 0.3),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(hjust = -0.01),
        legend.position="bottom")

######################################################################################
# 8. Create figures
#-------------------------------------------------------------------------------------

sex_figure <- ggplot(sex_res,aes(x=LPA, y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.35,0.05),limits=c(0.0,0.35)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  facet_wrap(sex ~ .,
             ncol=1,
             labeller=labeller(sex=c(Female="(a) Female",
                                         Male="(b) Male"))) +
  figure_theme

age_labs = c("(a) Ages 40-49","(b) Ages 50-59","(c) Ages 60-69")
names(age_labs) = c('1','2','3')

age_figure <- ggplot(age_res,aes(x=LPA, y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  facet_wrap(age ~ .,
             ncol=1,
             labeller=as_labeller(age_labs)) +
  figure_theme

sex_figure
age_figure

######################################################################################
# 9. Save each imputed analysis set
#-------------------------------------------------------------------------------------

saveRDS(results_sex,paste0(workdir,"Results/predictions by sex 20250302.rds"))
ggsave(paste0(workdir,"Results/figure by sex 20241121.png"),sex_figure)
saveRDS(results_age,paste0(workdir,"Results/predictions by age 20250302.rds"))
ggsave(paste0(workdir,"Results/figure by age 20241121.png"),age_figure)