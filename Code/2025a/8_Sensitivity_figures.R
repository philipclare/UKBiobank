######################################################################################
##   
## Analysis of light physical activity and survival
## Create figures for primary analysis
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
  workdir <- "C:/Users/pcla5984/The University of Sydney (Staff)/Susan Luo - data and R/"
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
# 2. Define functions to meld results and construct CIs
#-------------------------------------------------------------------------------------

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
# 3. Load analysis results
#-------------------------------------------------------------------------------------

s1_res <- readRDS(paste0(workdir,"Results/Model Results - S1 - Combined.rds"))
s2_res <- readRDS(paste0(workdir,"Results/Model Results - S2 - Combined.rds"))
s3_res <- readRDS(paste0(workdir,"Results/Model Results - S3 - Combined.rds"))
s4_res <- readRDS(paste0(workdir,"Results/Model Results - S4 - Combined.rds"))

######################################################################################
# 4. Pool results over imputations
#-------------------------------------------------------------------------------------

# 4.1. S1 analysis
s1_res <- meld_results(s1_res)

s1_res <- merge(pivot_longer(data.frame(s1_res[[1]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "est"),
                pivot_longer(data.frame(s1_res[[2]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "se"))

s1_res$LPA <- as.numeric(s1_res$LPA)
s1_res$MVPA <- factor(s1_res$MVPA,
                      levels=c("Tenth","Median","Ninetieth"),
                      labels=c("10th percentile","Median","90th percentile"))

s1_res <- construct_cis(s1_res)

# 4.2. S2 analysis
s2_res <- meld_results(s2_res)

s2_res <- merge(pivot_longer(data.frame(s2_res[[1]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "est"),
                pivot_longer(data.frame(s2_res[[2]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "se"))

s2_res$LPA <- as.numeric(s2_res$LPA)
s2_res$MVPA <- factor(s2_res$MVPA,
                      levels=c("Tenth","Median","Ninetieth"),
                      labels=c("10th percentile","Median","90th percentile"))

s2_res <- construct_cis(s2_res)

# 4.3. S3 analysis
s3_res <- meld_results(s3_res)

s3_res <- merge(pivot_longer(data.frame(s3_res[[1]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "est"),
                pivot_longer(data.frame(s3_res[[2]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "se"))

s3_res$LPA <- as.numeric(s3_res$LPA)
s3_res$MVPA <- factor(s3_res$MVPA,
                      levels=c("Tenth","Median","Ninetieth"),
                      labels=c("10th percentile","Median","90th percentile"))

s3_res <- construct_cis(s3_res)

# 4.4. S4 analysis
s4_res <- s4_res[[1]]

s4_res <- merge(pivot_longer(data.frame(s4_res[[1]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "est"),
                pivot_longer(data.frame(s4_res[[2]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "se"))

s4_res$LPA <- as.numeric(s4_res$LPA)
s4_res$MVPA <- factor(s4_res$MVPA,
                      levels=c("Tenth","Median","Ninetieth"),
                      labels=c("10th percentile","Median","90th percentile"))

s4_res <- construct_cis(s4_res)

######################################################################################
# 5. Define theme
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
# 6. Create figures
#-------------------------------------------------------------------------------------

# 6.1 S1
s1_figure <- ggplot(s1_res,aes(x=LPA,y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme

s1_figure

# 6.2 S2
s2_figure <- ggplot(s2_res,aes(x=LPA,y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme

s2_figure

# 6.3 S3
s3_figure <- ggplot(s3_res,aes(x=LPA,y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme

s3_figure

# 6.4 S4
s4_figure <- ggplot(s4_res,aes(x=LPA,y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme

s4_figure

######################################################################################
# 7. Save figures
#-------------------------------------------------------------------------------------

ggsave(paste0(workdir,"Results/figure S1 20250501.png"),s1_figure)
ggsave(paste0(workdir,"Results/figure S2 20250501.png"),s2_figure)
ggsave(paste0(workdir,"Results/figure S3 20250501.png"),s3_figure)
ggsave(paste0(workdir,"Results/figure S4 20250501.png"),s4_figure)
