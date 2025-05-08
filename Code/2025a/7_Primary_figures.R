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

results_lpa <- readRDS(paste0(workdir,"Results/Model Results - Primary - LPA.rds"))
results_mvpa <- readRDS(paste0(workdir,"Results/Model Results - Primary - MVPA.rds"))
results_both <- readRDS(paste0(workdir,"Results/Model Results - Primary - Combined.rds"))

######################################################################################
# 4. Pool results over imputations
#-------------------------------------------------------------------------------------

# 4.1. Main analysis
pr_res <- meld_results(results_both)

pr_res <- merge(pivot_longer(data.frame(pr_res[[1]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "est"),
                pivot_longer(data.frame(pr_res[[2]]),
                             cols = everything(),
                             names_to = c("MVPA","LPA"),
                             names_pattern = "(.*)_(.*)",
                             values_to = "se"))

pr_res$LPA <- as.numeric(pr_res$LPA)
pr_res$MVPA <- factor(pr_res$MVPA,
                      levels=c("Tenth","Median","Ninetieth"),
                      labels=c("10th percentile","Median","90th percentile"))

pr_res <- construct_cis(pr_res)

# 4.2. LPA analysis
lpa_res <- meld_results(results_lpa)

lpa_res <- merge(pivot_longer(data.frame(lpa_res[[1]]),
                              cols = everything(),
                              names_to = c("LPA"),
                              names_pattern = "X(.*)",
                              values_to = "est"),
                 pivot_longer(data.frame(lpa_res[[2]]),
                              cols = everything(),
                              names_to = c("LPA"),
                              names_pattern = "X(.*)",
                              values_to = "se"))

lpa_res$LPA <- as.numeric(lpa_res$LPA)

lpa_res <- construct_cis(lpa_res)

# 4.3. MVPA analysis
mvpa_res <- meld_results(results_mvpa)

mvpa_res <- merge(pivot_longer(data.frame(mvpa_res[[1]]),
                               cols = everything(),
                               names_to = c("MVPA"),
                               names_pattern = "X(.*)",
                               values_to = "est"),
                  pivot_longer(data.frame(mvpa_res[[2]]),
                               cols = everything(),
                               names_to = c("MVPA"),
                               names_pattern = "X(.*)",
                               values_to = "se"))

mvpa_res$MVPA <- as.numeric(mvpa_res$MVPA)

mvpa_res <- construct_cis(mvpa_res)

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

primary_figure <- ggplot(pr_res,aes(x=LPA,y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme

lpa_figure <- ggplot(lpa_res,aes(x=LPA,y=est)) + 
  geom_line(colour="#619CFF") + 
  geom_ribbon(aes(ymin=logit.low,ymax=logit.high),alpha=0.3,colour=NA,fill="#619CFF") +
  ylab("Probability of death (%)") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
mvpa_figure <- ggplot(mvpa_res,aes(x=MVPA,y=est)) + 
  geom_line(colour="#00BA38") + 
  geom_ribbon(aes(ymin=log.low,ymax=log.high),alpha=0.3,colour=NA,fill="#00BA38") +
  ylab("Probability of death (%)") +
  xlab("Minutes of moderate/vigorous activity per day") + 
  scale_x_continuous(breaks=seq(0,100,10)) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0.0,0.3,0.05),limits=c(0.0,0.3)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
secondary_figure <- ggarrange(lpa_figure,mvpa_figure,nrow=2)

primary_figure
secondary_figure

######################################################################################
# 7. Save figures
#-------------------------------------------------------------------------------------

ggsave(paste0(workdir,"Results/figure 1 20250501.png"),secondary_figure)
ggsave(paste0(workdir,"Results/figure 2 20250501.png"),primary_figure)
