######################################################################################
##   
## Analysis of light physical activity and survival
## Create relative risk and risk difference figures
## By: Philip Clare
## Date: 3/12/2024
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
libs <- c("Amelia","arm","dplyr","flexsurv","ggplot2","ggpubr","marginaleffects","msm","parallel","splines","survival","tidyr","rstpm2")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)


######################################################################################
# 2. Define functions to run model and get predictions in logit scale
#-------------------------------------------------------------------------------------

lpa <- lapply(seq(1,20), function (x) {
  readRDS(paste0(workdir,"Results/Artemis Output/lpa relative ",x,".rds"))
})

mvpa <- lapply(seq(1,20), function (x) {
  readRDS(paste0(workdir,"Results/Artemis Output/mvpa relative ",x,".rds"))
})

both <- lapply(seq(1,20), function (x) {
  readRDS(paste0(workdir,"Results/Artemis Output/both relative ",x,".rds"))
})

######################################################################################
# 5. Pool results over imputations
#-------------------------------------------------------------------------------------

b_lpa_rr <- do.call(rbind,lapply(lpa, function (x,y) {
    dat <- x[[3]]
  }))
se_lpa_rr <- do.call(rbind,lapply(lpa, function (x,y) {
    dat <- x[[4]]
  }))

b_lpa_rd <- do.call(rbind,lapply(lpa, function (x,y) {
    dat <- x[[5]]
  }))
se_lpa_rd <- do.call(rbind,lapply(lpa, function (x,y) {
    dat <- x[[6]]
  }))

rr_lpa <- t(do.call(rbind,mi.meld(q=as.matrix(b_lpa_rr),se=as.matrix(se_lpa_rr))))
rd_lpa <- t(do.call(rbind,mi.meld(q=as.matrix(b_lpa_rd),se=as.matrix(se_lpa_rd))))
res_lpa <- cbind(colnames(b_lpa_rr),rr_lpa,rd_lpa)
res_lpa <- data.frame(apply(res_lpa,2,as.numeric))
colnames(res_lpa) <- c("LPA","rr_est","rr_se","rd_est","rd_se")
res_lpa <- res_lpa %>% pivot_longer(c = 2:5,
                            names_to = c("type",".value"),
                            names_pattern = "(.+)_(.+)")

res_lpa$type <- factor(res_lpa$type)

b_mvpa_rr <- do.call(rbind,lapply(mvpa, function (x,y) {
  dat <- x[[3]]
}))
se_mvpa_rr <- do.call(rbind,lapply(mvpa, function (x,y) {
  dat <- x[[4]]
}))

b_mvpa_rd <- do.call(rbind,lapply(mvpa, function (x,y) {
  dat <- x[[5]]
}))
se_mvpa_rd <- do.call(rbind,lapply(mvpa, function (x,y) {
  dat <- x[[6]]
}))

rr_mvpa <- t(do.call(rbind,mi.meld(q=as.matrix(b_mvpa_rr),se=as.matrix(se_mvpa_rr))))
rd_mvpa <- t(do.call(rbind,mi.meld(q=as.matrix(b_mvpa_rd),se=as.matrix(se_mvpa_rd))))
res_mvpa <- cbind(colnames(b_mvpa_rr),rr_mvpa,rd_mvpa)
res_mvpa <- data.frame(apply(res_mvpa,2,as.numeric))
colnames(res_mvpa) <- c("MVPA","rr_est","rr_se","rd_est","rd_se")
res_mvpa <- res_mvpa %>% pivot_longer(c = 2:5,
                                      names_to = c("type",".value"),
                                      names_pattern = "(.+)_(.+)")

res_mvpa$type <- factor(res_mvpa$type)

b_both_rr <- lapply(seq(1,3), function (y) {
  do.call(rbind,lapply(both, function (x,y) {
    dat <- x[[3]][[y]]
  },y=y))
})
se_both_rr <- lapply(seq(1,3), function (y) {
  do.call(rbind,lapply(both, function (x,y) {
    dat <- x[[4]][[y]]
  },y=y))
})

b_both_rd <- lapply(seq(1,3), function (y) {
  do.call(rbind,lapply(both, function (x,y) {
    dat <- x[[5]][[y]]
  },y=y))
})
se_both_rd <- lapply(seq(1,3), function (y) {
  do.call(rbind,lapply(both, function (x,y) {
    dat <- x[[6]][[y]]
  },y=y))
})

meld_data <- do.call(rbind,lapply(seq(1,3), function (y) {
  rr <- t(do.call(rbind,mi.meld(q=as.matrix(b_both_rr[[y]]),se=as.matrix(se_both_rr[[y]]))))
  rd <- t(do.call(rbind,mi.meld(q=as.matrix(b_both_rd[[y]]),se=as.matrix(se_both_rd[[y]]))))
  res <- cbind(rep(y,nrow(rr)),colnames(b_both_rr[[y]]),rr,rd)
  res <- data.frame(apply(res,2,as.numeric))
  colnames(res) <- c("MVPA","LPA","rr_est","rr_se","rd_est","rd_se")
  res <- res %>% pivot_longer(c = 3:6,
                              names_to = c("type",".value"),
                              names_pattern = "(.+)_(.+)")
}))

meld_data$MVPA <- factor(meld_data$MVPA,
                         levels=c(1,2,3),
                         labels=c("10th percentile","Median","90th percentile"))
meld_data$type <- factor(meld_data$type)

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
scaleFUN <- function(x) sprintf("%.2f", x)

######################################################################################
# 6. Create figures
#-------------------------------------------------------------------------------------

lpa_rr_figure <- ggplot(res_lpa[which(res_lpa$type=="rr"),],aes(x=LPA,y=est)) + 
  geom_line(colour="#619CFF") + 
  geom_ribbon(aes(ymin=est-qnorm(0.975)*se,ymax=est+qnorm(0.975)*se),alpha=0.3,colour=NA,fill="#619CFF") +
  ylab("Relative risk") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels=scaleFUN, breaks=seq(0.0,1.4,0.2),limits=c(0.0,1.4)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
lpa_rd_figure <- ggplot(res_lpa[which(res_lpa$type=="rd"),],aes(x=LPA,y=est)) + 
  geom_line(colour="#619CFF") + 
  geom_ribbon(aes(ymin=est-qnorm(0.975)*se,ymax=est+qnorm(0.975)*se),alpha=0.3,colour=NA,fill="#619CFF") +
  ylab("Risk difference") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels=scaleFUN, breaks=seq(-0.25,0.1,0.05),limits=c(-0.25,0.1)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
lpa_comb_figure <- ggarrange(lpa_rr_figure,lpa_rd_figure,nrow=2)

ggsave(paste0(workdir,"Results/lpa rr 20250501.png"),lpa_rr_figure,width = 2100,height = 1100,units = "px")
ggsave(paste0(workdir,"Results/lpa rd 20250501.png"),lpa_rd_figure,width = 2100,height = 1100,units = "px")
ggsave(paste0(workdir,"Results/lpa rr and rd 20250501.png"),lpa_comb_figure,width = 2100,height = 2100,units = "px")

mvpa_rr_figure <- ggplot(res_mvpa[which(res_mvpa$type=="rr"),],aes(x=MVPA,y=est)) + 
  geom_line(colour="#00BA38") + 
  geom_ribbon(aes(ymin=est-qnorm(0.975)*se,ymax=est+qnorm(0.975)*se),alpha=0.3,colour=NA,fill="#00BA38") +
  ylab("Relative risk") +
  xlab("Minutes of moderate/vigorous activity per day") + 
  scale_x_continuous(breaks=seq(0,100,10)) +
  scale_y_continuous(labels=scaleFUN, breaks=seq(0.0,1.4,0.2),limits=c(0.0,1.4)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
mvpa_rd_figure <- ggplot(res_mvpa[which(res_mvpa$type=="rd"),],aes(x=MVPA,y=est)) + 
  geom_line(colour="#00BA38") + 
  geom_ribbon(aes(ymin=est-qnorm(0.975)*se,ymax=est+qnorm(0.975)*se),alpha=0.3,colour=NA,fill="#00BA38") +
  ylab("Risk difference") +
  xlab("Minutes of moderate/vigorous activity per day") + 
  scale_x_continuous(breaks=seq(0,100,10)) +
  scale_y_continuous(labels=scaleFUN, breaks=seq(-0.15,0.1,0.05),limits=c(-0.15,0.1)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
mvpa_comb_figure <- ggarrange(mvpa_rr_figure,mvpa_rd_figure,nrow=2)

ggsave(paste0(workdir,"Results/mvpa rr 20250501.png"),mvpa_rr_figure,width = 2100,height = 1100,units = "px")
ggsave(paste0(workdir,"Results/mvpa rd 20250501.png"),mvpa_rd_figure,width = 2100,height = 1100,units = "px")
ggsave(paste0(workdir,"Results/mvpa rr and rd 20250501.png"),mvpa_comb_figure,width = 2100,height = 2100,units = "px")

both_rr_figure <- ggplot(meld_data[which(meld_data$type=="rr"),],aes(x=LPA,y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=est-qnorm(0.975)*se,ymax=est+qnorm(0.975)*se, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Relative risk") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels=scaleFUN, breaks=seq(0.0,1.6,0.2),limits=c(0.0,1.6)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
both_rd_figure <- ggplot(meld_data[which(meld_data$type=="rd"),],aes(x=LPA,y=est, colour=MVPA)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=est-qnorm(0.975)*se,ymax=est+qnorm(0.975)*se, fill=MVPA),alpha=0.3,colour=NA) +
  ylab("Risk difference") +
  xlab("Minutes of light activity per day") + 
  scale_x_continuous(breaks=seq(0,600,60)) +
  scale_y_continuous(labels=scaleFUN, breaks=seq(-0.15,0.1,0.05),limits=c(-0.15,0.1)) +
  theme_light() + 
  theme(legend.position = "bottom") +
  figure_theme
both_comb_figure <- ggarrange(both_rr_figure,both_rd_figure,nrow=2)
both_comb_figure
ggsave(paste0(workdir,"Results/both rr 20250501.png"),both_rr_figure,width = 2100,height = 1100,units = "px")
ggsave(paste0(workdir,"Results/both rd 20250501.png"),both_rd_figure,width = 2100,height = 1100,units = "px")
ggsave(paste0(workdir,"Results/both rr and rd 20250501.png"),both_comb_figure,width = 2100,height = 2100,units = "px")
