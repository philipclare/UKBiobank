######################################################################################
##   
## Analysis of light physical activity and survival
## Multiple imputation
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
libs <- c("mice","miceadds","parallel","VIM")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)

rm(libs)
rm(missing)

# 1.3. Set argument either for testing or passed by scheduler on HPC
if (Sys.info()[['sysname']]=="Linux") {
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
} else if (Sys.info()[['sysname']]=="Windows" | Sys.info()[['sysname']]=="Darwin") {
  args <- 1 # for testing on local computer
}

args <- c(args[1],
          ceiling(args[1]/20),
          ceiling(args[1]-((ceiling(args[1]/20)-1)*20)))

# 1.4 Set seed
set.seed(235377)
seeds <- sample.int(100000, 20)
set.seed(seeds[args[3]])

######################################################################################
# 2. Load and merge data files
#-------------------------------------------------------------------------------------

# 2.1 Load data
data_list <- list(paste0(workdir,"Data/data_for_imputation.rds"),
                  paste0(workdir,"Data/data_for_imputation_sensitivity1.rds"),
                  paste0(workdir,"Data/data_for_imputation_sensitivity2.rds"),
                  paste0(workdir,"Data/data_for_imputation_sensitivity3.rds"))
dataread <- readRDS(file=data_list[[args[2]]])

# Sort from most to least missing
res<-summary(aggr(dataread))$missings
varorder <- res$Variable
res<-res[order(-res$Count),]
dataread <- dataread[,res$Variable]
aggr(dataread)

######################################################################################
# 3. Define Imputation Paramaters
#-------------------------------------------------------------------------------------

m <- 1 # Number of imputations, set to one for Katana so each node runs a single M
maxit <- 100; # Number of mice iterations

######################################################################################
# 4. Imputation
#-------------------------------------------------------------------------------------

# 4.1 Parallel imputation of temp data using parLapply
start_time1 <- Sys.time()
imputation <- mice(data=dataread,
                   m=m,
                   maxit=maxit,
                   defaultMethod = c("rf","rf","rf","rf"))

end_time1 <- Sys.time()

# 4.2 Extract imputation as data frame and return to original variable order
imp_data <- mids2datlist(imputation)[[1]]
imp_data <- imp_data[,varorder]

# 4.3 Calculate and report time taken
time_taken1 <- end_time1 - start_time1

cat('Imputation ', args[2], '-', args[3], 'with ', maxit, 'iterations took:', time_taken1, attr(time_taken1,"units"), ".","\n")

######################################################################################
# 5. Save results
#-------------------------------------------------------------------------------------

# 5.1 Save results as RData
save_list <- list("primary","s1","s2","s3")
saveRDS(imp_data,file=paste0(workdir,"Data/Imputed/Imputed_data_",save_list[[args[2]]],"_",args[3],".rds"))