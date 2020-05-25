#############################################################################
#
#	Central Europe Tree Growth Project
#
#       Annealing with Full Datasets
#
#       Model 4ABC
#
#       CDC,  April 28, 2020
# 
#############################################################################
library(likelihood)

computer <- "azure"

if(computer == "Charlie")
{
  dat_dir <- "\\\\canhamC6\\c\\Users\\canhamc\\Documents\\Current Manuscripts\\Arne - Carpathians\\Earth Day Files"
  out_dir <- "\\\\canhamC6\\c\\Users\\canhamc\\Documents\\Current Manuscripts\\Arne - Carpathians\\Earth Day Files\\Output"
  code_dir <- "\\\\canhamC6\\c\\Users\\canhamc\\Documents\\Current Manuscripts\\Arne - Carpathians\\Earth Day Files\\Code"
} else {
  dat_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/ABIALB/"
  out_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/ABIALB/"
  code_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/ABIALB/"
}

### Get model functions
setwd(code_dir)
source("Model Functions - EDF versions.R")

# focal species
spp_list <- c("PICABI","FAGSYL","ABIALB","ACEPSE")

index <- 3

### load working target data file
setwd(dat_dir)

load(paste(spp_list[index],"Full Working Dataset with All Precip Variables Nitrogen and Mean Climate.Rdata"))

### specify climate variables

working$temp_k <- working$tave_ann_hydro_k
working$temp_lag1_k <- working$tave_ann_hydro_k_lag1

### best precip variables

if (index %in% c(1,4))     # Picea and Acer use annual precip
{  working$precip <- working$ppt_ann_hydro_mm
   working$precip_lag1 <- working$ppt_ann_lag1_hydro_mm
}

if (index == 2)      # Fagus uses water deficit
{  working$precip <- working$hydro_yr_water_deficit_mm
   working$precip_lag1 <- working$hydro_yr_water_deficit_lag1_mm
}
if (index == 3)      # Abies uses seasonal precip
{  working$precip <- working$seas_prec_hydro_yr_mm
   working$precip_lag1 <- working$seas_prec_hydro_yr_lag1_mm
}

### Set annealing parameters

iterations <- 15000

setwd(out_dir)

load(paste(spp_list[index],"Model 4ABCN Full Results.Rdata")) 


var <- list(mean = "predicted", x = "incr_mm", siteplot = "stdcode", region = "k_new", log=T)

nr <- length(levels(as.factor(working$k_new)))

# set parameter limits
#if (index == 4) { par <- model_4ABC_results$best_pars } else
#                { par <- model_4ABC_full_results$best_pars  }
par <- model_4ABCN_full_results$best_pars

if (index == 1) { nspp <- 7 } else
                { nspp <- 5 } 

par$lambda <- par$lambda[1:nspp]

par_lo <- list(PG = rep(0,num_site_plot),
		sizeX0 = 0, sizeXb = 0,
		ageX0 = 0,
		int = 0, sigma = 0,
        	alpha = 0, beta = 0, C = 0, gamma = -2, D = 1, lambda = rep(0,nspp),

		temp.a = rep(0,nr), temp.lag.a = rep(0,nr),
		temp.b = rep(0.25,nr), temp.lag.b = rep(0.25,nr),
                temp.c = rep(min(working$temp_k),nr),
                temp.lag.c = rep(min(working$temp_k),nr),

                prec.a = rep(0,nr), prec.lag.a = rep(0,nr),
                prec.b = rep(0.5,nr), prec.lag.b = rep(0.5,nr),
                prec.c = rep(min(working$precip),nr), prec.lag.c = rep(min(working$precip),nr)    )

par_hi <- list(PG = rep(250,num_site_plot),
		sizeX0 = 100, sizeXb = 1,
		ageX0 = 1,
                int = 100, sigma = 2,
        	alpha = 4, beta = 3, C = 1000, gamma = 2, D = 3, lambda = rep(1,nspp),

		temp.a = rep(1,nr), temp.lag.a = rep(1,nr),
		temp.b = rep(5000,nr), temp.lag.b = rep(5000,nr),
                temp.c = rep(max(working$temp_k),nr),temp.lag.c = rep(max(working$temp_k),nr),

                prec.a = rep(1,nr), prec.lag.a = rep(1,nr),
                prec.b = rep(10000,nr), prec.lag.b = rep(10000,nr),
                prec.c = rep(max(working$precip),nr), prec.lag.c = rep(max(working$precip),nr)    )

#par$n.c <- 15
#par$n.b <- 20
par_lo$n.c <- 0
par_lo$n.b <- 0
par_hi$n.c <- 40
par_hi$n.b <- 200


model_4ABCN_full_azure_results <- anneal(model=model_4ABCN, par=par, var=var, source_data=working, par_lo=par_lo, par_hi=par_hi,
			  pdf=linear_dnorm_with_intercept, dep_var="incr_mm", max_iter=iterations)

setwd(out_dir)

write_results(model_4ABCN_full_azure_results, file = paste(spp_list[index],"Model 4ABCN Full Azure Results.txt", sep=" "), data=F, print_whole_hist=F)
save(model_4ABCN_full_azure_results, file=paste(spp_list[index],"Model 4ABCN Full Azure Results.Rdata",sep=" "))


