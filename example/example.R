
# Library packages ----------------------------------------------------------------------------

# install.packages("remotes")
# install.packages("devtools")
# remotes::install_github(repo="Bai-Li-NOAA/Age_Structured_Stock_Assessment_Model_Comparison", force=TRUE)
# library(ASSAMC)

setwd("C:/Users/bai.li/Documents/Age_Structured_Stock_Assessment_Model_Comparison_spatial/")
devtools::load_all()
## Need to install packages below:
## ASAPplots, r4ss, readxl, PBSadmb
#devtools::install_github("cmlegault/ASAPplots", build_vignettes = TRUE)
library(readxl)
library(PBSadmb)
library(ASAPplots)
library(r4ss)


# Working directory settings ------------------------------------------------------------------

maindir <- "C:/Users/bai.li/Documents/Age_Structured_Stock_Assessment_Model_Comparison_spatial/example"

# Basic simulation settings -------------------------------------------------------------------

om_sim_num <- 160 # total number of iterations per case
keep_sim_num <- 100 # number of kept iterations per case
figure_number <- 10 # number of individual iteration to plot

seed_num <- 9924

# Basic stock settings ------------------------------------------------------------------------------

year <- 1:30
num_stock <- 3
stocks <- vector(mode="list", length=num_stock)
names(stocks) <- paste("stock", 1:num_stock, sep="")

# Movement Settings ---------------------------------------------------------------------
movement_matrix <- lapply(1:length(year),
                          function(x)
                            matrix(c(0.68, 0.22, 0.10,
                                     0.24, 0.37, 0.39,
                                     0.08, 0.28, 0.64),
                                   ncol=3, byrow=T))


# Life-history settings -----------------------------------------------------------------------

ages <- 1:12   #Age structure of the popn

initial_equilibrium_F <- TRUE
median_R0 <- 1000000 #Average annual unfished recruitment (scales the popn)
median_h <- 0.75 #Steepness of the Beverton-Holt spawner-recruit relationship.
mean_R0 <- NULL
mean_h <- NULL
SRmodel <- 1 # 1=Beverton-Holt; 2=Ricker
M <- 0.2       #Age-invariant natural mortality

Linf <- 800	  #Asymptotic average length
K <- 0.18     	#Growth coefficient
a0 <- -1.36    #Theoretical age at size 0
a.lw <- 0.000000025  #Length-weight coefficient
b.lw <- 3.0     	    #Length-weight exponent
A50.mat <- 2.25      #Age at 50% maturity
slope.mat <- 3       #Slope of maturity ogive
pattern.mat <- 1     #Simple logistic maturity
female.proportion <- 0.5   #Sex ratio

#### Fleet settings ####
fleet_num <- 1

#CV of landings for OM
cv.L <- list()
cv.L$fleet1 <- 0.005

#Input CV of landings for EMs
input.cv.L <- list()
input.cv.L$fleet1 <- 0.01

#Annual sample size (nfish) of age comp samples
n.L <- list()
n.L$fleet1 <- 200

#Define fleet selectivity
sel_fleet <- list()

sel_fleet$fleet1$pattern <- 1
sel_fleet$fleet1$A50.sel1 <- 2
sel_fleet$fleet1$slope.sel1 <- 1

#### Survey settings ####
survey_num <- 1

#CV of surveys for OM
cv.survey <- list()
cv.survey$survey1 <- 0.1

#Input CV of surveys for EMs
input.cv.survey <- list()
input.cv.survey$survey1 <- 0.2

#Annual sample size (nfish) of age comp samples
n.survey <- list()
n.survey$survey1 <- 200

#Define survey selectivity
sel_survey <- list()

sel_survey$survey1$pattern <- 1
sel_survey$survey1$A50.sel1 <- 1.5
sel_survey$survey1$slope.sel1 <- 2

#### Other settings ####
logf_sd <- 0.2
f_dev_change <- FALSE
f_pattern <- 1
start_val <- 0.01
middle_val <- NULL
end_val <- 0.39
f_val <- NULL
start_year <- 1
middle_year <- NULL

logR_sd <- 0.6
r_dev_change <- TRUE

om_bias_cor <- TRUE
bias_cor_method <- "median_unbiased" #Options: "none", "median_unbiased", and "mean_unbiased"
em_bias_cor <- TRUE

stocks$stock1 <- save_stock_input(base_stock=TRUE)
stocks$stock2 <- save_stock_input(base_stock=FALSE,
                                  input_list=stocks$stock1,
                                  logR_sd=0.8)
stocks$stock3 <- save_stock_input(base_stock=FALSE,
                                  input_list=stocks$stock1)

#### Base Case (Case 0) ####
base_case_input <- save_initial_input(base_case = TRUE,
                                      input_list = stocks,
                                      case_name = "C0")
#### Run OM
run_om(input_list=base_case_input, show_iter_num=F)
