#' Function to read key output from OM and EMs
#' @name read_plot_data
#' @description Function to read key output from OM and EMs
#' @param em_names Names of estimation models
#' @param casedir Case working directory
#' @param keep_sim_num Final number of iterations saved for a case
#' @param adhoc_bias_cor Use ad hoc bias correction in recruitment? The default setting is FALSE.
#' @param SRmodel Spawner-recruit model. 1 = Beverton-Holt model; 2 = Ricker model
#'
#' @export
## Aggregate data of biomass, abundance, SSB, recruitment, F (apical F*selectivity), F multiplier, landings, and survey from models to matrix
read_plot_data <- function(em_names=NULL, casedir=NULL, keep_sim_num=NULL, adhoc_bias_cor=FALSE, SRmodel=1){

  ## OM
  subdir = "OM"
  load(file.path(casedir, "output", subdir, paste("OM", 1, ".RData", sep="")))

  om_biomass <- om_abundance <-
    om_ssb <- om_recruit <- om_Ftot <- om_Fmul <-
    om_landing <- om_survey <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

  om_msy <- om_fmsy <- om_ssbmsy <- matrix(NA, nrow=1, ncol=keep_sim_num)

  om_fratio <- om_ssbratio <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)
  om_agecomp <- list()

  om_landing_err <- om_survey_err <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)
  om_geomR0 <- om_arimR0 <-
    om_geomS0 <- om_arimS0 <-
    om_geomDf <- om_arimDf <- matrix(NA, nrow=1, ncol=keep_sim_num)

  for (om_sim in 1:keep_sim_num){

    load(file.path(casedir, "output", subdir, paste("OM", keep_sim_id[om_sim], ".RData", sep="")))

    om_biomass[,om_sim] <- om_output$biomass.mt
    om_abundance[,om_sim] <- om_output$abundance/1000
    om_ssb[,om_sim] <- om_output$SSB
    om_recruit[,om_sim] <- om_output$N.age[,1]/1000
    om_Ftot[,om_sim] <- apply(om_output$FAA, 1, max)
    om_Fmul[,om_sim] <- om_output$f
    om_landing[,om_sim] <- om_output$L.mt$fleet1
    om_survey[,om_sim] <- om_output$survey_index$survey1
    om_msy[, om_sim] <- om_output$msy$msy
    om_fmsy[, om_sim] <- round(om_output$msy$Fmsy, digits = 3)
    om_ssbmsy[, om_sim] <- om_output$msy$SSBmsy
    om_fratio[, om_sim] <- om_Ftot[, om_sim]/om_fmsy[om_sim]
    om_ssbratio[, om_sim] <- om_ssb[, om_sim]/om_ssbmsy[om_sim]
    om_agecomp[[om_sim]] <- apply(om_output$N.age/1000, 1, function(x) x/sum(x))
    #om_landing_err[,om_sim] <- em_input$L.obs$fleet1
    #om_survey_err[,om_sim] <- em_input$survey.obs$survey1
    om_geomR0[,om_sim] <- om_input$median_R0/1000
    om_arimR0[,om_sim] <- om_input$mean_R0/1000
    om_geomS0[,om_sim] <- om_input$median_R0*om_input$Phi.0
    om_arimS0[,om_sim] <- om_input$mean_R0*om_input$Phi.0
    om_geomDf[,om_sim] <- om_ssb[nrow(om_ssb),om_sim]/om_geomS0[,om_sim]
    om_arimDf[,om_sim] <- om_ssb[nrow(om_ssb),om_sim]/om_arimS0[,om_sim]

  }

  om_list <- list(om_biomass, om_abundance,
                  om_ssb, om_recruit, om_Ftot,
                  om_landing, om_survey,
                  om_msy, om_fmsy, om_ssbmsy,
                  om_fratio, om_ssbratio,
                  om_geomR0, om_arimR0,
                  om_geomS0, om_arimS0,
                  om_geomDf, om_arimDf,
                  om_agecomp)

  names(om_list) <- c("biomass", "abundance",
                      "ssb", "recruit", "Ftot",
                      "landing", "survey",
                      "msy", "fmsy", "ssbmsy",
                      "fratio", "ssbratio",
                      "geomR0", "arimR0",
                      "geomS0", "arimS0",
                      "geomDf", "arimDf",
                      "agecomp")

  om_list <<- om_list
  save(om_list, file=file.path(casedir, "output", "om_output.RData"))

  ## AMAK
  if ("AMAK" %in% em_names){

    amak_biomass <- amak_abundance <-
      amak_ssb <- amak_recruit <- amak_Ftot <- amak_Fmul <-
      amak_landing <- amak_survey <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    amak_msy <- amak_fmsy <- amak_ssbmsy <- matrix(NA, nrow=1, ncol=keep_sim_num)

    amak_fratio <- amak_ssbratio <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    amak_agecomp <- list()

    amak_geomR0 <- amak_arimR0 <-
      amak_geomS0 <- amak_arimS0 <-
      amak_geomDf <- amak_arimDf <- matrix(NA, nrow=1, ncol=keep_sim_num)

    subdir = "AMAK"
    for (om_sim in 1:keep_sim_num){

      setwd(file.path(casedir, "output", subdir, paste("s", keep_sim_id[om_sim], sep="")))
      amak_output <- readRep("For_R", suffix = ".rep")
      amak_std <- readRep("amak", suffix = ".std")

      amak_biomass[,om_sim] <- amak_output$TotBiom[which(amak_output$TotBiom[,1]>=om_input$year[1] & amak_output$TotBiom[,1]<=om_input$year[length(om_input$year)]),2]
      amak_abundance[,om_sim] <- apply(amak_output$N[,2:ncol(amak_output$N)], 1, sum)/1000
      amak_ssb[,om_sim] <- amak_output$SSB[which(amak_output$SSB[,1]>=om_input$year[1] &  amak_output$SSB[,1]<= om_input$year[length(om_input$year)]),2]
      amak_recruit[,om_sim] <- amak_output$R[,2]/1000
      amak_Ftot[,om_sim] <- apply(amak_output$TotF, 1, max)
      amak_Fmul[,om_sim] <- NA
      amak_landing[,om_sim] <- amak_output$Pred_catch_1
      amak_survey[,om_sim] <- amak_output$Obs_Survey_1[,3]

      if(adhoc_bias_cor==TRUE){
        amak_msy_adhoc=msy_calcs(steep=ifelse(SRmodel==2, log(4*amak_output$Steep[2]/(1-amak_output$Steep[2])), amak_output$Steep[2]),
                                 R0=exp(amak_std$value[which(amak_std$name=="log_Rzero")]),
                                 M=amak_output$M,
                                 wgt=amak_output$wt_a_pop,
                                 prop.f=om_input$proportion.female,
                                 selL=amak_output$sel_fsh_1[1,3:ncol(amak_output$sel_fsh_1)],
                                 selD=rep(0,om_input$nages),
                                 selZ=amak_output$sel_fsh_1[1,3:ncol(amak_output$sel_fsh_1)],
                                 mat.f=amak_output$mature_a/om_input$proportion.female,
                                 mat.m=NULL,
                                 sigma=amak_output$sigmar,
                                 maxF=om_output$msy$maxF,
                                 step=om_output$msy$step,
                                 om_bias_cor=TRUE,
                                 bias_cor_method="median_unbiased",
                                 SRmodel=SRmodel
                                 )
        amak_msy[, om_sim] <- amak_msy_adhoc$msy
        amak_fmsy[, om_sim] <- round(amak_msy_adhoc$Fmsy, digits = 3)
        amak_ssbmsy[, om_sim] <- amak_msy_adhoc$SSBmsy
      } else {
          amak_msy[, om_sim] <- amak_std$value[which(amak_std$name=="MSY")]
          amak_fmsy[, om_sim] <- round(amak_std$value[which(amak_std$name=="Fmsy")], digits = 3)
          amak_ssbmsy[, om_sim] <- amak_std$value[which(amak_std$name=="Bmsy")]
      }

      amak_fratio[, om_sim] <- amak_Ftot[, om_sim]/amak_fmsy[om_sim]
      amak_ssbratio[, om_sim] <- amak_ssb[,om_sim]/amak_ssbmsy[om_sim]
      amak_agecomp[[om_sim]] <- apply(amak_output$N[,2:ncol(amak_output$N)]/1000, 1, function(x) x/sum(x))

      amak_geomR0[, om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
      amak_geomS0[, om_sim] <- amak_geomR0[, om_sim]*amak_output$phizero*1000
      amak_geomDf[, om_sim] <- amak_ssb[nrow(amak_ssb),om_sim]/amak_geomS0[, om_sim]

      SRparms <- convertSRparms(R0=exp(amak_std$value[which(amak_std$name=="log_Rzero")]),
                                h=ifelse(SRmodel==2, log(4*amak_output$Steep[2]/(1-amak_output$Steep[2])), amak_output$Steep[2]),
                                phi=amak_output$phizero,
                                sigmaR=amak_output$sigmar,
                                mean2med=FALSE,
                                model=SRmodel)
      amak_arimR0[, om_sim] <- SRparms$R0BC/1000
      amak_arimS0[, om_sim] <- SRparms$S0BC
      amak_arimDf[, om_sim] <- amak_ssb[nrow(amak_ssb),om_sim]/amak_arimS0[, om_sim]
    }
      amak_list <- list(amak_biomass, amak_abundance,
                        amak_ssb, amak_recruit, amak_Ftot,
                        amak_landing, amak_survey,
                        amak_msy, amak_fmsy, amak_ssbmsy,
                        amak_fratio, amak_ssbratio,
                        amak_geomR0, amak_arimR0,
                        amak_geomS0, amak_arimS0,
                        amak_geomDf, amak_arimDf,
                        amak_agecomp)
      names(amak_list) <- c("biomass", "abundance",
                            "ssb", "recruit", "Ftot",
                            "landing", "survey",
                            "msy", "fmsy", "ssbmsy",
                            "fratio", "ssbratio",
                            "geomR0", "arimR0",
                            "geomS0", "arimS0",
                            "geomDf", "arimDf",
                            "agecomp")
      amak_list <<- amak_list
      save(amak_list, file=file.path(casedir, "output", "amak_output.RData"))
  }

  ## ASAP
  if ("ASAP" %in% em_names){

    asap_biomass <- asap_abundance <-
      asap_ssb <- asap_recruit <- asap_Ftot <- asap_Fmul <-
      asap_landing <- asap_survey <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    asap_msy <- asap_fmsy <- asap_ssbmsy <- matrix(NA, nrow=1, ncol=keep_sim_num)

    asap_fratio <- asap_ssbratio <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    asap_agecomp <- list()

    asap_geomR0 <- asap_arimR0 <-
      asap_geomS0 <- asap_arimS0 <-
      asap_geomDf <- asap_arimDf <- matrix(NA, nrow=1, ncol=keep_sim_num)

    subdir = "ASAP"

    for (om_sim in 1:keep_sim_num){

      asap_output <- dget(file.path(casedir, "output", subdir, paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
      setwd(file.path(casedir, "output", subdir, paste("s", keep_sim_id[om_sim], sep="")))
      asap_std <- readRep("asap3", suffix = ".std")

      asap_biomass[,om_sim] <- asap_output$tot.jan1.B
      asap_abundance[,om_sim] <- apply(asap_output$N.age, 1, sum)
      asap_ssb[,om_sim] <- asap_output$SSB
      asap_recruit[,om_sim] <- asap_output$N.age[,1]
      asap_Ftot[,om_sim] <- apply(asap_output$fleet.FAA$FAA.directed.fleet1, 1, max)
      asap_Fmul[,om_sim] <- asap_output$fleet.Fmult
      asap_landing[,om_sim] <- asap_output$catch.pred
      asap_survey[,om_sim] <- asap_output$index.pred$ind01

      if(adhoc_bias_cor==TRUE) {
        asap_msy_adhoc <- msy_calcs(steep=asap_output$SR.parms$SR.steepness,
                                    R0=asap_output$SR.parms$SR.R0*1000,
                                    M=asap_output$M.age,
                                    wgt=asap_output$WAA.mats$WAA.jan1[1,]/1000,
                                    prop.f=om_input$proportion.female,
                                    selL=asap_output$fleet.sel.mats$sel.m.fleet1[1,],
                                    selD=rep(0,om_input$nages),
                                    selZ=asap_output$fleet.sel.mats$sel.m.fleet1[1,],
                                    mat.f=asap_output$maturity[1,]/om_input$proportion.female,
                                    mat.m=NULL,
                                    sigma=sqrt(log(asap_output$control.parms$recruit.cv[1]^2+1)),
                                    maxF=om_output$msy$maxF,
                                    step=om_output$msy$step,
                                    om_bias_cor=TRUE,
                                    bias_cor_method="median_unbiased",
                                    SRmodel=SRmodel)
        asap_msy[, om_sim] <- asap_msy_adhoc$msy
        asap_fmsy[, om_sim] <- round(asap_msy_adhoc$Fmsy, digits = 3)
        asap_ssbmsy[, om_sim] <- asap_msy_adhoc$SSBmsy
      } else {
       asap_msy[, om_sim] <- asap_std$value[which(asap_std$name=="MSY")]
       asap_fmsy[, om_sim] <- round(asap_std$value[which(asap_std$name=="Fmsy_report")], digits = 3)
       asap_ssbmsy[, om_sim] <- asap_std$value[which(asap_std$name=="SSBmsy_report")]
      }

      asap_fratio[, om_sim] <- asap_Ftot[, om_sim]/asap_fmsy[om_sim]
      asap_ssbratio[, om_sim] <- asap_ssb[,om_sim]/asap_ssbmsy[om_sim]
      asap_agecomp[[om_sim]] <- apply(asap_output$N.age, 1, function(x) x/sum(x))

      asap_geomR0[, om_sim] <- asap_output$SR.parms$SR.R0
      asap_geomS0[, om_sim] <- asap_output$SR.parms$SR.S0
      asap_geomDf[, om_sim] <- asap_ssb[nrow(asap_ssb),om_sim]/asap_geomS0[, om_sim]

      SRparms <- convertSRparms(R0=asap_output$SR.parms$SR.R0,
                                h=asap_output$SR.parms$SR.steepness,
                                phi=asap_output$SR.parms$SR.SPR0,
                                sigmaR=sqrt(log(asap_output$control.parms$recruit.cv[1]^2+1)),
                                mean2med=FALSE,
                                model=SRmodel)
      asap_arimR0[, om_sim] <- SRparms$R0BC
      asap_arimS0[, om_sim] <- SRparms$S0BC
      asap_arimDf[, om_sim] <- asap_ssb[nrow(asap_ssb),om_sim]/asap_arimS0[, om_sim]

    }

    asap_list <- list(asap_biomass, asap_abundance,
                      asap_ssb, asap_recruit, asap_Ftot,
                      asap_landing, asap_survey,
                      asap_msy, asap_fmsy, asap_ssbmsy,
                      asap_fratio, asap_ssbratio,
                      asap_geomR0, asap_arimR0,
                      asap_geomS0, asap_arimS0,
                      asap_geomDf, asap_arimDf,
                      asap_agecomp)

    names(asap_list) <- c("biomass", "abundance",
                          "ssb", "recruit", "Ftot",
                          "landing", "survey",
                          "msy", "fmsy", "ssbmsy",
                          "fratio", "ssbratio",
                          "geomR0", "arimR0",
                          "geomS0", "arimS0",
                          "geomDf", "arimDf",
                          "agecomp")

    asap_list <<- asap_list
    save(asap_list, file=file.path(casedir, "output", "asap_output.RData"))
  }

  ## BAM
  if("BAM" %in% em_names){
    bam_biomass <- bam_abundance <-
      bam_ssb <- bam_recruit <- bam_Ftot <- bam_Fmul <-
      bam_landing <- bam_survey <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    bam_msy <- bam_fmsy <- bam_ssbmsy <- matrix(NA, nrow=1, ncol=keep_sim_num)

    bam_fratio <- bam_ssbratio <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    bam_agecomp <- list()

    bam_geomR0 <- bam_arimR0 <-
      bam_geomS0 <- bam_arimS0 <-
      bam_geomDf <- bam_arimDf <- matrix(NA, nrow=1, ncol=keep_sim_num)

    subdir = "BAM"

    for (om_sim in 1:keep_sim_num){

      bam_output <- dget(file.path(casedir, "output", subdir, paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))

      bam_biomass[,om_sim] <- bam_output$t.series$B[1:om_input$year[length(om_input$year)]]
      bam_abundance[,om_sim] <- bam_output$t.series$N[1:om_input$year[length(om_input$year)]]/1000
      bam_ssb[,om_sim] <- bam_output$t.series$SSB[1:om_input$year[length(om_input$year)]]
      bam_recruit[,om_sim] <- bam_output$t.series$recruits[1:om_input$year[length(om_input$year)]]/1000
      bam_Ftot[,om_sim] <- bam_output$t.series$F.full[1:om_input$year[length(om_input$year)]]
      bam_Fmul[,om_sim] <- NA
      bam_landing[,om_sim] <- bam_output$t.series$L.fleet1.pr[1:om_input$year[length(om_input$year)]]
      bam_survey[,om_sim] <- bam_output$t.series$U.survey1.pr[1:om_input$year[length(om_input$year)]]
      bam_msy[, om_sim] <- bam_output$parms$msy.mt
      bam_fmsy[, om_sim] <- round(bam_output$parms$Fmsy, digits = 3)
      bam_ssbmsy[, om_sim] <- bam_output$parms$SSBmsy
      bam_fratio[, om_sim] <- bam_output$t.series$F.Fmsy[1:om_input$year[length(om_input$year)]]
      bam_ssbratio[, om_sim] <- bam_output$t.series$SSB.SSBmsy[1:om_input$year[length(om_input$year)]]
      bam_agecomp[[om_sim]] <- apply(bam_output$N.age[1:om_input$nyr,]/1000, 1, function(x) x/sum(x))
      bam_geomR0[, om_sim] <- bam_output$parms$R0/1000
      bam_arimR0[, om_sim] <- bam_output$parms$R.virgin.bc/1000
      bam_geomS0[, om_sim] <- bam_output$parms$SSB0
      bam_arimS0[, om_sim] <- bam_output$parms$R.virgin.bc*bam_output$parms[[grep("Phi0", names(bam_output$parms), value=TRUE)]]
      bam_geomDf[, om_sim] <- bam_ssb[nrow(bam_ssb),om_sim]/bam_geomS0[, om_sim]
      bam_arimDf[, om_sim] <- bam_ssb[nrow(bam_ssb),om_sim]/bam_arimS0[, om_sim]

    }

    bam_list <- list(bam_biomass, bam_abundance,
                     bam_ssb, bam_recruit, bam_Ftot,
                     bam_landing, bam_survey,
                     bam_msy, bam_fmsy, bam_ssbmsy,
                     bam_fratio, bam_ssbratio,
                     bam_geomR0, bam_arimR0,
                     bam_geomS0, bam_arimS0,
                     bam_geomDf, bam_arimDf,
                     bam_agecomp)

    names(bam_list) <- c("biomass", "abundance",
                         "ssb", "recruit", "Ftot",
                         "landing", "survey",
                         "msy", "fmsy", "ssbmsy",
                         "fratio", "ssbratio",
                         "geomR0", "arimR0",
                         "geomS0", "arimS0",
                         "geomDf", "arimDf",
                         "agecomp")
    bam_list <<- bam_list
    save(bam_list, file=file.path(casedir, "output", "bam_output.RData"))
  }

  ## SS
  if ("SS" %in% em_names){
    library(r4ss)
    ss_biomass <- ss_abundance <-
      ss_ssb <- ss_recruit <- ss_Ftot <- ss_Fmul <-
      ss_landing <- ss_survey <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    ss_msy <- ss_fmsy <- ss_ssbmsy <- matrix(NA, nrow=1, ncol=keep_sim_num)

    ss_fratio <- ss_ssbratio <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)

    ss_agecomp <- list()

    ss_geomR0 <- ss_arimR0 <-
      ss_geomS0 <- ss_arimS0 <-
      ss_geomDf <- ss_arimDf <- matrix(NA, nrow=1, ncol=keep_sim_num)

    subdir = "SS"
    for (om_sim in 1:keep_sim_num){
      ss_output <- SS_output(
        dir=file.path(casedir, "output", subdir, paste("s", keep_sim_id[om_sim], sep="")),
        ncols = 300,
        verbose=F, printstats=F)

      setwd(file.path(casedir, "output", subdir, paste("s", keep_sim_id[om_sim], sep="")))
      ss_std <- readRep("ss", suffix = ".std")
      msy_std <- ss_std[ss_std$name=="Mgmt_quant",]

      ss_biomass[,om_sim] <- ss_output$timeseries$`SmryBio_SX:1_GP:1`[which(ss_output$timeseries$Yr>=om_input$year[1] & ss_output$timeseries$Yr<=om_input$year[length(om_input$year)])]
      ss_abundance[,om_sim] <- ss_output$timeseries$`SmryNum_SX:1_GP:1`[which(ss_output$timeseries$Yr>=om_input$year[1] & ss_output$timeseries$Yr<=om_input$year[length(om_input$year)])]
      ss_ssb[,om_sim] <- ss_output$timeseries$SpawnBio[which(ss_output$timeseries$Yr>=om_input$year[1] & ss_output$timeseries$Yr<=om_input$year[length(om_input$year)])]
      ss_recruit[,om_sim] <- ss_output$natage_annual_2_with_fishery[which(ss_output$natage_annual_2_with_fishery$Yr>=om_input$year[1] & ss_output$natage_annual_2_with_fishery$Yr<=om_input$year[length(om_input$year)]),5]
      ss_Ftot[,om_sim] <- ss_output$timeseries$`F:_1`[which(ss_output$timeseries$Yr>=om_input$year[1] & ss_output$timeseries$Yr<=om_input$year[length(om_input$year)])]
      ss_Fmul[,om_sim] <- NA
      ss_landing[,om_sim] <- ss_output$timeseries$`sel(B):_1`[which(ss_output$timeseries$Yr>=om_input$year[1] & ss_output$timeseries$Yr<=om_input$year[length(om_input$year)])]
      ss_survey[,om_sim] <- ss_output$cpue$Exp[which(ss_output$cpue$Fleet==2)]
      ss_msy[, om_sim] <- msy_std[15, "value"]
      ss_fmsy[, om_sim] <- round(msy_std[14, "value"], digits = 3)
      ss_ssbmsy[, om_sim] <- msy_std[12, "value"]
      ss_fratio[, om_sim] <- ss_output$Kobe$F.Fmsy[which(ss_output$Kobe$Yr>=om_input$year[1] & ss_output$Kobe$Yr<=om_input$year[length(om_input$year)])]
      ss_ssbratio[, om_sim] <- ss_output$Kobe$B.Bmsy[which(ss_output$Kobe$Yr>=om_input$year[1] & ss_output$Kobe$Yr<=om_input$year[length(om_input$year)])]
      ss_agecomp[[om_sim]] <- apply(ss_output$natage_annual_2_with_fishery[1:om_input$nyr, 5:ncol(ss_output$natage_annual_2_with_fishery)], 1, function(x) x/sum(x))

      ss_arimR0[, om_sim] <- ss_output$timeseries$Recruit_0[which(ss_output$timeseries$Era=="INIT")]
      ss_arimS0[, om_sim] <- ss_output$timeseries$SpawnBio[which(ss_output$timeseries$Era=="INIT")]
      ss_arimDf[, om_sim] <- ss_ssb[nrow(ss_ssb),om_sim]/ss_arimS0[, om_sim]

      if(SRmodel==1){
        SRparms <- convertSRparms(R0=ss_output$timeseries$Recruit_0[which(ss_output$timeseries$Era=="INIT")],
                                  h=ss_output$parameters$Value[ss_output$parameters$Label=="SR_BH_steep"],
                                  phi=ss_arimS0[, om_sim]/ss_arimR0[, om_sim],
                                  sigmaR=ss_output$sigma_R_in,
                                  mean2med=TRUE,
                                  model=SRmodel)
      }

      if(SRmodel==2){
        SRparms <- convertSRparms(R0=ss_output$timeseries$Recruit_0[which(ss_output$timeseries$Era=="INIT")],
                                  h=ss_output$parameters$Value[ss_output$parameters$Label=="SR_Ricker_beta"],
                                  phi=ss_arimS0[, om_sim]/ss_arimR0[, om_sim],
                                  sigmaR=ss_output$sigma_R_in,
                                  mean2med=TRUE,
                                  model=SRmodel)
      }

      ss_geomR0[, om_sim] <- SRparms$R0BC
      ss_geomS0[, om_sim] <- SRparms$S0BC
      ss_geomDf[, om_sim] <- ss_ssb[nrow(ss_ssb),om_sim]/ss_geomS0[, om_sim]

    }
    ss_list <- list(ss_biomass, ss_abundance,
                    ss_ssb, ss_recruit, ss_Ftot,
                    ss_landing, ss_survey,
                    ss_msy, ss_fmsy, ss_ssbmsy,
                    ss_fratio, ss_ssbratio,
                    ss_geomR0, ss_arimR0,
                    ss_geomS0, ss_arimS0,
                    ss_geomDf, ss_arimDf,
                    ss_agecomp)

    names(ss_list) <- c("biomass", "abundance",
                        "ssb", "recruit", "Ftot",
                        "landing", "survey",
                        "msy", "fmsy", "ssbmsy",
                        "fratio", "ssbratio",
                        "geomR0", "arimR0",
                        "geomS0", "arimS0",
                        "geomDf", "arimDf",
                        "agecomp")
    ss_list <<- ss_list
    save(ss_list, file=file.path(casedir, "output", "ss_output.RData"))
  }

  ## MAS
  if ("MAS" %in% em_names){

    mas_biomass <- mas_abundance <- mas_ssb <- mas_recruit <- mas_Ftot <- mas_Fmul <- mas_landing <- mas_survey <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)
    mas_msy <- mas_fmsy <- mas_ssbmsy <- matrix(NA, nrow=1, ncol=keep_sim_num)
    mas_fratio <- mas_ssbratio <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)
    mas_agecomp <- list()

    mas_geomR0 <- mas_arimR0 <-
      mas_geomS0 <- mas_arimS0 <-
      mas_geomDf <- mas_arimDf <- matrix(NA, nrow=1, ncol=keep_sim_num)

    subdir = "MAS"
    for (om_sim in 1:keep_sim_num){
      mas_output <- read_json(file.path(casedir, "output",  subdir, paste("s", keep_sim_id[om_sim], sep=""), paste("s", keep_sim_id[om_sim], ".json", sep="")))
      popdy<-mas_output$population_dynamics
      pop<-popdy$populations[[1]]
      flt<-popdy$fleets[[1]]
      srvy<-popdy$surveys[[1]]

      mas_biomass[,om_sim] <- unlist(pop$undifferentiated$biomass$values)
      mas_abundance[,om_sim] <- unlist(pop$undifferentiated$abundance$values)
      mas_ssb[,om_sim] <- unlist(pop$undifferentiated$spawning_stock_biomass$values)
      mas_recruit[,om_sim] <- unlist(pop$undifferentiated$recruits$values)
      mas_Ftot[,om_sim] <- unlist(pop$undifferentiated$fishing_mortality$values)
      mas_Fmul[,om_sim] <- NA
      mas_landing[,om_sim] <- unlist(flt$undifferentiated$catch_biomass$values)
      mas_survey[,om_sim] <- unlist(srvy$undifferentiated$survey_biomass$values)
      mas_msy[, om_sim] <- pop$MSY$msy
      mas_fmsy[, om_sim] <- round(pop$MSY$F_msy, digits = 3)
      mas_ssbmsy[, om_sim] <- pop$MSY$SSB_msy
      mas_fratio[, om_sim] <- mas_Ftot[, om_sim]/mas_fmsy[om_sim]
      mas_ssbratio[, om_sim] <- mas_ssb[,om_sim]/mas_ssbmsy[om_sim]
      mas_agecomp[[om_sim]] <- apply(matrix(unlist(pop$undifferentiated$numbers_at_age$values), nrow=popdy$nyears, ncol=popdy$nages, byrow = T), 1, function(x) x/sum(x))

      parameter <- unlist(mas_output$estimated_parameters$parameters)
      parameter_table <- as.data.frame(matrix(parameter, ncol = 3, byrow = TRUE))
      colnames(parameter_table) <- c(
        "Parameter",
        "Value",
        "Gradient"
      )

      mas_geomR0[, om_sim] <- exp(as.numeric(parameter_table$Value[parameter_table$Parameter == "log_R0_1"]))
      mas_geomS0[, om_sim] <- 0
      mas_geomDf[, om_sim] <- 0
      mas_arimR0[, om_sim] <- exp(as.numeric(parameter_table$Value[parameter_table$Parameter == "log_R0_1"]))
      mas_arimS0[, om_sim] <- 0
      mas_arimDf[, om_sim] <- 0

    }
    mas_list <- list(mas_biomass, mas_abundance,
                     mas_ssb, mas_recruit, mas_Ftot,
                     mas_landing, mas_survey,
                     mas_msy, mas_fmsy, mas_ssbmsy,
                     mas_fratio, mas_ssbratio,
                     mas_geomR0, mas_arimR0,
                     mas_geomS0, mas_arimS0,
                     mas_geomDf, mas_arimDf,
                     mas_agecomp)
    names(mas_list) <- c("biomass", "abundance",
                         "ssb", "recruit", "Ftot",
                         "landing", "survey",
                         "msy", "fmsy", "ssbmsy",
                         "fratio", "ssbratio",
                         "geomR0", "arimR0",
                         "geomS0", "arimS0",
                         "geomDf", "arimDf",
                         "agecomp")
    mas_list <<- mas_list
    save(mas_list, file=file.path(casedir, "output", "mas_output.RData"))
  }

  ## FIMS
  if ("FIMS" %in% em_names){

    fims_biomass <- fims_abundance <- fims_ssb <- fims_recruit <- fims_Ftot <- fims_Fmul <- fims_landing <- fims_survey <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)
    fims_msy <- fims_fmsy <- fims_ssbmsy <- matrix(NA, nrow=1, ncol=keep_sim_num)
    fims_fratio <- fims_ssbratio <- matrix(NA, nrow=om_input$nyr, ncol=keep_sim_num)
    fims_agecomp <- list()

    fims_geomR0 <- fims_arimR0 <-
      fims_geomS0 <- fims_arimS0 <-
      fims_geomDf <- fims_arimDf <- matrix(NA, nrow=1, ncol=keep_sim_num)

    subdir = "FIMS"
    for (om_sim in 1:keep_sim_num){

      load(file.path(casedir, "output",  subdir, paste("s", keep_sim_id[om_sim], sep=""), paste("s", keep_sim_id[om_sim], ".RData", sep="")))

      fims_biomass[,om_sim] <- report$biomass[1:om_input$nyr]
      fims_naa <- matrix(report$naa[1:(om_input$nyr*om_input$nages)], nrow = om_input$nyr, byrow = TRUE)
      fims_abundance[,om_sim] <-  apply(fims_naa, 1, sum)/1000
      fims_ssb[,om_sim] <- report$ssb[1:om_input$nyr]
      fims_recruit[,om_sim] <- fims_naa[,1]/1000
      fims_Ftot[,om_sim] <- report$F_mort
      fims_Fmul[,om_sim] <- report$F_mort
      fims_landing[,om_sim] <- report$expected_index[,1]
      fims_survey[,om_sim] <- em_input$survey.obs$survey1
      fims_msy[, om_sim] <- om_output$msy$msy*1.01
      fims_fmsy[, om_sim] <- round(om_output$msy$Fmsy, digits = 3)*1.01
      fims_ssbmsy[, om_sim] <- om_output$msy$SSBmsy*1.01
      fims_fratio[, om_sim] <- om_Ftot[, om_sim]/om_fmsy[om_sim]
      fims_ssbratio[, om_sim] <- om_ssb[, om_sim]/om_ssbmsy[om_sim]
      fims_agecomp[[om_sim]] <- apply(fims_naa/1000, 1, function(x) x/sum(x))
      fims_geomR0[, om_sim] <- om_input$median_R0/1000
      fims_geomS0[, om_sim] <- om_input$median_R0*om_input$Phi.0
      fims_geomDf[, om_sim] <- om_ssb[nrow(om_ssb),om_sim]/om_geomS0[,om_sim]
      fims_arimR0[, om_sim] <- om_input$mean_R0/1000
      fims_arimS0[, om_sim] <- om_input$mean_R0*om_input$Phi.0
      fims_arimDf[, om_sim] <- om_ssb[nrow(om_ssb),om_sim]/om_arimS0[,om_sim]

    }
    fims_list <- list(fims_biomass, fims_abundance,
                      fims_ssb, fims_recruit, fims_Ftot,
                      fims_landing, fims_survey,
                      fims_msy, fims_fmsy, fims_ssbmsy,
                      fims_fratio, fims_ssbratio,
                      fims_geomR0, fims_arimR0,
                      fims_geomS0, fims_arimS0,
                      fims_geomDf, fims_arimDf,
                      fims_agecomp)
    names(fims_list) <- c("biomass", "abundance",
                          "ssb", "recruit", "Ftot",
                          "landing", "survey",
                          "msy", "fmsy", "ssbmsy",
                          "fratio", "ssbratio",
                          "geomR0", "arimR0",
                          "geomS0", "arimS0",
                          "geomDf", "arimDf",
                          "agecomp")
    fims_list <<- fims_list
    save(fims_list, file=file.path(casedir, "output", "fims_output.RData"))
  }
}

