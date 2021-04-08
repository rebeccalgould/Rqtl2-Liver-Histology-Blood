# Script created by Rebecca Gould, UGA
# Using Rdata environment created by 1-Rqtl2-run.R script


####################################################
## Load in necessary packages
####################################################

library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (pander)
library (writexl)
library (RSQLite)


##################################################################
## Evaluate gluthatione synthesis and recycling genes
##################################################################

## Liver GSH ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_LiverGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gpx1 Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_LiverGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_LiverGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gclc Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_LiverGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_LiverGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gclm Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_LiverGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_LiverGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gss Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_LiverGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_LiverGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Gsr Position -- Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


## Liver GSSG ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_LiverGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gpx1 Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_LiverGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_LiverGSSG_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gclc Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_LiverGSSG_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSSG_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_LiverGSSG_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gclm Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_LiverGSSG_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_LiverGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gss Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_LiverGSSG_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSSG_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_LiverGSSG_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSSG, main = "Gsr Position -- Liver GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  

## Liver Total Glutathione ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_LiverTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Gpx1 Position -- Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_LiverTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_LiverTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Gclc Position -- Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_LiverTotalGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Gclm Position -- Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_LiverTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Gss Position -- Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_LiverTotalGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Gsr Position -- Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  

## Liver GSH/GSSG Ratio ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_LiverGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gpx1 Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_LiverGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gclc Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_LiverGSH_GSSGRatio_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gclm Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_LiverGSH_GSSGRatio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gss Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_LiverGSH_GSSGRatio_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Gsr Position -- Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  
## Liver Glutathione Redox Potential ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_LiverRedoxPotentialGSSG2GSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Gpx1 Position -- Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_LiverRedoxPotentialGSSG2GSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr9, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Gclc Position -- Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_LiverRedoxPotentialGSSG2GSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Gclm Position -- Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_LiverRedoxPotentialGSSG2GSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Gss Position -- Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_LiverRedoxPotentialGSSG2GSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Gsr Position -- Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  
  
  
  