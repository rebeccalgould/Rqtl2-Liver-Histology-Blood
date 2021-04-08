# Script created by Rebecca Gould, UGA


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



####################################################
## Read in the control file (gm.json)
####################################################

  R01_GSH_DO_QTLdata <- read_cross2(file = "~/Liver-Glutathione/data/control.json")

  

####################################################
## Genotype probabilities and allele probabilities - calculated by the Jackson Laboratory
####################################################
  
  probs <- readRDS("~/Rqtl2-Glutathione-Genetics/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")


  
####################################################
## Variant files
####################################################

  # set the download timeout to 900 seconds from the default 60 seconds to help with large file downloads
  options(timeout=900) 
  # download the data files needed for SNP association mapping and obtaining genes in QTL intervals
  download.file(url="https://ndownloader.figshare.com/files/18533342", destfile="./data/cc_variants.sqlite") 
  download.file(url="https://ndownloader.figshare.com/files/24607961", destfile="./data/mouse_genes.sqlite")
  download.file(url="https://ndownloader.figshare.com/files/24607970", destfile="./data/mouse_genes_mgi.sqlite")
  
  #Will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice. Make sure they are the most updated versions!
  query_variants <- create_variant_query_func("./data/cc_variants.sqlite")
  query_genes_mgi <- create_gene_query_func("./data/mouse_genes_mgi.sqlite")
  query_genes <- create_gene_query_func("./data/mouse_genes.sqlite")

  

####################################################
## Calculating kinship
####################################################

  kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE, cores = 10)

  

####################################################
## Read in pheno file and adjust for R/qtl2
####################################################

  #to have the phenotype file for reference - can be used when plotting the data to see if it needs to be transformed
  pheno <- read.csv(file = "./data/pheno_covar.csv", header = TRUE)
  
  #make row names the ID of each sample
  rownames(pheno) <- pheno$id
  
  #change sex to numeric variables
  pheno$sex[pheno$sex == "M"] <- 1
  pheno$sex[pheno$sex == "F"] <- 0
  
  #both added covariates must be numeric, not characters 
  pheno$sex <- as.numeric(pheno$sex)
  pheno$generation <- as.numeric(pheno$generation)  


  
####################################################
## Transform data
####################################################
  
  # RankZ function
  rankZ <- function(x) {x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))}
  
  # Rank Z transformations
  pheno$zLiverGSH = rankZ(pheno$Liver_GSH)
  pheno$zLiverGSSG = rankZ(pheno$Liver_GSSG)
  pheno$zLiverTotalGSH = rankZ(pheno$Liver_Total_GSH)
  pheno$zLiverGSH_GSSGRatio = rankZ(pheno$Liver_GSH_GSSG_Ratio)
  pheno$zLiverRedoxPotentialGSSG2GSH = rankZ(pheno$Liver_Redox_Potential_GSSG_2GSH)
  pheno$zLiverNADH = rankZ(pheno$Liver_NADH)
  pheno$zLiverNADP = rankZ(pheno$Liver_NADP)
  pheno$zLiverNADPH = rankZ(pheno$Liver_NADPH)
  pheno$zLiverNADP_NADPHRatio = rankZ(pheno$Liver_NADP_NADPH_Ratio)


  
####################################################
## Add covariates
####################################################

  #adding sex and generation as covariates
  sexgen = model.matrix(~ sex + generation, data = pheno)[,-1]

  

####################################################
## Liver GSH
## Plot Genome Scans with Permutation Tests
####################################################

  qtlscan_LiverGSH <- scan1(genoprobs = probs, pheno = pheno["zLiverGSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverGSH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH"], addcovar = sexgen, n_perm = 1000, cores=10)

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH = summary(perm_LiverGSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH", ylim = c(0,11))
  abline(h = threshold_LiverGSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksGSH <- find_peaks(scan1_output = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksGSH <- find_peaks(scan1_output = qtlscan_LiverGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

#Liver GSH --- Chromosome 14 
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = 14
  coef_blup_LiverGSH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,20)
  plot_coefCC(x = coef_blup_LiverGSH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH, main = "Liver GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

  #using pmap (Mbp)
  chr = 14
  variants_LiverGSH_chr14 <- query_variants(chr, 21, 25)
  out_snps_LiverGSH_chr14 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                         chr = chr, start = 21, end = 25, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSH_chr14$lod, out_snps_LiverGSH_chr14$snpinfo, main = "Liver GSH SNPs")
  LiverGSH_Genes_MGI_chr14 <- query_genes_mgi(chr = chr, start = 21, end = 25)
  plot(out_snps_LiverGSH_chr14$lod, out_snps_LiverGSH_chr14$snpinfo, drop_hilit=1.5, genes = LiverGSH_Genes_MGI_chr14, main = "Liver GSH Genes MGI")
       
       
       
####################################################
## Liver GSSG
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverGSSG <- scan1(genoprobs = probs, pheno = pheno["zLiverGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverGSSG <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSSG = summary(perm_LiverGSSG, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSSG", ylim = c(0,11))
  abline(h = threshold_LiverGSSG, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksGSSG <- find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksGSSG <- find_peaks(scan1_output = qtlscan_LiverGSSG, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  
  
####################################################
## Liver Total Glutathione
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverTotalGSH<- scan1(genoprobs = probs, pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverTotalGSH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverTotalGSH"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverTotalGSH = summary(perm_LiverTotalGSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Total GSH", ylim = c(0,11))
  abline(h = threshold_LiverTotalGSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
 
  #using gmap (cM)
  gmap_peaksTotalGSH <- find_peaks(scan1_output = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksTotalGSH <- find_peaks(scan1_output = qtlscan_LiverTotalGSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverTotalGSH, alpha = 0.2), peakdrop = 1.0, prob = 0.95, expand2markers = FALSE)
  
  
# Liver Total Glutathione --- Chromosome 14 
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 14
  coef_blup_LiverTotalGSH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,20)
  plot_coefCC(x = coef_blup_LiverTotalGSH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverTotalGSH, main = "Liver Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 14
  variants_LiverTotalGSH_chr14 <- query_variants(chr, 21, 25)
  out_snps_LiverTotalGSH_chr14 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                            chr = chr, start = 21, end = 25, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverTotalGSH_chr14$lod, out_snps_LiverTotalGSH_chr14$snpinfo, main = "Liver Total GSH SNPs")
  LiverTotalGSH_Genes_MGI_chr14 <- query_genes_mgi(chr = chr, start = 21, end = 25)
  plot(out_snps_LiverTotalGSH_chr14$lod, out_snps_LiverTotalGSH_chr14$snpinfo, drop_hilit=1.5, genes = LiverTotalGSH_Genes_MGI_chr14, main = "Liver Total GSH Genes MGI")
  
  
  
####################################################
## Liver GSH/GSSG Ratio
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverGSH_GSSGRatio<- scan1(genoprobs = probs, pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverGSH_GSSGRatio <- scan1perm(genoprobs = probs, pheno = pheno["zLiverGSH_GSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH_GSSGRatio = summary(perm_LiverGSH_GSSGRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH/GSSG Ratio", ylim = c(0,11))
  abline(h = threshold_LiverGSH_GSSGRatio, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_LiverGSH_GSSGRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

# Liver GSH/GSSG --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 16
  coef_blup_LiverGSH_GSSGRatio_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,9)
  plot_coefCC(x = coef_blup_LiverGSH_GSSGRatio_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_GSSGRatio, main = "Liver GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 16
  variants_LiverGSH_GSSGRatio_chr16 <- query_variants(chr, 7, 11)
  out_snps_LiverGSH_GSSGRatio_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                 chr = chr, start = 7, end = 11, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSH_GSSGRatio_chr16$lod, out_snps_LiverGSH_GSSGRatio_chr16$snpinfo, main = "Liver GSH/GSSG SNPs")
  LiverGSH_GSSGRatio_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = 7, end = 11)
  plot(out_snps_LiverGSH_GSSGRatio_chr16$lod, out_snps_LiverGSH_GSSGRatio_chr16$snpinfo, drop_hilit=1.5, genes = LiverGSH_GSSGRatio_Genes_MGI_chr16, main = "Liver GSH/GSSG Genes MGI")

  
  
####################################################
## Liver Glutathione Redox Potential
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverRedoxPotentialGSSG2GSH<- scan1(genoprobs = probs, pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverRedoxPotentialGSSG2GSH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], addcovar = sexgen, n_perm = 1000, cores=10)

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverRedoxPotentialGSSG2GSH = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Redox Potential GSSG/2GSH", ylim = c(0,11))
  abline(h = threshold_LiverRedoxPotentialGSSG2GSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksLiverRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksLiverRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

# Liver Glutathione Redox Potential --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 16
  coef_blup_LiverRedoxPotentialGSSG2GSH_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,12)
  plot_coefCC(x = coef_blup_LiverRedoxPotentialGSSG2GSH_chr16, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverRedoxPotentialGSSG2GSH, main = "Liver Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 16
  variants_LiverRedoxPotentialGSSG2GSH_chr16 <- query_variants(chr, 7, 11)
  out_snps_LiverRedoxPotentialGSSG2GSH_chr16 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                          chr = chr, start = 7, end = 11, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverRedoxPotentialGSSG2GSH_chr16$lod, out_snps_LiverRedoxPotentialGSSG2GSH_chr16$snpinfo, main = "Liver Redox Potential GSSG/2GSH SNPs")
  LiverRedoxPotentialGSSG2GSH_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = 7, end = 11)
  plot(out_snps_LiverRedoxPotentialGSSG2GSH_chr16$lod, out_snps_LiverRedoxPotentialGSSG2GSH_chr16$snpinfo, drop_hilit=1.5, genes = LiverRedoxPotentialGSSG2GSH_Genes_MGI_chr16, main = "Liver Redox Potential GSSG/2GSH Genes MGI")
  


####################################################
## Liver NADH
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverNADH <- scan1(genoprobs = probs, pheno = pheno["zLiverNADH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverNADH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverNADH"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADH = summary(perm_LiverNADH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADH", ylim = c(0,11))
  abline(h = threshold_LiverNADH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksNADH <- find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksNADH <- find_peaks(scan1_output = qtlscan_LiverNADH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)


  
####################################################
## Liver NADP
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverNADP <- scan1(genoprobs = probs, pheno = pheno["zLiverNADP"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverNADP <- scan1perm(genoprobs = probs, pheno = pheno["zLiverNADP"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP = summary(perm_LiverNADP, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP", ylim = c(0,11))
  abline(h = threshold_LiverNADP, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksNADP <- find_peaks(scan1_output = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksNADP <-  find_peaks(scan1_output = qtlscan_LiverNADP, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADP, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

# Liver NADP --- Chromosome 3
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 3
  coef_blup_LiverNADP_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42,55)
  plot_coefCC(x = coef_blup_LiverNADP_chr3, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 3
  start = peaksNADP[peaksNADP$chr ==  chr,"ci_lo"]
  end = peaksNADP[peaksNADP$chr == chr, "ci_hi"] 

  variants_LiverNADP_chr3 <- query_variants(chr, start - 1, end + 1)
  out_snps_LiverNADP_chr3 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = start-1, end = end+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_chr3$lod, out_snps_LiverNADP_chr3$snpinfo, main = "Liver NADP SNPs")
  LiverNADP_Genes_MGI_chr3 <- query_genes_mgi(chr = chr, start = start-1, end = end+1)
  plot(out_snps_LiverNADP_chr3$lod, out_snps_LiverNADP_chr3$snpinfo, drop_hilit=1.5, genes = LiverNADP_Genes_MGI_chr3, main = "Liver NADP Genes MGI")
  
# Liver NADP --- Chromosome 8
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 8
  coef_blup_LiverNADP_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(25,35)
  plot_coefCC(x = coef_blup_LiverNADP_chr8, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP, main = "Liver NADP BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 8
  variants_LiverNADP_chr8 <- query_variants(chr, 60, 66)
  out_snps_LiverNADP_chr8 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = 60, end = 66, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_chr8$lod, out_snps_LiverNADP_chr8$snpinfo, main = "Liver NADP SNPs")
  LiverNADP_Genes_MGI_chr8 <- query_genes_mgi(chr = chr, start = 60, end = 66)
  plot(out_snps_LiverNADP_chr8$lod, out_snps_LiverNADP_chr8$snpinfo, drop_hilit=1.5, genes = LiverNADP_Genes_MGI_chr8, main = "Liver NADP Genes MGI")
  
  

####################################################
## Liver NADPH
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverNADPH <- scan1(genoprobs = probs, pheno = pheno["zLiverNADPH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverNADPH <- scan1perm(genoprobs = probs, pheno = pheno["zLiverNADPH"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADPH = summary(perm_LiverNADPH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADPH", ylim = c(0,11))
  abline(h = threshold_LiverNADPH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksNADPH <- find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksNADPH <- find_peaks(scan1_output = qtlscan_LiverNADPH, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADPH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

# Liver NADPH --- Chromosome 12
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 12
  coef_blup_LiverNADPH_chr12 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADPH_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(5,15)
  plot_coefCC(x = coef_blup_LiverNADPH_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH, main = "Liver NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 12
  variants_LiverNADPH_chr12 <- query_variants(chr, 27, 30.5)
  out_snps_LiverNADPH_chr12 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                         chr = chr, start = 27, end = 30.5, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADPH_chr12$lod, out_snps_LiverNADPH_chr12$snpinfo, main = "Liver NADPH SNPs")
  
  LiverNADPH_Genes_MGI_chr12 <- query_genes_mgi(chr = chr, start = 27, end = 30.5)
  plot(out_snps_LiverNADPH_chr12$lod, out_snps_LiverNADPH_chr12$snpinfo, drop_hilit=1.5, genes = LiverNADPH_Genes_MGI_chr12, main = "Liver NADPH Genes MGI")
  
  

####################################################
## Liver NADP/NADPH Ratio
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverNADP_NADPHRatio <- scan1(genoprobs = probs, pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverNADP_NADPHRatio <- scan1perm(genoprobs = probs, pheno = pheno["zLiverNADP_NADPHRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_NADPHRatio = summary(perm_LiverNADP_NADPHRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP/NADPH Ratio", ylim = c(0,11))
  abline(h = threshold_LiverNADP_NADPHRatio, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksNADP_NADPHRatio <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  #using pmap (Mbp)
  peaksNADP_NADPHRatio <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPHRatio, map = R01_GSH_DO_QTLdata$pmap, threshold = summary(perm_LiverNADP_NADPHRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

# Liver NADP/NADPH Ratio --- Chromosome 12
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 12
  coef_blup_LiverNADP_NADPHRatio_chr12 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(5,15)
  plot_coefCC(x = coef_blup_LiverNADP_NADPHRatio_chr12, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPHRatio, main = "Liver NADP/NADPH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 12
  variants_LiverNADP_NADPHRatio_chr12 <- query_variants(chr, 27, 30.5)
  out_snps_LiverNADP_NADPHRatio_chr12 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                   chr = chr, start = 27, end =30.5, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPHRatio_chr12$lod, out_snps_LiverNADP_NADPHRatio_chr12$snpinfo, main = "Liver NADP/NADPH SNPs")
  LiverNADP_NADPHRatio_Genes_MGI_chr12 <- query_genes_mgi(chr = chr, start = 27, end = 30.5)
  plot(out_snps_LiverNADP_NADPHRatio_chr12$lod, out_snps_LiverNADP_NADPHRatio_chr12$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPHRatio_Genes_MGI_chr12, main = "Liver NADP/NADPH Genes MGI")

  
  
####################################################
## Export all QTL with LOD scores > 6 and all genes in QTL intervals
####################################################

qtl_gmap <- find_peaks(scans, map = R01_GSH_DO_QTLdata$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
qtl_gmap

write_xlsx(list(  "GSH chr14" = LiverGSH_Genes_MGI_chr14, 
                  "Total GSH chr14" = LiverTotalGSH_Genes_MGI_chr14, 
                  "GSH GSSG Ratio chr16" = LiverGSH_GSSGRatio_Genes_MGI_chr16), 
           "GlutathioneGenesMGI - rankZ sexgen.xlsx")

write_xlsx(list(  "NADP chr3" = LiverNADP_Genes_MGI_chr3,
                  "NADP chr8" = LiverNADP_Genes_MGI_chr8,
                  "NADPH chr12" = LiverNADPH_Genes_MGI_chr12,
                  "NADP NADPH Ratio chr12" = LiverNADP_NADPHRatio_Genes_MGI_chr12),
           "NADSystemsGenesMGI - rankZ sexgen.xlsx")



####################################################
## Export all QTL with LOD scores > 6 and all genes in QTL intervals
####################################################

qtl_gmap <- find_peaks(scans, map = R01_GSH_DO_QTLdata$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
qtl_pmap <- find_peaks(scans, map = R01_GSH_DO_QTLdata$pmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

qtl_gmap$marker.id <- find_marker(map = R01_GSH_DO_QTLdata$gmap, chr = qtl_gmap$chr, pos = qtl_gmap$pos)
qtl_pmap$marker.id <- find_marker(map = R01_GSH_DO_QTLdata$pmap, chr = qtl_pmap$chr, pos = qtl_pmap$pos)


write_xlsx(list("QTL List RankZ SexGen - cM" = qtl_gmap,
                "QTL List RankZ SexGen - Mbp" = qtl_pmap),
           "QTL List.xlsx")




