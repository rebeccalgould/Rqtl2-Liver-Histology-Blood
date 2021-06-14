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

  control <- read_cross2(file = "~/Rqtl2-Liver-Histology-Blood/data/control.json")


  
####################################################
## Genotype probabilities and allele probabilities - calculated by the Jackson Laboratory
####################################################
  
  probs <- readRDS("~/Rqtl2-Liver-Histology-Blood/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")


  
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
  pheno$zAST = rankZ(pheno$AST)
  pheno$zALT = rankZ(pheno$ALT)
  pheno$zASTALTRatio = rankZ(pheno$AST_ALT_Ratio)
  pheno$zSteatosis = rankZ(pheno$Liver_Steatosis)
  pheno$zBallooning = rankZ(pheno$Liver_Ballooning)
  pheno$zFibrosis = rankZ(pheno$Liver_Fibrosis)
  pheno$zInflammation = rankZ(pheno$Liver_Inflammation)
  pheno$zLiverWeight = rankZ(pheno$Liver_Weight)
  pheno$zBodyWeight = rankZ(pheno$Final_Weight)
  pheno$zLiverWeightBodyWeight = rankZ(pheno$Liver_Weight_Final_Weight_Ratio)


  
####################################################
## Add covariates
####################################################

  #adding sex and generation as covariates
  sexgen = model.matrix(~ sex + generation, data = pheno)[,-1]

  

####################################################
## ALT
## Plot Genome Scans with Permutation Tests
####################################################

  qtlscan_ALT <- scan1(genoprobs = probs, pheno = pheno["zALT"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_ALT <- scan1perm(genoprobs = probs, pheno = pheno["zALT"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_ALT = summary(perm_ALT, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_ALT, map = control$gmap,  main = "Genome Scan for ALT", ylim = c(0,11))
  abline(h = threshold_ALT, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksALT <- find_peaks(scan1_output = qtlscan_ALT, map = control$gmap, threshold = summary(perm_ALT, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksALT <- find_peaks(scan1_output = qtlscan_ALT, map = control$pmap, threshold = summary(perm_ALT, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
       
       
####################################################
## AST
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_AST <- scan1(genoprobs = probs, pheno = pheno["zAST"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_AST <- scan1perm(genoprobs = probs, pheno = pheno["zAST"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_AST = summary(perm_AST, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_AST, map = control$gmap,  main = "Genome Scan for AST", ylim = c(0,11))
  abline(h = threshold_AST, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksAST <- find_peaks(scan1_output = qtlscan_AST, map = control$gmap, threshold = summary(perm_AST, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksAST <- find_peaks(scan1_output = qtlscan_AST, map = control$pmap, threshold = summary(perm_AST, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
# AST --- Chromosome 2 
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 2
  coef_blup_AST_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores=2)
  plot_coefCC(x = coef_blup_AST_chr2, map = control$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(1,20)
  plot_coefCC(x = coef_blup_AST_chr2, map = control$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 2
  start = peaksAST[peaksAST$chr ==  chr,"ci_lo"]
  end = peaksAST[peaksAST$chr == chr, "ci_hi"]
  variants_AST_chr2 <- query_variants(chr, start - 1, end + 1)
  out_snps_AST_chr2 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                 chr = chr, start = start - 1, end = end + 1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_AST_chr2$lod, out_snps_AST_chr2$snpinfo, main = "AST SNPs")
  AST_Genes_MGI_chr2 <- query_genes_mgi(chr = chr, start = start - 1, end = end + 1)
  plot(out_snps_AST_chr2$lod, out_snps_AST_chr2$snpinfo, drop_hilit=1.5, genes = AST_Genes_MGI_chr2, main = "AST Genes MGI")
  
# AST --- Chromosome 16 
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 16
  coef_blup_AST_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores=2)
  plot_coefCC(x = coef_blup_AST_chr16, map = control$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(25,45)
  plot_coefCC(x = coef_blup_AST_chr16, map = control$gmap, scan1_output = qtlscan_AST, main = "AST BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 2
  start = peaksAST[peaksAST$chr ==  chr,"ci_lo"]
  end = peaksAST[peaksAST$chr == chr, "ci_hi"]
  variants_AST_chr16 <- query_variants(chr, start -1, end + 1)
  out_snps_AST_chr16 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zAST"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                  chr = chr, start = start - 1, end = end + 1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_AST_chr16$lod, out_snps_AST_chr16$snpinfo, main = "AST SNPs")
  
  AST_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = start - 1, end = end + 1)
  plot(out_snps_AST_chr16$lod, out_snps_AST_chr16$snpinfo, drop_hilit=1.5, genes = AST_Genes_MGI_chr16, main = "AST Genes MGI")
  
  
####################################################
## AST/ALT Ratio
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_ASTALTRatio <- scan1(genoprobs = probs, pheno = pheno["zASTALTRatio"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_ASTALTRatio <- scan1perm(genoprobs = probs, pheno = pheno["zASTALTRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_ASTALTRatio = summary(perm_ASTALTRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_ASTALTRatio, map = control$gmap,  main = "Genome Scan for AST/ALT Ratio", ylim = c(0,11))
  abline(h = threshold_ASTALTRatio, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksASTALTRatio <- find_peaks(scan1_output = qtlscan_ASTALTRatio, map = control$gmap, threshold = summary(perm_ASTALTRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksASTALTRatio <- find_peaks(scan1_output = qtlscan_ASTALTRatio, map = control$pmap, threshold = summary(perm_ASTALTRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  
  
####################################################
## Liver Steatosis
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_Steatosis <- scan1(genoprobs = probs, pheno = pheno["zSteatosis"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_Steatosis <- scan1perm(genoprobs = probs, pheno = pheno["zSteatosis"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  Xcovar = get_x_covar(control)
  perm_strata = mat2strata(Xcovar)
  perm_X_Steatosis <- scan1perm(genoprobs = probs, pheno = pheno["zSteatosis"], addcovar = sexgen, n_perm = 1000, perm_Xsp = TRUE, perm_strata = perm_strata, chr_lengths = chr_lengths(control$gmap), cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Steatosis = summary(perm_Steatosis, alpha = c(0.2, 0.1, 0.05))
  threshold_X_Steatosis = summary(perm_X_Steatosis, alpha = c(0.2, 0.1, 0.05))
  
  plot_scan1(x = qtlscan_Steatosis, map = control$gmap,  main = "Genome Scan for Steatosis (Autosome vs X)", ylim = c(0,11))
  segments(x0 = 0, y0 = threshold_X_Steatosis$A, x1 = 1695, y1 =   threshold_X_Steatosis$A, col = c("purple", "red", "blue"), "dashed", lwd = 2)
  segments(x0 = 1695, y0 = threshold_X_Steatosis$X, x1 = 2000, y1 = threshold_X_Steatosis$X, col = c("purple", "red", "blue"), "dashed", lwd = 2)
  
  #using gmap (cM)
  gmap_peaksSteatosis <- find_peaks(scan1_output = qtlscan_Steatosis, map = control$gmap, threshold = summary(perm_Steatosis, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksSteatosis <- find_peaks(scan1_output = qtlscan_Steatosis, map = control$pmap, threshold = summary(perm_Steatosis, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
# Steatosis --- Chromosome 2 
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 18
  coef_blup_Steatosis_chr18 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_Steatosis_chr18, map = control$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,20)
  plot_coefCC(x = coef_blup_Steatosis_chr18, map = control$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 18
  variants_Steatosis_chr18 <- query_variants(chr, 16, 18)
  out_snps_Steatosis_chr18 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                        chr = chr, start = 16, end = 18, keep_all_snps = TRUE)
  plot_snpasso(out_snps_Steatosis_chr18$lod, out_snps_Steatosis_chr18$snpinfo, main = "Steatosis SNPs")
  
  Steatosis_Genes_MGI_chr18 <- query_genes_mgi(chr = chr, start = 16, end = 18)
  plot(out_snps_Steatosis_chr18$lod, out_snps_Steatosis_chr18$snpinfo, drop_hilit=1.5, genes = Steatosis_Genes_MGI_chr18, main = "Steatosis Genes MGI")
  
# Steatosis --- Chromosome X
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "X"
  coef_blup_Steatosis_chrX <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_Steatosis_chrX, map = control$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(10,35)
  plot_coefCC(x = coef_blup_Steatosis_chrX, map = control$gmap, scan1_output = qtlscan_Steatosis, main = "Steatosis BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "X"
  variants_Steatosis_chrX <- query_variants(chr, 46, 53)
  out_snps_Steatosis_chrX <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zSteatosis"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = "X", start = 46, end = 53, keep_all_snps = TRUE)
  plot_snpasso(out_snps_Steatosis_chrX$lod, out_snps_Steatosis_chrX$snpinfo, main = "Steatosis SNPs")
  
  Steatosis_Genes_MGI_chrX <- query_genes_mgi(chr = chr, start = 46, end = 53)
  plot(out_snps_Steatosis_chrX$lod, out_snps_Steatosis_chrX$snpinfo, drop_hilit=1.5, genes = Steatosis_Genes_MGI_chrX, main = "Steatosis Genes MGI")
  
  

####################################################
## Liver Hydropic Degeneration (Ballooning)
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_Ballooning <- scan1(genoprobs = probs, pheno = pheno["zBallooning"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_Ballooning <- scan1perm(genoprobs = probs, pheno = pheno["zBallooning"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Ballooning = summary(perm_Ballooning, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_Ballooning, map = control$gmap,  main = "Genome Scan for Ballooning", ylim = c(0,11))
  abline(h = threshold_Ballooning, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksBallooning <- find_peaks(scan1_output = qtlscan_Ballooning, map = control$gmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksBallooning <- find_peaks(scan1_output = qtlscan_Ballooning, map = control$pmap, threshold = summary(perm_Ballooning, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

  
  
####################################################
## Inflammation
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_Inflammation <- scan1(genoprobs = probs, pheno = pheno["zInflammation"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_Inflammation <- scan1perm(genoprobs = probs, pheno = pheno["zInflammation"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_Inflammation = summary(perm_Inflammation, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_Inflammation, map = control$gmap,  main = "Genome Scan for Inflammation", ylim = c(0,11))
  abline(h = threshold_Inflammation, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksInflammation <- find_peaks(scan1_output = qtlscan_Inflammation, map = control$gmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksInflammation <- find_peaks(scan1_output = qtlscan_Inflammation, map = control$pmap, threshold = summary(perm_Inflammation, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  


####################################################
## Liver Weight
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverWeight <- scan1(genoprobs = probs, pheno = pheno["zLiverWeight"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverWeight <- scan1perm(genoprobs = probs, pheno = pheno["zLiverWeight"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverWeight = summary(perm_LiverWeight, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverWeight, map = control$gmap,  main = "Genome Scan for Liver Weight", ylim = c(0,11))
  abline(h = threshold_LiverWeight, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksLiverWeight <- find_peaks(scan1_output = qtlscan_LiverWeight, map = control$gmap, threshold = summary(perm_LiverWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksLiverWeight <- find_peaks(scan1_output = qtlscan_LiverWeight, map = control$pmap, threshold = summary(perm_LiverWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

# Liver Weight --- Chromosome 11
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "11"
  coef_blup_LiverWeight_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverWeight_chr11, map = control$gmap, scan1_output = qtlscan_LiverWeight, main = "Liver Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(15,35)
  plot_coefCC(x = coef_blup_LiverWeight_chr11, map = control$gmap, scan1_output = qtlscan_LiverWeight, main = "Liver Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "11"
  variants_LiverWeight_chr11 <- query_variants(chr,37.5, 41)
  out_snps_LiverWeight_chr11 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                          chr = chr, start = 37.5, end = 41, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverWeight_chr11$lod, out_snps_LiverWeight_chr11$snpinfo, main = "Liver Weight SNPs")
  
  LiverWeight_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = 37.5, end = 41)
  plot(out_snps_LiverWeight_chr11$lod, out_snps_LiverWeight_chr11$snpinfo, drop_hilit=1.5, genes = LiverWeight_Genes_MGI_chr11, main = "Liver Weight Genes MGI")

  
  
####################################################
## Body Weight
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_BodyWeight <- scan1(genoprobs = probs, pheno = pheno["zBodyWeight"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_BodyWeight <- scan1perm(genoprobs = probs, pheno = pheno["zBodyWeight"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_BodyWeight = summary(perm_BodyWeight, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_BodyWeight, map = control$gmap,  main = "Genome Scan for Body Weight", ylim = c(0,11))
  abline(h = threshold_BodyWeight, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksBodyWeight <- find_peaks(scan1_output = qtlscan_BodyWeight, map = control$gmap, threshold = summary(perm_BodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksBodyWeight <- find_peaks(scan1_output = qtlscan_BodyWeight, map = control$pmap, threshold = summary(perm_BodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  

  
####################################################
## Liver Weight/Body Weight Ratio
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_LiverWeightBodyWeight <- scan1(genoprobs = probs, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  perm_LiverWeightBodyWeight <- scan1perm(genoprobs = probs, pheno = pheno["zLiverWeightBodyWeight"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverWeightBodyWeight = summary(perm_LiverWeightBodyWeight, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_LiverWeightBodyWeight, map = control$gmap,  main = "Genome Scan for Liver Weight/Body Weight", ylim = c(0,11))
  abline(h = threshold_LiverWeightBodyWeight, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksLiverWeightBodyWeight <- find_peaks(scan1_output = qtlscan_LiverWeightBodyWeight, map = control$gmap, threshold = summary(perm_LiverWeightBodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksLiverWeightBodyWeight <- find_peaks(scan1_output = qtlscan_LiverWeightBodyWeight, map = control$pmap, threshold = summary(perm_LiverWeightBodyWeight, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
# Liver Weight/Body Weight --- Chromosome 11
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "11"
  coef_blup_LiverWeightBodyWeight_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr11, map = control$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(15,40)
  plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr11, map = control$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "11"
  variants_LiverWeightBodyWeight_chr11 <- query_variants(chr, 54.5, 58)
  out_snps_LiverWeightBodyWeight_chr11 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                    chr = chr, start = 54.5, end = 58, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverWeightBodyWeight_chr11$lod, out_snps_LiverWeightBodyWeight_chr11$snpinfo, main = "Liver Weight/Body Weight SNPs")
  
  LiverWeightBodyWeight_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = 54.5, end = 58)
  plot(out_snps_LiverWeightBodyWeight_chr11$lod, out_snps_LiverWeightBodyWeight_chr11$snpinfo, drop_hilit=1.5, genes = LiverWeightBodyWeight_Genes_MGI_chr11, main = "Liver Weight/Body Weight Genes MGI")

# Liver Weight/Body Weight --- Chromosome 12
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "12"
  coef_blup_LiverWeightBodyWeight_chr12 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr12, map = control$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(50,70)  
  plot_coefCC(x = coef_blup_LiverWeightBodyWeight_chr12, map = control$gmap, scan1_output = qtlscan_LiverWeightBodyWeight, main = "Liver Weight/Body Weight BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "12"
  variants_LiverWeightBodyWeight_chr12 <- query_variants(chr, 111.5, 114)
  out_snps_LiverWeightBodyWeight_chr12 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["zLiverWeightBodyWeight"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                    chr = chr, start = 111.5, end = 114, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverWeightBodyWeight_chr12$lod, out_snps_LiverWeightBodyWeight_chr12$snpinfo, main = "Liver Weight/Body Weight SNPs")
  
  LiverWeightBodyWeight_Genes_MGI_chr12 <- query_genes_mgi(chr = chr, start = 111.5, end = 114)
  plot(out_snps_LiverWeightBodyWeight_chr12$lod, out_snps_LiverWeightBodyWeight_chr12$snpinfo, drop_hilit=1.5, genes = LiverWeightBodyWeight_Genes_MGI_chr12, main = "Liver Weight/Body Weight Genes MGI")
  
  
  
####################################################
## Export genes
####################################################

write_xlsx(list(  "AST chr2" = AST_Genes_MGI_chr2,
                  "AST chr16" = AST_Genes_MGI_chr16),
           "BloodValuesGenesMGI - rankZ sexgen.xlsx")

write_xlsx(list(  "Steatosis chr18" = Steatosis_Genes_MGI_chr18,
                  "Steatosis chrX" = Steatosis_Genes_MGI_chrX),
           "HistologyGenesMGI - rankZ sexgen.xlsx")


####################################################
## Export all QTL with LOD scores > 6 
####################################################

scans <- cbind(qtlscan_AST, qtlscan_ALT, qtlscan_ASTALTRatio, qtlscan_Ballooning, qtlscan_Steatosis, qtlscan_Inflammation, qtlscan_LiverWeight, qtlscan_BodyWeight, qtlscan_LiverWeightBodyWeight)

qtl_gmap <- find_peaks(scans, map = control$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
qtl_pmap <- find_peaks(scans, map = control$pmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

qtl_gmap$marker.id <- find_marker(map = control$gmap, chr = qtl_gmap$chr, pos = qtl_gmap$pos)
qtl_pmap$marker.id <- find_marker(map = control$pmap, chr = qtl_pmap$chr, pos = qtl_pmap$pos)


write_xlsx(list("QTL List RankZ SexGen - cM" = qtl_gmap,
                "QTL List RankZ SexGen - Mbp" = qtl_pmap),
           "QTL List.xlsx")


####################################################
## Heritability Calculations
####################################################

kinship_lmm <- calc_kinship(probs = probs, use_allele_probs = TRUE, cores = 2)

herit_AST <- est_herit(pheno["zAST"], kinship_lmm, sexgen, cores = 2)
herit_ALT <- est_herit(pheno["zALT"], kinship_lmm, sexgen, cores = 2)
herit_ASTALTRatio <- est_herit(pheno["zASTALTRatio"], kinship_lmm, sexgen, cores = 2)
herit_Steatosis <- est_herit(pheno["zSteatosis"], kinship_lmm, sexgen, cores = 2)
herit_Ballooning <- est_herit(pheno["zBallooning"], kinship_lmm, sexgen, cores = 2)
herit_Fibrosis <- est_herit(pheno["zFibrosis"], kinship_lmm, sexgen, cores = 2)
herit_Inflammation <- est_herit(pheno["zInflammation"], kinship_lmm, sexgen, cores = 2)
herit_LiverWeight <- est_herit(pheno["zLiverWeight"], kinship_lmm, sexgen, cores = 2)
herit_BodyWeight <- est_herit(pheno["zBodyWeight"], kinship_lmm, sexgen, cores = 2)
herit_LiverWeightBodyWeight <- est_herit(pheno["zLiverWeightBodyWeight"], kinship_lmm, sexgen, cores = 2)


