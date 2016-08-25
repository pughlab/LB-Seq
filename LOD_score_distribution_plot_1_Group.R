###############################################################################################
#
#                     ctDNA Data Filtering - MEGAPLOT of LOD scores in all samples
#
# Author: Olena Kis
# Last Modified: July 6, 2016
# Date Created: May 26, 2015
###############################################################################################


###############################
#         Variables
###############################

group.name <- "1GROUP"

CWD <- "/Users/okis/Data_Analysis/No_DS_MYL/LOD_distribution_plot/analysis/2Groups"
setwd(CWD)

file.path <- "/Users/okis/Data_Analysis/No_DS_MYL/LOD_distribution_plot/maf/Training"
file.ext <- ".maf"


# Strand bias filter
    LODfwd_rev_ratio_thr <- 5  # highest allowable fold difference between forward and reverse LOD scores
# Tumor LOD score filter
    LOD_Zscore_thresh <- 20  # determined using ROC_curve_LOD_score.R script
# Threshold allele frequency in 1000 Genomes data used to call a variant as a germline SNP 
    X1000Gen_thresh <- 0.001  ## Filters based on >= 0.1% frequency in normal population 
# Specify the minimum fraction of total reads which has to be used as supporting reads (Alt + Ref count) to call a mutations
    Fraction_supp_reads <- 0.1 # the sum of ref_allele_counts and alt_allele counts should be at least 10% of total reads

    summary_vect <- c("Sample_ID",
                      "Genome_Change",
                      "Hugo_Symbol",
                      "Protein_Change",
                      "UniProt_AApos",
                      "Prot_Change_UniProt",
                      "Variant_Classification",
                      "cDNA_Change",
                      "tumor_f",
                      "t_lod_fstar",
                      "modified_Z_score",
                      "t_alt_count",
                      "t_ref_count",
                      "total_pairs",
                      "fraction_supp_reads",
                      "judgement",
                      "t_lod_fstar_forward",
                      "t_lod_fstar_reverse",
                      "skew1",
                      "skew2",
                      "Tumor_Sample_Barcode",
                      "dbSNP_Val_Status",
                      "dbSNP_RS",
                      "X1000gp3_AF",
                      "X1000gp3_AMR_AF",
                      "X1000gp3_AFR_AF",
                      "X1000gp3_EUR_AF",
                      "X1000gp3_EAS_AF",
                      "X1000gp3_SAS_AF")
    
###############################
#         Functions           #
###############################

# Function to filter mutations with strand bias (LOD fwd / LOD rev >= |5|)
    strand_bias_function <- function(filtered_ctDNA_maf){
      filtered_ctDNA_maf <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew1 < LODfwd_rev_ratio_thr,]
      filtered_ctDNA_maf <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew1 >= 0,]
      filtered_ctDNA_maf <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew2 < LODfwd_rev_ratio_thr,]
      filtered_ctDNA_maf <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew2 >= 0,]
      filtered_ctDNA_maf <- filtered_ctDNA_maf[complete.cases(filtered_ctDNA_maf$Start_position),] #removes NA data
      return(filtered_ctDNA_maf)
    }
    strand_bias_function_rm <- function(filtered_ctDNA_maf){
      removed_ctDNA_maf_1 <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew1 >= LODfwd_rev_ratio_thr,]
      removed_ctDNA_maf_2 <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew1 < 0,]
      removed_ctDNA_maf_3 <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew2 >= LODfwd_rev_ratio_thr,]
      removed_ctDNA_maf_4 <- filtered_ctDNA_maf[filtered_ctDNA_maf$skew2 < 0,]
      removed_ctDNA_maf_strand_bias <- rbind(removed_ctDNA_maf_1, removed_ctDNA_maf_2, removed_ctDNA_maf_3, 
                                             removed_ctDNA_maf_4)
      removed_ctDNA_maf_strand_bias <- unique(removed_ctDNA_maf_strand_bias)
      return(removed_ctDNA_maf_strand_bias)
    }
    
# Function for removing germline SNPs based on 1000Genomes Project allele frequency in specific ethnic populations
    ctDNA_filter_SNP_eth <- function(filtered_ctDNA_maf){
      SNP_maf1 <- subset(filtered_ctDNA_maf, as.numeric(as.character(X1000gp3_AF)) >= X1000Gen_thresh)
      SNP_maf2 <- subset(filtered_ctDNA_maf, as.numeric(as.character(X1000gp3_AFR_AF)) >= X1000Gen_thresh)
      SNP_maf3  <- subset(filtered_ctDNA_maf, as.numeric(as.character(X1000gp3_AMR_AF)) >= X1000Gen_thresh)
      SNP_maf4  <- subset(filtered_ctDNA_maf, as.numeric(as.character(X1000gp3_EUR_AF)) >= X1000Gen_thresh)
      SNP_maf5  <- subset(filtered_ctDNA_maf, as.numeric(as.character(X1000gp3_EAS_AF)) >= X1000Gen_thresh)
      SNP_maf6  <- subset(filtered_ctDNA_maf, as.numeric(as.character(X1000gp3_SAS_AF)) >= X1000Gen_thresh)
      ethnic_SNP_maf <- rbind(SNP_maf1,SNP_maf2,SNP_maf3,SNP_maf4,SNP_maf5,SNP_maf6)
      ethnic_SNP_maf <- unique(ethnic_SNP_maf)
      return(ethnic_SNP_maf)
    }
    

###############################
#           Main
###############################

#1. TRAINING cohort with BM

# Combined all maf files from ctDNA data directory into one file
filenames <- list.files(path = file.path, pattern = file.ext,
                        full.names = T, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)
    all_ctDNA_maf <- NULL
    all_pop_data <- NULL      
    
    for (file in filenames){
      ctDNA_maf <- read.table(file, header = TRUE, sep = "\t", quote = "",
                              comment.char = "#", stringsAsFactors =FALSE)
      ## remove all calls REJECTED by MuTect as likely sequencer errors
      ## (only needed if REJECT mutations were included in the MuTect output file)
      ctDNA_maf <- ctDNA_maf[ctDNA_maf$judgement == "KEEP",]
      ## remove mutation calls made using a small fraction of total reads
      ctDNA_maf$total_supp_reads <- ctDNA_maf$t_alt_count + ctDNA_maf$t_ref_count
      ctDNA_maf$fraction_supp_reads <- ctDNA_maf$total_supp_reads / ctDNA_maf$total_pairs
      ctDNA_maf <- ctDNA_maf[ctDNA_maf$fraction_supp_reads > Fraction_supp_reads,]
      ## simplify sample name to study subject ID (MYL-XXX)
      name <- basename(file)
      name <- gsub(".call_stats.maf","",name)
      name <- gsub("-01.bam","",name)    #*** change based on .maf file names ***
      name <- gsub("-02.bam","",name)
      name <- gsub("-03.bam","",name)
      name <- gsub("-04.bam","",name)
      name <- gsub("-14.bam","",name)
      name <- gsub("-4.bam","",name)
      name <- gsub(".bam","",name)
      name <- gsub(".processed","",name)
      name <- gsub("myl","MYL",name)
      name <- gsub("_Tumor","",name)     #*** change based on .maf file names ***
      name <- gsub("_hyb2","",name)
      name <- gsub("_ctDNA","",name)
      name <- gsub("-ctDNA","",name)
      name <- gsub(" (83 ng)","",name)
      name <- gsub("-reseq","",name)
      name <- gsub("-BM","",name)
      name <- gsub("B1_","",name)
      name <- gsub("B2_","",name)        #*** change based on .maf file names ***
      name <- gsub("B3_","",name)
      name <- gsub("B4_","",name)
      name <- gsub("B5_","",name)
      name <- gsub("B6_","",name)
      name <- gsub("B7_","",name)
      name <- gsub("B8_","",name)
      name <- gsub("B9_","",name)
      name <- gsub("B10_","",name)
      ctDNA_maf$row_names <- paste(name,ctDNA_maf$Genome_Change,sep="_")
      ctDNA_maf$Sample_ID <- paste(name)
      ## Convert tumor LOD scores to modified Z-scores, the number of Median 
      ## Absolute Deviation (MADs) from the median LOD score in each sample
      LOD_scores <- ctDNA_maf$t_lod_fstar
      median_LOD <- as.numeric(median(LOD_scores))
      MAD_LOD  <- as.numeric(mad(LOD_scores))
      LOD_threshold <- as.numeric(median_LOD + LOD_Zscore_thresh*MAD_LOD)
      ctDNA_maf$modified_Z_score <- (ctDNA_maf$t_lod_fstar - median_LOD)/MAD_LOD
      total_mut <- as.numeric(nrow(ctDNA_maf))
      pop_data <- c(name, total_mut, median_LOD, MAD_LOD, LOD_threshold)
      all_pop_data <- rbind(all_pop_data, pop_data)
      all_ctDNA_maf <- rbind(all_ctDNA_maf, ctDNA_maf) 
    }
    
colnames(all_pop_data) <- paste(c("sample_ID", "Total_Mutations", "Median_LOD", "MAD", "LOD_Threshold"))
write.table(all_pop_data, file = paste("temp","txt",sep="."),
            row.names=FALSE, append = FALSE,na = "NA", quote = FALSE, sep = "\t", col.names = TRUE)
all_pop_data <- read.table("temp.txt", header = TRUE, sep = "\t", quote = "",
                          comment.char = "#", stringsAsFactors =FALSE)

# Create a "Protein_Change_UniProt" column for a different annotation (used by COSMIC) *** important for EGFR 

      substrRight <- function(x, n){
        substr(x, nchar(x)-n+1, nchar(x))
      }
      x <- all_ctDNA_maf$Protein_Change
      all_ctDNA_maf$Protein_Change_start <- substr(x, 1, 3) 
      all_ctDNA_maf$Protein_Change_end <- substrRight(x, 1)
      all_ctDNA_maf$Prot_Change_UniProt <- paste(all_ctDNA_maf$Protein_Change_start, all_ctDNA_maf$UniProt_AApos,
                                           all_ctDNA_maf$Protein_Change_end ,sep="")
            
# Add additional columns in the dataframe for downstreem filtering
      all_ctDNA_maf$skew1 <- all_ctDNA_maf$t_lod_fstar_forward / all_ctDNA_maf$t_lod_fstar_reverse
      all_ctDNA_maf$skew2 <- all_ctDNA_maf$t_lod_fstar_reverse / all_ctDNA_maf$t_lod_fstar_forward
      filtered_ctDNA_maf <- all_ctDNA_maf
      
# STEP 1: Remove mutations with strand bias (calls with X-fold difference between Forward and Reverse t_lod_fstar)
      removed_ctDNA_maf_strand_bias <- strand_bias_function_rm(all_ctDNA_maf)
      filtered_ctDNA_maf <- strand_bias_function(all_ctDNA_maf)
      count_filtered_maf_1 <- nrow(filtered_ctDNA_maf)
      count_removed_strand_bias <- nrow(removed_ctDNA_maf_strand_bias)
      
# STEP 2: Apply LOD score filter: keep only the mutations with t_lod_fstar >= THRESHOLD
      removed_ctDNA_maf_LOD <- filtered_ctDNA_maf[filtered_ctDNA_maf$modified_Z_score < LOD_Zscore_thresh,]
      filtered_ctDNA_maf <- filtered_ctDNA_maf[filtered_ctDNA_maf$modified_Z_score >= LOD_Zscore_thresh,]
      count_filtered_maf_2 <- nrow(filtered_ctDNA_maf)
      count_removed_LOD_thr <- nrow(removed_ctDNA_maf_LOD)
      
# STEP 3: Filtering validated germline SNPs using Oncotator version 1.5.3.0 annotation data
      removed_ctDNA_maf_valSNP <- filtered_ctDNA_maf[!filtered_ctDNA_maf$dbSNP_Val_Status == "",]
      removed_ctDNA_maf_valSNP <- subset(filtered_ctDNA_maf, as.numeric(as.character(X1000gp3_AF)) >= X1000Gen_thresh)
      val_SNPs_list <- removed_ctDNA_maf_valSNP$Genome_Change
      filtered_ctDNA_maf <- subset(filtered_ctDNA_maf, !(Genome_Change %in% val_SNPs_list))
      
# STEP 4: Filtering additional SNPs with 1000 Genomes data AF >= X1000Gen_thresh in specific ethnic groups
      ethnic_SNP_maf <- ctDNA_filter_SNP_eth(filtered_ctDNA_maf)
      ethnic_SNPs <- ethnic_SNP_maf$Genome_Change 
      filtered_ctDNA_maf <- subset(filtered_ctDNA_maf, !(Genome_Change %in% ethnic_SNPs))
      count_filtered_maf_3 <- nrow(filtered_ctDNA_maf)
      SNP_maf <- rbind(removed_ctDNA_maf_valSNP, ethnic_SNP_maf)
      count_removed_SNPs <- nrow(SNP_maf)
      
# Create filtering data summary
      filter_type    <- c('pre-filtering','Strand Bias (LOD fwd/rev rario)',
                       "Below LOD Threshold",'Germline SNPs')
      removed_count  <- c(0, count_removed_strand_bias,
                      count_removed_LOD_thr, count_removed_SNPs)
      remaining_count <- c(nrow(all_ctDNA_maf), count_filtered_maf_1,  
                           count_filtered_maf_2, count_filtered_maf_3)
      filtering_summary <- data.frame(filter_type, removed_count, remaining_count)

# Plot the distribution of LOD scores and Normalized LOD scores

positive_LOD_maf <- all_ctDNA_maf[all_ctDNA_maf$t_lod_fstar > 0,]
#positive_filtered_maf <- filtered_ctDNA_maf[filtered_ctDNA_maf$t_lod_fstar > 0,]
positive_strand_bias <- removed_ctDNA_maf_strand_bias[removed_ctDNA_maf_strand_bias$t_lod_fstar > 0,]

# Create a dataframe with the list of samples without any somatic mutations identified to have a list of samples without any mutations identified!!!
unique_samples <- NULL
without_mut_ctDNA_maf <- NULL
unique_samples <- all_ctDNA_maf[ !duplicated(all_ctDNA_maf$Sample_ID), ]
samples_with_mut <- filtered_ctDNA_maf$Sample_ID
without_mut_ctDNA_maf <- subset(unique_samples, !(Sample_ID %in% samples_with_mut))
without_mut_ctDNA_maf$t_lod_fstar <- 0
without_mut_ctDNA_maf$modified_Z_score <- 0

filtered_maf <- rbind(filtered_ctDNA_maf,without_mut_ctDNA_maf)
order.sample.ID = order(filtered_maf$Sample_ID)
filtered_maf = filtered_maf[order.sample.ID,]
order.sample.ID = order(SNP_maf$Sample_ID)
SNP_maf = SNP_maf[order.sample.ID,]
order.sample.ID = order(positive_strand_bias$Sample_ID)
strand_bias = positive_strand_bias[order.sample.ID,] 
order.sample.ID = order(removed_ctDNA_maf_LOD$Sample_ID)
unremoved_artifacts = removed_ctDNA_maf_LOD[order.sample.ID,] 
order.sample.ID = order(all_pop_data$sample_ID)
all_pop_data = all_pop_data[order.sample.ID,] 


##########################################
#     LOD Score Distribution Plot        #
##########################################

pdf(file = paste("LOD_distribution_plot",group.name,"pdf",sep="."), width = 11, height = 7)
range <- c(log10(min(positive_LOD_maf$t_lod_fstar)), log10(max(positive_LOD_maf$t_lod_fstar)))
stripchart(log10(strand_bias$t_lod_fstar) ~ strand_bias$Sample_ID, 
           plot = TRUE, vertical = TRUE, frame.plot = TRUE, ylab="Tumor LOD score", yaxt="n", 
           method = "overplot", ylim=range, cex = 1, cex.axis = 0.8, 
           par(las=2, mai = c(1.2, 0.9, 0.5, 0.3)), pch = 1, col="#BCEE68")
stripchart(log10(unremoved_artifacts$t_lod_fstar) ~ unremoved_artifacts$Sample_ID, 
           vertical = TRUE, ylab="Tumor LOD score", yaxt="n", method = "overplot", 
           frame.plot = TRUE, ylim=range, add = TRUE, pch = 1, col="#CDC8B1")
stripchart(log10(SNP_maf$t_lod_fstar) ~ SNP_maf$Sample_ID, vertical = TRUE, 
           ylab="Tumor LOD score", yaxt="n", method = "overplot", 
           frame.plot = TRUE, ylim=range, add = TRUE, pch = 16, col="red")
stripchart(log10(filtered_maf$t_lod_fstar) ~ filtered_maf$Sample_ID, 
           vertical = TRUE, ylab="Tumor LOD score", yaxt="n", method = "overplot", 
           frame.plot = TRUE, ylim=range, add = TRUE, pch = 16, col="blue")
bp <- boxplot(log10(positive_LOD_maf$t_lod_fstar) ~ positive_LOD_maf$Sample_ID,
              ylab="Tumor LOD score", yaxt="n",whisklty = 0, staplelty = 0, add = TRUE, xaxt="n", outline=FALSE,
              ylim=range, col="#ffffbf")
stripchart(log10(all_pop_data$LOD_Threshold) ~ all_pop_data$sample_ID, 
           vertical = TRUE, ylab="Tumor LOD score", yaxt="n", method = "overplot", 
           frame.plot = TRUE, ylim=range, add = TRUE, pch = "-", cex = 1.2, col="black")
axis(2, las=2, cex.axis = 0.8,
     at=log10(c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,
                0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
               200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,
                8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,
               200000,300000,400000,500000,600000,700000,800000,900000,1000000)),
     labels=c("0.01","","","","","","","","","0.1","","","","","","","","","1",
              "","","","","","","","","10","","","","","","","","","100",
              "","","","","","","","","1,000","","","","","","","","","10,000",
              "","","","","","","","","100,000","","","","","","","","","1,000,000"))
dev.off()




pdf(file = paste("LEGEND_and_Somatic_mutations",group.name,"pdf",sep="."), width = 10, height = 7)
range <- c(log10(min(positive_LOD_maf$t_lod_fstar)), log10(max(positive_LOD_maf$t_lod_fstar)))
stripchart(log10(filtered_maf$t_lod_fstar) ~ filtered_maf$Sample_ID, 
           plot = TRUE, vertical = TRUE, frame.plot = TRUE, ylab="Tumor LOD score", yaxt="n", 
           method = "overplot", ylim=range, cex = 1, cex.axis = 0.82, 
           par(las=2, mai = c(1.2, 1, 0.5, 0.5)), pch = 16, col="blue")
legend("bottomleft", pch = c(1,1), col= c("#BCEE68","#CDC8B1"),
       legend=c("sequence context-dependent artifacts", "likely artifacts (below LOD threshold)"), 
       cex = 0.85, pt.cex = 1.5, bty="n")
legend("bottom", pch = c(16,16), col= c("red","blue"),
       legend=c("germline polymorphisms","somatic mutation calls"), 
       cex = 0.85, pt.cex = 1.5, bty="n")
legend("bottomright", pch = "-", col="black",
       legend="sample-specific LOD threshold", 
       cex = 0.85, pt.cex = c(1.5,1), bty="n")
axis(2, las=2, cex.axis = 0.8,
     at=log10(c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,
                0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
                200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,
                8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,
                200000,300000,400000,500000,600000,700000,800000,900000,1000000)),
     labels=c("0.01","","","","","","","","","0.1","","","","","","","","","1",
              "","","","","","","","","10","","","","","","","","","100",
              "","","","","","","","","1,000","","","","","","","","","10,000",
              "","","","","","","","","100,000","","","","","","","","","1,000,000"))
dev.off()


##########################################
#          DATA SUMMARY TABLES           #
##########################################

write.table(all_pop_data, file = paste("LOD_distribution_data",group.name,"txt",sep="."),
            row.names=FALSE, append = FALSE,na = "NA", quote = FALSE, sep = "\t", col.names = TRUE)

write.table(filtering_summary, file = paste("ALL_filtering_summary",group.name,"txt",sep="."),
            row.names=FALSE, append = FALSE,na = "NA", quote = FALSE, sep = "\t", col.names = TRUE)

filtered_summary_maf <- filtered_maf[,summary_vect]
write.table(filtered_summary_maf, file = paste("ALL_somatic_mutations",group.name,"txt",sep="."),
            row.names=FALSE, append = FALSE, na = "NA", quote = FALSE, sep = "\t", col.names = TRUE)

file.remove(file = paste(CWD,"/temp",".txt",sep=""))
