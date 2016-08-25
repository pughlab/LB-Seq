############################################################################
#
#                     ctDNA Indel Data Filtering
#
# Author: Olena Kis
# Last Modified: Aug 25, 2016
# Date Created: June 3, 2015

###############################
#         Variables
###############################

project.name <- "MYELSTONE_ctDNA"
setwd("/Users/okis/Data_Analysis/Indels/MYELSTONE/analysis/ctDNA") 
file.path <- "/Users/okis/Data_Analysis/Indels/MYELSTONE/maf/ctDNA" 
file.ext <- ".maf"

blacklist <- "/Users/okis/Data_Analysis/file_path/pancan_mutation_blacklist.v10.hg19.txt"

###############################
#         Arguments           #
###############################

# Set the threshold for low coverage filter (minumum depth of coverage required to keep the data)
Cov.thresh <- 5000

# Set the threshold for low confidance skewed indels calls (minumum number of calls required 
# in the forward and reverse direction to support the indel)
Skew.thresh <- 10

summary_vect <- c("Matched_Norm_Sample_Barcode",
                  "Hugo_Symbol",
                  "Genome_Change",
                  "Variant_Type",
                  "Variant_Classification",
                  "Chromosome",
                  "Start_position",
                  "End_position",
                  "Reference_Allele",
                  "Tumor_Seq_Allele1",
                  "Tumor_Seq_Allele2",
                  "CGC_Chr.Band",
                  "AC_consensus",
                  "AC_any",
                  "tumor_f",
                  "tumor_f_any",
                  "depth_across_samples",
                  "indel_forward",
                  "indel_reverse",
                  "ref_forward",
                  "ref_reverse",
                  "dbSNP_Val_Status",
                  "dbSNP_RS",
                  "COSMIC_overlapping_mutations")

###############################
#         Functions           #
###############################

# Data summary for filtered maf results
summary_function <- function(all_ctDNA_indels,summary_vect){
  summary_ctDNA_indels<- all_ctDNA_indels[,summary_vect]
  return(summary_ctDNA_indels) 
}


###############################
#           Main              #
###############################


# Combined all maf files from ctDNA data directory into one file
filenames <- list.files(path = file.path, pattern = file.ext,
                        full.names = T, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)
all_ctDNA_indels <- NULL

for (file in filenames){
  ctDNA_indels <- read.table(file, header = TRUE, sep = "\t", quote = "",
                          comment.char = "#", stringsAsFactors =FALSE)
  name <- basename(file)
  name <- gsub(".bam.vcf.maf","",name)  
  name <- gsub("_ctDNA","",name) 
  name <- gsub("-ctDNA","",name) 
  all_ctDNA_indels <- rbind(all_ctDNA_indels, ctDNA_indels) 
}

## Add 'Allele Fraction' and "tumor_f to indel file
allele_count <-unlist(strsplit(all_ctDNA_indels$allele_count,","))
allele_count_consensus <-allele_count[seq(from=1,to=length(allele_count),by=2)]
allele_count_any <-allele_count[seq(from=2,to=length(allele_count),by=2)]
all_ctDNA_indels$AC_consensus <- as.numeric(allele_count_consensus)
all_ctDNA_indels$AC_any <- as.numeric(allele_count_any)
all_ctDNA_indels$tumor_f <-all_ctDNA_indels$AC_consensus/all_ctDNA_indels$depth_across_samples
all_ctDNA_indels$tumor_f_any  <-all_ctDNA_indels$AC_any/all_ctDNA_indels$depth_across_samples

## Add forward and reverse reads for Indel-supporti and reference

supporting_reads <-unlist(strsplit(all_ctDNA_indels$SC,","))
indel_forward <-supporting_reads[seq(from=1,to=length(supporting_reads),by=4)]
indel_reverse <-supporting_reads[seq(from=2,to=length(supporting_reads),by=4)]
ref_forward <-supporting_reads[seq(from=3,to=length(supporting_reads),by=4)]
ref_reverse <-supporting_reads[seq(from=4,to=length(supporting_reads),by=4)]
all_ctDNA_indels$indel_forward <- as.numeric(indel_forward)
all_ctDNA_indels$indel_reverse <- as.numeric(indel_reverse)
all_ctDNA_indels$ref_forward <- as.numeric(ref_forward)
all_ctDNA_indels$ref_reverse <- as.numeric(ref_reverse)

write.table(all_ctDNA_indels, file = paste("All_indels",project.name,"txt",sep="."),
            row.names=FALSE, append = FALSE,na = "NA", quote = FALSE, sep = "\t", col.names = TRUE)

summary_ctDNA_indels <- summary_function(all_ctDNA_indels, summary_vect)

write.table(summary_ctDNA_indels, file = paste("Indels_summary",project.name,"txt",sep="."),
            row.names=FALSE, append = FALSE,na = "NA", quote = FALSE, sep = "\t", col.names = TRUE)

# Filter out validated indels

filtered_ctDNA_indels <- all_ctDNA_indels[all_ctDNA_indels$dbSNP_Val_Status == "",]

nrow(all_ctDNA_indels)
nrow(filtered_ctDNA_indels)

# Filter out low coverage data

filtered_ctDNA_indels <- filtered_ctDNA_indels[filtered_ctDNA_indels$depth_across_samples >= Cov.thresh,]

# Filter out skewed data (if have < 10 reads from forward or reverse strand)

filtered_ctDNA_indels <- filtered_ctDNA_indels[filtered_ctDNA_indels$indel_forward >= Skew.thresh,]
filtered_ctDNA_indels <- filtered_ctDNA_indels[filtered_ctDNA_indels$indel_reverse >= Skew.thresh,]


PanCan_mutations <- read.table(blacklist, header = TRUE,  
                               sep = "\t",quote = "", comment.char = "#", stringsAsFactors =FALSE)
PanCan_mutations$locus <- paste("g.chr",PanCan_mutations$chr,":",PanCan_mutations$start,
                                PanCan_mutations$ref_allele,">",PanCan_mutations$newbase, sep="")
ctDNA_filter_PanCan <- PanCan_mutations$locus
filtered_ctDNA_indels <- filtered_ctDNA_indels[!filtered_ctDNA_indels$Genome_Change %in% ctDNA_filter_PanCan,]

summary_filtered_ctDNA_indels <- summary_function(filtered_ctDNA_indels, summary_vect)

write.table(summary_filtered_ctDNA_indels, file = paste("Summary_filtered_indels",project.name,"txt",sep="."),
            row.names=FALSE, append = FALSE,na = "NA", quote = FALSE, sep = "\t", col.names = TRUE)