
library(igraph)
library(stringr)
library(ggplot2)
options(scipen = 999)
library(matrixStats)
library(RColorBrewer)
library("plyr")
library(dplyr)
library(GenomicRanges)

distances <- c(100000,250000,500000,1000000)
num_permutations <- 1000

for (m in 1:length(distances)){
target <- c("Smarca4","Smarca5","Ino80","Chd4","Chd2","Chd3","SRCAP","Atrx","Hells","Chd1","Btaf1","Chd7","Chd8","Chd9","Chd1l","Chd5","Ercc6","Ercc6l","Hltf","Rad54b","Rad54l","Rad54l2","Shprh","SmarcaD","SmarcaL","Ttf2","Zranb3","Ep400","Smarca2")
PT_table <- NULL
for (i in 1:length(target)) {
rm(list=ls(pattern="mRNA_DHS_"))
#i=1
name <- target[i]

mRNAs_up <- read.delim(paste0("~/TTseq_DEseq_res/",name,"_promoters_up.bed"), header = F, sep = "\t", quote = "")
colnames(mRNAs_up)<- c("chrom","start","end","name","score","strand")
DHS_up <- read.delim(paste0("~/TTseq_DEseq_res/",name,"_DHSs_up.bed"), header = F, sep = "\t", quote = "")
colnames(DHS_up)<- c("chrom","start","end","name","score","strand")


gr_df1 <- GRanges(
  seqnames = Rle(mRNAs_up$chrom),
  ranges = IRanges(start = mRNAs_up$start, end = mRNAs_up$end), metadata = DataFrame(region_names = mRNAs_up$name))
gr_df2 <- GRanges(
  seqnames = Rle(DHS_up$chrom),
  ranges = IRanges(start = DHS_up$start, end = DHS_up$end), metadata = DataFrame(region_names = DHS_up$name))

overlaps_up <- findOverlaps(gr_df1, gr_df2, maxgap = distances[m])

overlapping_regions_gr1 <- as.data.frame(gr_df1[queryHits(overlaps_up)])
overlapping_regions_gr2 <- as.data.frame(gr_df2[subjectHits(overlaps_up)])

num_overlaps_up <- length(unique(queryHits(overlaps_up)))
#overlapping_regions_up <- as.data.frame(gr_df1[queryHits(overlaps_enhancers)])

mRNAs_down <- read.delim(paste0("~/TTseq_DEseq_res/",name,"_promoters_down.bed"), header = F, sep = "\t", quote = "")
colnames(mRNAs_down)<- c("chrom","start","end","name","score","strand")
DHS_down <- read.delim(paste0("~/TTseq_DEseq_res/",name,"_DHSs_down.bed"), header = F, sep = "\t", quote = "")
colnames(DHS_down)<- c("chrom","start","end","name","score","strand")


gr_df3 <- GRanges(
  seqnames = Rle(mRNAs_down$chrom),
  ranges = IRanges(start = mRNAs_down$start, end = mRNAs_down$end), metadata = DataFrame(region_names = mRNAs_down$name))
gr_df4 <- GRanges(
  seqnames = Rle(DHS_down$chrom),
  ranges = IRanges(start = DHS_down$start, end = DHS_down$end), metadata = DataFrame(region_names = DHS_down$name))

overlaps_down <- findOverlaps(gr_df3, gr_df4, maxgap = distances[m])
#num_overlaps_down <- length(unique(queryHits(overlaps_down)))
overlapping_regions_gr3 <- as.data.frame(gr_df3[queryHits(overlaps_down)])
overlapping_regions_gr4 <- as.data.frame(gr_df4[subjectHits(overlaps_down)])

num_overlaps <- length(overlapping_regions_gr3[,1])+length(overlapping_regions_gr1[,1])

genome_bed <- read.delim(paste0("~/chromosome_lengths.txt"), header = F, sep = "\t", quote = "")
genome_bed <- genome_bed[genome_bed$V1 %in% c(mRNAs_up$chrom,DHS_up$chrom,mRNAs_down$chrom,DHS_down$chrom),]
# Initialize an empty list for chromosome_ranges
chromosome_ranges <- list()

# Loop through each row of the genome data frame
for (k in 1:nrow(genome_bed)) {
  # Add a named element to the chromosome_ranges list
  chromosome_ranges[[genome_bed$V1[k]]] <- c(1, as.numeric(genome_bed$V2[k]))
}

generate_random_coordinates <- function(chrom, start, end, length) {
  chr_range <- chromosome_ranges[[chrom]]
  new_start <- sample(max(chr_range[1], start):(min(chr_range[2], end - length + 1)), 1)
  new_end <- new_start + length - 1
  return(c(new_start, new_end))
}


# Null distribution
null_distribution <- numeric(num_permutations)
num_overlaps_exp <- NULL
# Permutation test
for (j in 1:num_permutations) {
 # j=1
  df <- DHS_up
  df_2 <- df

# Ensure that the length column is numeric
df$length <- as.numeric(df$end - df$start + 1)

# Initialize empty vectors to store new start and end coordinates
new_starts <- numeric(nrow(df))
new_ends <- numeric(nrow(df))

# Use a for loop to generate new coordinates
for (k in seq_along(new_starts)) {
  #k=1
  chrom <- df$chrom[k]  # Use the correct column name for the chromosome
  new_coordinates <- generate_random_coordinates(chrom[k], df$start[k], df$end[k], df$length[k])
  new_starts[k] <- new_coordinates[1]
  new_ends[k] <- new_coordinates[2]
}

# Add new start and end coordinates to the dataframe
df_2$start <- new_starts
df_2$end <- new_ends
  
  gr_DHS_up <- GRanges(
  seqnames = Rle(df_2$chrom),
  ranges = IRanges(start = df_2$start, end = df_2$end)
)
  df <- mRNAs_up
  df_2 <- df

# Ensure that the length column is numeric
df$length <- as.numeric(df$end - df$start + 1)

# Initialize empty vectors to store new start and end coordinates
new_starts <- numeric(nrow(df))
new_ends <- numeric(nrow(df))

# Use a for loop to generate new coordinates
for (k in seq_along(new_starts)) {
  #k=1
  chrom <- df$chrom[k]  # Use the correct column name for the chromosome
  new_coordinates <- generate_random_coordinates(chrom[k], df$start[k], df$end[k], df$length[k])
  new_starts[k] <- new_coordinates[1]
  new_ends[k] <- new_coordinates[2]
}

# Add new start and end coordinates to the dataframe
df_2$start <- new_starts
df_2$end <- new_ends
  
  gr_mRNA_up <- GRanges(
  seqnames = Rle(df_2$chrom),
  ranges = IRanges(start = df_2$start, end = df_2$end)
)
  
overlaps <- findOverlaps(gr_mRNA_up, gr_DHS_up, maxgap = distances[m])

overlapping_regions_up <- as.data.frame(gr_mRNA_up[queryHits(overlaps)])

  df <- DHS_down
  df_2 <- df

# Ensure that the length column is numeric
df$length <- as.numeric(df$end - df$start + 1)

# Initialize empty vectors to store new start and end coordinates
new_starts <- numeric(nrow(df))
new_ends <- numeric(nrow(df))

# Use a for loop to generate new coordinates
for (k in seq_along(new_starts)) {
  #k=1
  chrom <- df$chrom[k]  # Use the correct column name for the chromosome
  new_coordinates <- generate_random_coordinates(chrom[k], df$start[k], df$end[k], df$length[k])
  new_starts[k] <- new_coordinates[1]
  new_ends[k] <- new_coordinates[2]
}

# Add new start and end coordinates to the dataframe
df_2$start <- new_starts
df_2$end <- new_ends
  
  gr_DHS_down <- GRanges(
  seqnames = Rle(df_2$chrom),
  ranges = IRanges(start = df_2$start, end = df_2$end)
)
  df <- mRNAs_down
  df_2 <- df

# Ensure that the length column is numeric
df$length <- as.numeric(df$end - df$start + 1)

# Initialize empty vectors to store new start and end coordinates
new_starts <- numeric(nrow(df))
new_ends <- numeric(nrow(df))

# Use a for loop to generate new coordinates
for (k in seq_along(new_starts)) {
  #k=1
  chrom <- df$chrom[k]  # Use the correct column name for the chromosome
  new_coordinates <- generate_random_coordinates(chrom[k], df$start[k], df$end[k], df$length[k])
  new_starts[k] <- new_coordinates[1]
  new_ends[k] <- new_coordinates[2]
}
# Add new start and end coordinates to the dataframe
df_2$start <- new_starts
df_2$end <- new_ends
  
  gr_mRNA_down <- GRanges(
  seqnames = Rle(df_2$chrom),
  ranges = IRanges(start = df_2$start, end = df_2$end)
)
  
overlaps_down <- findOverlaps(gr_mRNA_down, gr_DHS_down, maxgap = distances[m])

overlapping_regions_down <- as.data.frame(gr_mRNA_down[queryHits(overlaps_down)])

null_distribution[j] <- length(overlapping_regions_down[,1])+ length(overlapping_regions_up[,1])
#num_overlaps_exp <- c(num_overlaps_exp,length(unique(queryHits(overlaps))))
}

#print(hist(null_distribution, breaks = 100))
mean(null_distribution)
p_value <- sum(null_distribution >= num_overlaps) / num_permutations
o_e <- log2(num_overlaps/mean(null_distribution))
f <- (num_overlaps/mean(null_distribution))
t <- c(name,num_overlaps,mean(null_distribution),f,o_e,p_value)
PT_table <- rbind(PT_table,t)
}

PT_table <- as.data.frame(PT_table)
PT_table[,2:6] <- PT_table[,2:6] %>% mutate_all(as.numeric)
colnames(PT_table) <- c("target","obs_EP_pairs","mean_exp_EP_pairs","enrichment","log2_enrichment","P")

write.table(PT_table, file= paste0("~/permutation_table_randomized_elements_",distances[m],"_table.txt"), sep="\t", quote=F, row.names=F,col.names=F)
}

