rm(list = ls(all.names = TRUE))
library(ggplot2)
library(dplyr)
library(Biostrings)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
###########################################################
########### GET INTRONS PROXIMAL TO STOP CODONS ###########
###########################################################
#Get a set of intron coordinates (within 50bp of stop) for use in identifying deletions that span an intron/exon junction
#This will be used to flag variants that have potentially confounded outcomes. 

#################################################
########## Get Ensembl transcript data ##########
#################################################
d1 <- read.table('/data/full_havanna.txt', header = TRUE, sep = '\t')

#Find all introns
intron_frame<-subset(d1, biotype == "protein_coding" & (V2 == "havana" | V2 == "ensembl_havana") & (V3 == "three_prime_utr" | V3 == "stop_codon" |V3 == "CDS"))
intron_frame<- intron_frame[order(intron_frame$ensembl_transcript_id, intron_frame$V4), ] 
intron_frame$rel_pos_frame <- ave(intron_frame$V4, intron_frame$ensembl_transcript_id, FUN = seq_along)

intron_frame<-intron_frame %>%
  group_by(ensembl_transcript_id) %>%
  mutate(intron_size = V4 - lag(V5, default = V5[1]))
intron_frame$intron_size<-intron_frame$intron_size-1
intron_frame$intron_size[which(intron_frame$rel_pos_frame == 1 | intron_frame$intron_size == 0)]<-NA
intron_frame<-na.omit(intron_frame)

#Calculate the lengths of Introns within genes
intron_frame$intron_end<-intron_frame$V4-1
intron_frame$intron_start<-(intron_frame$V4 - intron_frame$intron_size)


###########################################
########## FIND PROXIMAL INTRONS ##########
###########################################
#Find Introns that are within 50bp of the start or end of a stop codon
stop_coords <- read.table('/data/MAPS_bed_coord_rescue_stops_full.txt', header = TRUE, sep = '\t')[,c(1,3,4)]
intron_frame<-merge(intron_frame,stop_coords, all.x=TRUE)
int_match<-subset(intron_frame,  (between(intron_end, (end - 50),(end+50)) | between(intron_start, (end - 50),(end+50)) | between(end, intron_start,intron_end)))
int_match<-int_match[c("ensembl_transcript_id","V3","intron_start","intron_end")]
colnames(int_match)<-c("ensembl_transcript_id","Intron_type","intron_start","intron_end")
int_match<-unique(int_match)
write.table(int_match, file='/data/Reference_introns.txt', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

