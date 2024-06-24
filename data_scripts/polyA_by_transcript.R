rm(list = ls(all.names = TRUE))
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

#####################################
########## GET TRANSCRIPTS ##########
#####################################
#get a full set of havanna/havanna_ensembl protein coding transcripts
ens_data <- read.table('/Users/alexmg/shiny_apps/stop_start_app/data/full_havanna.txt', header = TRUE, sep = '\t')
ens_data$ensembl_gene_id<- str_match(ens_data$V9, "gene_id\\s*(.*?);\\s*gene_version")[,2]
stops<-subset(ens_data, V3 == "stop_codon")
utr<-subset(ens_data, V3 == "three_prime_utr")


##################################################
############ GET TRANSCRIPT DETAILS ##############
##################################################
utr<-utr[c(1,4,5,11:13)]
colnames(utr)[1:3]<-c("chrom","exon_start","exon_end")

full_utr<-utr[c("chrom","ensembl_transcript_id","ensembl_gene_id")]#Get unique gene/transcript combinations
full_utr<-unique(full_utr)
  
#Get the start and end genomic coordinates per transcript 
fusta<-utr %>%
  group_by(ensembl_transcript_id) %>%
  summarise(min = min(exon_start, na.rm=TRUE))
fusto<-utr %>%
  group_by(ensembl_transcript_id) %>%
  summarise(max = max(exon_end, na.rm=TRUE))

full_utr<-merge(full_utr, fusta)
full_utr<-merge(full_utr, fusto)
colnames(full_utr)[c(4,5)]<-c("start","end")
remove(fusta,fusto)


######################################
########### GET POLYA SITES ##########
######################################
setwd("/Users/alexmg/shiny_apps/stop_start_app/input/")
polya<-read.csv(file = "success_LiftOver_polyA_DB_3_chr_position_hg19_hg38.txt", sep = ',', header = TRUE, na.strings=c("","NA")) 
#PolyAs are already mapped to ENS gene ID, so to get them to a transcript level we need only compare with those already seen in the gene
polya$PAS_rank<-NA
polya$Feature<-NA
mapped_primary_polyas<-polya[0,c(1,4,5,7,14:17,21:24)]
#Find sites per transcript
pb=txtProgressBar(min=0,max = nrow(full_utr), initial = 1)
for(i in 1:nrow(full_utr)){
  setTxtProgressBar(pb,i)
  temp<-polya[0,]
  gene_polya<-polya[which(polya$Ensemble.ID==full_utr$ensembl_gene_id[i]),]
  if(nrow(gene_polya != 0)){
    for(j in 1:nrow(gene_polya)){
      if(between(gene_polya$hg38_Position[j], full_utr$start[i], full_utr$end[i])){
        temp<-rbind(temp,gene_polya[j,])
        colnames(temp)
      }
    }
    temp<-unique(temp)
    t2 <- temp %>% mutate(PAS_rank = dense_rank(dplyr::desc(Mean.RPM)))
    if(nrow(t2 > 0)){
      t2<-t2[,c(1,4,5,7,14:17,21:24)]
      t2$Feature<-full_utr$ensembl_transcript_id[i]
      mapped_primary_polyas[(nrow(mapped_primary_polyas)+1),]<-t2[which(t2$PAS_rank == 1),]
    }
  }
}
close(pb) #Close the progress bar
mapped_primary_polyas<-unique(mapped_primary_polyas)


##############################################################
########## FIND PRIMARY POLYAs PRIOR TO RESCUE STOP ##########
##############################################################
stop<- read.table('/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops.txt', header = TRUE, sep = '\t')
stop<-subset(stop, !(is.na(genomic_start)))
ref_stops<-read.csv(file = "/Users/alexmg/shiny_apps/stop_start_app/data/basic_reference_stops.txt", sep = '\t', header = TRUE, na.strings=c("","NA")) #transcript ID, chrom, start, end and gene symbol for all stops (stops with introns span 2 lines)

stop$polyA_pos<-"stop_prior_to_primary_polyA"
pb=txtProgressBar(min=0,max = nrow(stop), initial = 1)
for(i in 1:nrow(stop)){
  setTxtProgressBar(pb,i)
  if(length(mapped_primary_polyas$hg38_Position[which(mapped_primary_polyas$Feature == stop$ensembl_transcript_id[i])])==0){
    stop$polyA_pos[i]<-"Primary polyadenylation site unknown"
  }
  if(length(mapped_primary_polyas$hg38_Position[which(mapped_primary_polyas$Feature == stop$ensembl_transcript_id[i])])!=0){
    if(between(mapped_primary_polyas$hg38_Position[which(mapped_primary_polyas$Feature == stop$ensembl_transcript_id[i])], 
               min(c(ref_stops$start[which(ref_stops$ensembl_transcript_id == stop$ensembl_transcript_id[i])],stop$genomic_start[i])),
               max(c(ref_stops$start[which(ref_stops$ensembl_transcript_id == stop$ensembl_transcript_id[i])],stop$genomic_start[i])))){
      stop$polyA_pos[i]<-mapped_primary_polyas$hg38_Position[which(mapped_primary_polyas$Feature == stop$ensembl_transcript_id[i])]
    }
  }
}
close(pb)
temp<-stop


###################################################################
########## FILL IN DATA FOR TRANSCRIPTS WITH NO ALT STOP ##########
###################################################################
stop2<- read.table('/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops.txt', header = TRUE, sep = '\t')
stop2<-merge(stop2, stop, all.x=TRUE)
stop2$polyA_pos[which(stop2$flags %in% c("No stop in this frame","No stop in any frame"))]<-stop2$flags[which(stop2$flags %in% c("No stop in this frame","No stop in any frame"))]
stop2$polyA_flag<-stop2$polyA_pos
stop2$polyA_flag[which(!(stop2$polyA_flag %in% c("No stop in this frame","No stop in any frame","Primary polyadenylation site unknown","stop_prior_to_primary_polyA")))]<-"Likely Polyadenylation related NSD"

write.table(stop2, file='/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops_polyA.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)


###################################################
######### GET DISTRIBUTION OF EXTENSIONS ########## 
###################################################
stop2<-read.csv(file = "/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops_polyA.txt", sep = '\t', header = TRUE, na.strings=c("","NA")) #transcript ID, position, position in stop, motif for all NTs in each stop codon
stop3<-stop2[which(stop2$frame==0),]
stop4<-stop3[which(stop3$polyA_flag == "stop_prior_to_primary_polyA"),]
p <- ggplot(stop4, aes(x=percent_protein_gain)) + 
  geom_density()+
  theme_classic()
p

stop4$percentile<-ntile(stop4$percent_protein_gain, 10)
stop4 %>%
  group_by(percentile) %>%
  summarise(max = max(percent_protein_gain, na.rm=TRUE))
stop_percentile<-merge(stop3,stop4,all.x=TRUE)

########################################################
########## REFORMAT FOR CONSERVATION ANALYSIS ##########
########################################################
maps<-read.csv(file = "/Users/alexmg/shiny_apps/stop_start_app/data/MAPS_bed_coord_rescue_stops_full.txt", sep = '\t', header = TRUE, na.strings=c("","NA"))[,c(1,4,6,7)] #transcript ID, position, position in stop, motif for all NTs in each stop codon
colnames(maps)<-c("ensembl_transcript_id","motif_allele_position","motif_allele","reference_motif")
m2<-merge(maps,stop_percentile, by=c("ensembl_transcript_id"),all.x=TRUE)
nrow(as.data.frame(unique(maps$ensembl_transcript_id)))
#m2<-m2[,c(1:4,8:10,23:27)]

m2<-subset(m2, !(is.na(chrom)))
m2$percentile[which(m2$polyA_flag %in% c("No stop in this frame","Likely Polyadenylation related NSD","No stop in any frame"))]<-0

#m2<-na.omit(m2) #Drop the stops that don't have a 3'utr
write.table(m2, file='/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops_cons.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)


######################################
########## ADD LOEUF SCORES ##########
######################################
stop3<- read.table('/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops_polyA.txt', header = TRUE, sep = '\t')
loeuf <- read.table('/Users/alexmg/shiny_apps/stop_start_app/input/supplementary_dataset_11_full_constraint_metrics.tsv.gz', header = TRUE, sep = '\t')[c(2,37)]
colnames(loeuf)<-c("ensembl_transcript_id","LOEUF_decile")
stop3<-merge(stop3,loeuf,all.x=TRUE)
write.table(stop3, file='/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops_loeuf.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

######################################################
########## ADD LOEUF and PolyA to MAPS file ##########
######################################################
maps_before <- read.table('/Users/alexmg/shiny_apps/stop_start_app/data/MAPS_bed_coord_rescue_stops_full.txt', header = TRUE, sep = '\t')
stop3 <- read.table('/Users/alexmg/shiny_apps/stop_start_app/data/Closest_stops_loeuf.txt', header = TRUE, sep = '\t')
stop3$flags[which(stop3$polyA_flag == "Likely Polyadenylation related NSD")]<-"Likely Polyadenylation related NSD"
s3<-subset(stop3, frame == "0")
s3<-s3[,c(1,11,24)]
maps_before<-merge(maps_before, s3, all.x=TRUE)
maps_before$frame[which(maps_before$flags %in% c("No stop in any frame","No stop in this frame", "Likely Polyadenylation related NSD"))]<-"NSD"
maps_before$frame[which(maps_before$frame == "IN_FRAME")]<-"RESCUE_STOP_FOUND"
colnames(maps_before)[11]<-"NSD"

################################################################################
###### GET A BED DETAIL FILE OF REFERENCE STOP CODONS ON A PER-BASE BASIS ######
################################################################################
#.bed detail file of reference stop coordinates, where the info field contains:1. The reference_allele, 2. The stop motif_allele, 3. The stop motif (as seen in mRNA), 4. The Strand, 5. The position in the stop motif (e.g. 1st/2nd/3rd),6. The amino acid extension, 7. The predicted outcome of each substitution (DIS=disrupting/PRES=preserving) in the following format [A=PRES|T=DIS|C=DIS|G=REF], 8. Ensembl transcript ID, 9. Frame - IN_FRAME/NO_IN_FRAME_STOP/No stop in any frame
maps_before$attributes<-paste0(maps_before$reference_allele,";",maps_before$motif_allele,";",maps_before$motif,";",maps_before$strand,";",maps_before$position_in_motif,";",maps_before$aa_extension,";",maps_before$Disrupt_Preserve,";",maps_before$ensembl_transcript_id,";",maps_before$NSD,";",maps_before$LOEUF_decile,sep="")
stops_bed<-maps_before[,c(2,3,4,15)]
stops_bed$chrom[stops_bed$chrom == "X"]<-23
stops_bed$chrom[stops_bed$chrom == "Y"]<-24
stops_bed$start<-as.numeric(as.character(stops_bed$start))
stops_bed$chrom<-as.numeric(as.character(stops_bed$chrom))
stops_bed<-stops_bed[order(stops_bed$start),]
stops_bed<-stops_bed[order(stops_bed$chrom),]
stops_bed$chrom[stops_bed$chrom == "23"]<-"X"
stops_bed$chrom[stops_bed$chrom == "24"]<-"Y"
stops_bed$chrom<-paste("chr",stops_bed$chrom, sep="")
write.table(stops_bed, file='/Users/alexmg/shiny_apps/stop_start_app/data/MAPS_bed_detail_rescue_stops_PolyA_LOEUF.bed', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
