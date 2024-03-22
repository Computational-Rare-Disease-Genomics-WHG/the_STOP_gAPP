rm(list = ls(all.names = TRUE))
library(ggplot2)
library(dplyr)
library(Biostrings)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
####################################################
########## GET REFERENCE STOP COORDINATES ##########
####################################################
#get a full set of havanna/havanna_ensembl protein coding transcripts in chrs 1:22/X/Y
stops <- read.table('/input/Homo_sapiens.GRCh38.110.gtf', header = FALSE, sep = '\t')
stops<-subset(stops, V1 %in% c(1:22,"X","Y"))
stops$biotype<- str_match(stops$V9, "gene_biotype\\s*(.*?);\\s*transcript_")[,2]
stops<-subset(stops, biotype %in% c("protein_coding") & (V2 == "havana" | V2 == "ensembl_havana"))
stops$ensembl_transcript_id<- str_match(stops$V9, "transcript_id\\s*(.*?);\\s*transcript_version")[,2]
utr<-subset(stops,V3 == "three_prime_utr")
stops<-subset(stops, V3 == "stop_codon")


####################################################
########## GET REFERENCE GENOME SEQUENCE ###########
####################################################
#get genome and limit to nuclear chromosomes
fasta<- readDNAStringSet('/input/GRCh38_latest_genomic.fna', format="fasta", skip=0L, seek.first.rec=FALSE, use.names=TRUE)
seqName <- as.data.frame(names(fasta))
colnames(seqName)[1]<-"header"
seqName <- subset(seqName, grepl("NC", header) & !(grepl("mito", header)))
fasta<-fasta[seqName$header]
## assign new seq names  by mapping fasta seq name to data frame names
seqName$chrom<- str_match(names(fasta), "chromosome\\s*(.*?),\\s*GRC")[,2]
names(fasta) <- seqName[match(names(fasta) , seqName$header) , "chrom"]


####################################################
############# GET REFERENCE STOP MOTIFS ############
####################################################
stops$motif<-""
stops$motif[which(stops$V7=="-")]<-reverseComplement(subseq(fasta[stops$V1[which(stops$V7=="-")]], start=stops$V4[which(stops$V7=="-")], end=stops$V5[which(stops$V7=="-")]))
stops$motif[which(stops$V7=="+")]<-subseq(fasta[stops$V1[which(stops$V7=="+")]], start=stops$V4[which(stops$V7=="+")], end=stops$V5[which(stops$V7=="+")])
stops$reference_sequence<-""
stops$reference_sequence<-subseq(fasta[stops$V1], start=stops$V4, end=stops$V5)
stops$reference_sequence<-as.character(stops$reference_sequence)
remove(fasta,seqName)
gc()

#Get a list of transcripts without UTRs, and exclude them from core transcript list
no_utr<-subset(stops, !(ensembl_transcript_id %in% utr$ensembl_transcript_id))

#############################################################
####### GET STATS FOR FIRST IN-FRAME DOWNSTREAM STOP ########
#############################################################
new_stops<-read.table('/data/Closest_stops.txt', header = TRUE, sep = '\t')[c(20,16,1,7,11)]
nrow(as.data.frame(subset(stops, !(ensembl_transcript_id %in% new_stops$ensembl_transcript_id) & !(ensembl_transcript_id %in% no_utr$ensembl_transcript_id))))
new_stops<-subset(new_stops, frame == "0" |frame == "No Stop in any frame" |(frame == "0" & flags == "No stop in this frame"))
stops<-merge(stops, new_stops, all.x=TRUE)
stops$percent_protein_gain[which(stops$ensembl_transcript_id %in% no_utr$ensembl_transcript_id)]<-"NO_UTR"
stops$frame[which(stops$ensembl_transcript_id %in% no_utr$ensembl_transcript_id)]<-"NO_UTR"
stops$frame[which(stops$frame =="0")]<-"IN_FRAME"
stops$frame[which(stops$flags == "No stop in this frame")]<-"NO_IN_FRAME_STOP"
stops<-stops[,c(1,2,5,6,8,12:16)]
remove(new_stops)
gc()
stops2<-stops
colnames(stops2)[2:4]<-c("chrom","start","end")
write.table(stops2[,c(1:4,10)], file='/data/basic_reference_stops.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

##################################################################
####### HANDLE REFERENCE STOPS THAT ARE SPLIT BY AN INTRON #######
##################################################################
split_stops<-stops[which(nchar(stops$motif)!=3),]
stops<-stops[which(nchar(stops$motif)==3),]

pos<-split_stops[which(split_stops$V7=="+"),]
neg<-split_stops[which(split_stops$V7=="-"),]
pos<- pos[order(pos$ensembl_transcript_id, pos$V4), ] 
neg<- neg[order(neg$ensembl_transcript_id,-neg$V4),]
pos$rel_pos_frame <- ave(pos$V4, pos$ensembl_transcript_id, FUN = seq_along)
neg$rel_pos_frame <- ave(neg$V4, neg$ensembl_transcript_id, FUN = seq_along)
split_stops<-rbind(pos,neg)
remove(pos,neg)
gc()

#Calculate start and end coordinates of each exon in RNA
split_stops$m_start<-""
split_stops$m_end<-""
split_stops$exon_sum<-(split_stops$V5-split_stops$V4)+1
split_stops$m_end <- ave(split_stops$exon_sum, split_stops$ensembl_transcript_id, FUN=cumsum)
split_stops$m_start<-(split_stops$m_end - split_stops$exon_sum)+1

#Stitch together stop motif
split_stops_motif<-split_stops %>%
  group_by(ensembl_transcript_id) %>%
  summarise(full_stop_motif=paste(motif,collapse=''))
split_stops<-merge(split_stops,split_stops_motif,all.x=TRUE)

#get the single base and double base components
split_singles<-subset(split_stops, m_start == split_stops$m_end)
split_doubles<-subset(split_stops, m_start != split_stops$m_end)

#get correct details for single bases
split_singles$start_genomic<-split_singles$V4 - 1 #strand doesn't matter here
split_singles$end_genomic<-split_singles$V5
split_singles$reference_nt<-split_singles$reference_sequence
split_singles$motif_nt<-split_singles$motif
split_singles$motif<-split_singles$full_stop_motif
split_singles$motif_position<-split_singles$m_start #position in the mRNA
split_singles<-split_singles[,c(1,2,16:19,6,5,20,8,9)]

#get correct details for double bases (first process the 'left' base)
split_doubles2<-split_doubles
split_doubles$motif_position<-""
split_doubles$start_genomic<-split_doubles$V4 - 1
split_doubles$end_genomic<-split_doubles$V4
split_doubles$reference_nt<-substr(split_doubles$reference_sequence, 1, 1)
split_doubles$motif_nt[which(split_doubles$V7=="+")]<-substr(split_doubles$motif[which(split_doubles$V7=="+")],1,1) #position in the mRNA
split_doubles$motif_nt[which(split_doubles$V7=="-")]<-substr(split_doubles$motif[which(split_doubles$V7=="-")],2,2) #position in the mRNA
split_doubles$motif<-split_doubles$full_stop_motif
split_doubles$motif_position[which(split_doubles$V7=="+")]<-split_doubles$m_start[which(split_doubles$V7=="+")] #position in the mRNA
split_doubles$motif_position[which(split_doubles$V7=="-")]<-split_doubles$m_end[which(split_doubles$V7=="-")]
split_doubles<-split_doubles[,c(1,2,16:19,6,5,20,8,9)]

#(then process the 'right' base)
split_doubles2$motif_position<-""
split_doubles2$start_genomic<-split_doubles2$V5 - 1
split_doubles2$end_genomic<-split_doubles2$V5
split_doubles2$reference_nt<-substr(split_doubles2$reference_sequence, 2, 2)
split_doubles2$motif_nt[which(split_doubles2$V7=="+")]<-substr(split_doubles2$motif[which(split_doubles2$V7=="+")],2,2) #position in the mRNA
split_doubles2$motif_nt[which(split_doubles2$V7=="-")]<-substr(split_doubles2$motif[which(split_doubles2$V7=="-")],1,1) #position in the mRNA
split_doubles2$motif<-split_doubles2$full_stop_motif
split_doubles2$motif_position[which(split_doubles2$V7=="+")]<-split_doubles2$m_end[which(split_doubles2$V7=="+")] #position in the mRNA
split_doubles2$motif_position[which(split_doubles2$V7=="-")]<-split_doubles2$m_start[which(split_doubles2$V7=="-")]
split_doubles2<-split_doubles2[,c(1,2,16:19,6,5,20,8,9)]

######################################################
##### Split complete stops by position, and label ####
######################################################
#Generate all positions in each stop
stops$motif_position<-""
pos1<-stops
pos2<-stops
pos3<-stops
#define position in stop (coordinates as reported on the forward strand, stop position as recorded in mRNA)
pos1$motif_position<-"1"
pos1$motif_position[which(pos1$V7=="-")]<-3
pos2$motif_position<-"2"
pos3$motif_position<-"3"
pos3$motif_position[which(pos3$V7=="-")]<-1

#rejig coordinates (Base0 half open)
pos1$end_genomic<-as.numeric(pos1$V4)
pos1$start_genomic<-as.numeric(pos1$V4)-1
pos2$start_genomic<-as.numeric(pos2$V4)
pos2$end_genomic<-as.numeric(pos2$V4)+1
pos3$start_genomic<-as.numeric(pos3$V4)+1
pos3$end_genomic<-as.numeric(pos3$V5)
stops<-rbind(pos1,pos2,pos3)
remove(pos1,pos2,pos3)

#Get the stop NT 
stops$motif_nt<-substr(stops$motif, stops$motif_position, stops$motif_position)

#Get the reference NT
stops$reference_nt<-stops$motif_nt
stops$reference_nt[which(stops$V7=="-" & stops$motif_position == 1)]<-substr(stops$reference_sequence[which(stops$V7=="-" & stops$motif_position == 1)], 3, 3)
stops$reference_nt[which(stops$V7=="-" & stops$motif_position == 2)]<-substr(stops$reference_sequence[which(stops$V7=="-" & stops$motif_position == 1)], 2, 2)
stops$reference_nt[which(stops$V7=="-" & stops$motif_position == 3)]<-substr(stops$reference_sequence[which(stops$V7=="-" & stops$motif_position == 3)], 1, 1)
stops<-stops[,c(1,2,13,12,15,14,6,5,11,8,9)]
stops<-rbind(stops,split_singles,split_doubles,split_doubles2)
remove(split_singles,split_doubles,split_doubles2,split_stops,split_stops_motif)

##################################################################################
################## GENERATE A SET OF STOP PRESERVING CONDITIONS ##################
##################################################################################
stops$Disrupt_Preserve<-""
stops$Disrupt_Preserve[which(stops$reference_nt == "A")]<-"[A=REF|T=DIS|C=DIS|G=DIS]"
stops$Disrupt_Preserve[which(stops$reference_nt == "T")]<-"[A=DIS|T=REF|C=DIS|G=DIS]"
stops$Disrupt_Preserve[which(stops$reference_nt == "C")]<-"[A=DIS|T=DIS|C=REF|G=DIS]"
stops$Disrupt_Preserve[which(stops$reference_nt == "G")]<-"[A=DIS|T=DIS|C=DIS|G=REF]"
### TAA
#A to G at position 2 or 3 [T to C on reverse]
stops$Disrupt_Preserve[which(stops$motif == "TAA" & stops$V7 == "+" & (stops$motif_position == "2" | stops$motif_position == "3"))]<-"[A=REF|T=DIS|C=DIS|G=PRES]"
stops$Disrupt_Preserve[which(stops$motif == "TAA" & stops$V7 == "-" & (stops$motif_position == "2" | stops$motif_position == "3"))]<-"[A=DIS|T=REF|C=PRES|G=DIS]"
### TAG
#G to A at position 3 [C to T reverse]
stops$Disrupt_Preserve[which(stops$motif == "TAG" & stops$V7 == "+" & stops$motif_position == "3")]<-"[A=PRES|T=DIS|C=DIS|G=REF]"
stops$Disrupt_Preserve[which(stops$motif == "TAG" & stops$V7 == "-" & stops$motif_position == "3")]<-"[A=DIS|T=PRES|C=REF|G=DIS]"
### TGA
#G to A at position 2 [C to T reverse]
stops$Disrupt_Preserve[which(stops$motif == "TGA" & stops$V7 == "+" & stops$motif_position == "2")]<-"[A=PRES|T=DIS|C=DIS|G=REF]"
stops$Disrupt_Preserve[which(stops$motif == "TGA" & stops$V7 == "-" & stops$motif_position == "2")]<-"[A=DIS|T=PRES|C=REF|G=DIS]"
colnames(stops)<-c("ensembl_transcript_id", "chrom","start","end","reference_allele","motif_allele","motif","strand","position_in_motif","percent_protein_gain","frame","Disrupt_Preserve")


#######################################################
##### ANNOTATE STOPS WITH NO DOWNSTREAM 'RESCUE' ######
#######################################################
no_stop<-subset(stops,frame == "No Stop in any frame")
stops$frame[which(stops$ensembl_transcript_id %in% no_stop$ensembl_transcript_id)]<-"No Stop in any frame"
write.table(stops, file='/data/MAPS_bed_coord_rescue_stops_full.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)


################################################################################
###### GET A BED DETAIL FILE OF REFERENCE STOP CODONS ON A PER-BASE BASIS ######
################################################################################
#.bed detail file of reference stop coordinates, where the info field contains:1. The reference_allele, 2. The stop motif_allele, 3. The stop motif (as seen in mRNA), 4. The Strand, 5. The position in the stop motif (e.g. 1st/2nd/3rd),6. The percent_protein_gain, 7. The predicted outcome of each substitution (DIS=disrupting/PRES=preserving) in the following format [A=PRES|T=DIS|C=DIS|G=REF], 8. Ensembl transcript ID, 9. Frame - IN_FRAME/NO_IN_FRAME_STOP/No stop in any frame
stops$attributes<-paste0(stops$reference_allele,";",stops$motif_allele,";",stops$motif,";",stops$strand,";",stops$position_in_motif,";",stops$percent_protein_gain,";",stops$Disrupt_Preserve,";",stops$ensembl_transcript_id,";",stops$frame,sep="")
stops_bed<-stops[,c(2,3,4,13)]
stops_bed$chrom[stops_bed$chrom == "X"]<-23
stops_bed$chrom[stops_bed$chrom == "Y"]<-24
stops_bed<-stops_bed[order(as.numeric(as.character(stops_bed$start))),]
stops_bed<-stops_bed[order(as.numeric(as.character(stops_bed$chrom))),]
stops_bed$chrom[stops_bed$chrom == "23"]<-"X"
stops_bed$chrom[stops_bed$chrom == "24"]<-"Y"
stops_bed$chrom<-paste("chr",stops_bed$chrom, sep="")
write.table(stops_bed, file='/data/MAPS_bed_detail_rescue_stops.bed', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

#################################################################
###### GET A LIST OF PRESERVING SINGLE NUCLEOTIDE VARIANTS ######
#################################################################
#Create a quick ref file for stop preserving variants for the stop app
stops$Pres<-"NONE"
### TAA
#A to G at position 2 or 3 [T to C on reverse]
stops$Pres[which(stops$motif == "TAA" & stops$strand == "+" & (stops$position_in_motif == "2" | stops$position_in_motif == "3"))]<-"G"
stops$Pres[which(stops$motif == "TAA" & stops$strand == "-" & (stops$position_in_motif == "2" | stops$position_in_motif == "3"))]<-"C"
### TAG
#G to A at position 3 [C to T reverse]
stops$Pres[which(stops$motif == "TAG" & stops$strand== "+" & stops$position_in_motif == "3")]<-"A"
stops$Pres[which(stops$motif == "TAG" & stops$strand == "-" & stops$position_in_motif == "3")]<-"T"
### TGA
#G to A at position 2 [C to T reverse]
stops$Pres[which(stops$motif == "TGA" & stops$strand == "+" & stops$position_in_motif == "2")]<-"A"
stops$Pres[which(stops$motif == "TGA" & stops$strand == "-" & stops$position_in_motif == "2")]<-"T"
stops<-stops[which(stops$Pres != "NONE"),c(1,2,4,5,14)]

write.table(stops, file='/data/stop_preserving_variants.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

