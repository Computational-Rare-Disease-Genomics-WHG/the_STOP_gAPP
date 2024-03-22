rm(list = ls(all.names = TRUE))
library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(beepr)
library("Biostrings")
library(ggplot2)
options(scipen=1000000)
gc()
start.time <- Sys.time()
#Timed at 28.34 Mins

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
####################################################
############### GET TRANSCRIPT DATA ################
####################################################
#get a full set of havanna/havanna_ensembl protein coding transcripts
ens_data <- read.table('/input/Homo_sapiens.GRCh38.110.gtf', header = FALSE, sep = '\t')
ens_data$biotype<- str_match(ens_data$V9, "gene_biotype\\s*(.*?);\\s*transcript_")[,2]
ens_data<-subset(ens_data, biotype %in% c("protein_coding") & (V2 == "havana" | V2 == "ensembl_havana"))
ens_data<-subset(ens_data, V1 %in% c(1:22,"X","Y"))
ens_data$ensembl_transcript_id<- str_match(ens_data$V9, "transcript_id\\s*(.*?);\\s*transcript_version")[,2]
ens_data$gene_symbol<- str_match(ens_data$V9, "gene_name\\s*(.*?);\\s*gene_source")[,2]

#Some transcripts have unusual stop codons e.g. ACT in ENST00000419164, It looks like when this happens, no stop is recorded in the data, and instead the stop is part of the last CDS exon. For the puropses of the App we do not want to retain these transcripts. 
stops<-subset(ens_data, V3 == "stop_codon")
ens_data<-subset(ens_data, ensembl_transcript_id %in% stops$ensembl_transcript_id)
utr<-subset(ens_data, V3 == "three_prime_utr") #Get the 3' UTR exons
#Store a reference set of transcripts for downstream scripts
write.table(ens_data, file='/data/full_havanna.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
#Get a list of transcripts without UTRs, and exclude them from core transcript list
no_utr<-subset(ens_data, !(ensembl_transcript_id %in% utr$ensembl_transcript_id) & V3=="transcript")[,c(1:3,7,11)]
ens_data<-subset(ens_data, !(ensembl_transcript_id %in% no_utr$ensembl_transcript_id))
cds<-subset(ens_data, V3 == "CDS") #Get the CDS exons


####################################################
############## ADD CDS BASE/AA COUNTS ##############
####################################################
#Calculate the total number of AAs in the CDS
cds$exon_sum<-(cds$V5-cds$V4)+1 #Calculate the total number of bases in each CDS exon including first and last 
ens_transcripts<-cds %>% 
  group_by(ensembl_transcript_id) %>% 
  summarise(CDS_length_nt = sum(exon_sum))
ens_transcripts$CDS_length_aa<-ens_transcripts$CDS_length_nt/3 #This excludes the reference stop
remove(cds)


##########################################################################
############ GET REFERENCE GENOME SEQUENCES AND 'REBUILD' RNA ############
##########################################################################
#get NCBI reference genome sequence and limit to nuclear chromosomes
fasta<- readDNAStringSet('/input/GRCh38_latest_genomic.fna', format="fasta", skip=0L, seek.first.rec=FALSE, use.names=TRUE)
seqName <- as.data.frame(names(fasta))
colnames(seqName)[1]<-"header"
seqName <- subset(seqName, grepl("NC", header) & !(grepl("mito", header)))
fasta<-fasta[seqName$header]
seqName$chrom<- str_match(names(fasta), "chromosome\\s*(.*?),\\s*GRC")[,2]## assign new sequence names  by mapping fasta seq name to data frame names
names(fasta) <- seqName[match(names(fasta) , seqName$header) , "chrom"]

#Add 3' exon sequences to data
utr$seq<-""
utr$seq[which(utr$V7=="-")]<-reverseComplement(subseq(fasta[utr$V1[which(utr$V7=="-")]], start=utr$V4[which(utr$V7=="-")], end=utr$V5[which(utr$V7=="-")]))
utr$seq[which(utr$V7=="+")]<-subseq(fasta[utr$V1[which(utr$V7=="+")]], start=utr$V4[which(utr$V7=="+")], end=utr$V5[which(utr$V7=="+")])
remove(fasta,seqName)
gc()

#determine exon order per transcript
colnames(utr)[1:9]<-c("chrom", "source", "type", "dna_start", "dna_end", "phase", "strand", "other", "attributes")
pos<-utr[which(utr$strand=="+"),]
neg<-utr[which(utr$strand=="-"),]
pos<- pos[order(pos$ensembl_transcript_id, pos$dna_end), ] 
neg<- neg[order(neg$ensembl_transcript_id,-neg$dna_end),]
pos$rel_pos_frame <- ave(pos$dna_end, pos$ensembl_transcript_id, FUN = seq_along)
neg$rel_pos_frame <- ave(neg$dna_end, neg$ensembl_transcript_id, FUN = seq_along)
utr<-rbind(pos,neg)
remove(pos,neg,ens_data)
gc()

#Calculate position in RNA of each 3' exon 
utr$rna_start<-""
utr$rna_end<-""
utr$exon_sum<-(utr$dna_end-utr$dna_start)+1
utr$rna_end <- ave(utr$exon_sum, utr$ensembl_transcript_id, FUN=cumsum)
utr$rna_start<-(utr$rna_end - utr$exon_sum)+1

#Combine all exon sequences to create a 'complete' RNA sequence
full_utr<-utr %>%
  group_by(ensembl_transcript_id) %>%
  summarise(full_UTR_seq=paste(seq,collapse=''))

#Calculate the total length of UTR (RNA)
full_utr$utr_length<-nchar(full_utr$full_UTR_seq)


####################################################
############## FIND DOWNSTREAM STOPS ###############
####################################################
#Define stop codons
forward_stops<-PDict(c("TAA","TAG","TGA"))

#Get all stops (in any frame) in the 3'UTR mRNA 
full_utr$TAA<-""
full_utr$TAG<-""
full_utr$TGA<-""
get_stops<-function(x,y){
  paste(as.data.frame(matchPDict(forward_stops,DNAStringSet(x)[[1]])[[y]])[,c(1)],collapse="|")
}
full_utr$TAA<-lapply(full_utr$full_UTR_seq,get_stops, y=1) #Note these each take ~10 minutes to run
full_utr$TAG<-lapply(full_utr$full_UTR_seq, get_stops,y=2)
full_utr$TGA<-lapply(full_utr$full_UTR_seq, get_stops,y=3)
beep(sound = 2) #Notify of completion

#Store a copy of original data
fu2<-full_utr

#Unnest rescue stop positions, and combine all stops into one column with an additional column as motif
taa<-full_utr[,c(1:4)]
tag<-full_utr[,c(1:3,5)]
tga<-full_utr[,c(1:3,6)]

taa<-taa %>% mutate(TAA = strsplit(as.character(TAA), "\\|")) %>% tidyr::unnest(TAA)
taa<-subset(taa, TAA != "|")
taa$motif<-"TAA"
colnames(taa)[4]<-"stop_RNA_start"
taa<-unique(taa)

tag<-tag %>% mutate(TAG = strsplit(as.character(TAG), "\\|")) %>% tidyr::unnest(TAG)
tag<-subset(tag, TAG != "|")
tag$motif<-"TAG"
colnames(tag)[4]<-"stop_RNA_start"
tag<-unique(tag)

tga<-tga %>% mutate(TGA = strsplit(as.character(TGA), "\\|")) %>% tidyr::unnest(TGA)
tga<-subset(tga, TGA != "|")
tga$motif<-"TGA"
colnames(tga)[4]<-"stop_RNA_start"
tga<-unique(tga)

full_utr<-rbind(taa,tag,tga)
remove(taa,tag,tga,forward_stops,get_stops)
gc()

full_utr<-as.data.frame(full_utr)
full_utr<-full_utr[,c(1,3:5)]


##########################################################################
############ GET GENOMIC COORDINATES OF DOWNSTREAM STOPS #################
##########################################################################
utr<-utr[,c(1,4,5,7,11,12,14:17)]
full_utr<-merge(full_utr, utr, by="ensembl_transcript_id", all.x=TRUE)
full_utr<-subset(full_utr, between(as.numeric(stop_RNA_start), rna_start, rna_end)) #exons containing a stop
colnames(full_utr)[3]<-"Stop_start_position"

#Get Genomic coordinates of the first position in each stop
full_utr$genomic_start<-apply(full_utr, 1, function(t){
  if(t["strand"]=="+"){   
    gp<-(as.numeric(t["Stop_start_position"]) - as.numeric(t["rna_start"])) + as.numeric(t["dna_start"])
  }
  if(t["strand"]=="-"){ 
    gp<-(as.numeric(t["rna_end"]) - as.numeric(t["Stop_start_position"])) + as.numeric(t["dna_start"])
  }
  return(gp)
})

#Add back transcripts without stops, record new stop as the end of the transcript and add non-stop decay flag, and populate as many of the remaiing columns as possible
full_utr$flags<-""
fu2<-fu2[,c(1,3)]
fu2<-unique(fu2)
fu2<-subset(fu2, !(ensembl_transcript_id %in% full_utr$ensembl_transcript_id))
fu2$Stop_start_position<-fu2$utr_length
fu2$flags<-"No stop in any frame"
fu2<-merge(fu2,utr,all.x=TRUE)
fu2[setdiff(names(full_utr), names(fu2))] <- NA
full_utr<-rbind(full_utr,fu2)
nrow(as.data.frame(unique(full_utr$ensembl_transcript_id)))
remove(fu2,utr)
gc()


##############################################################
############ CALCULATE ADDITIONAL ANNOTATION DATA ############
##############################################################
#Add CDS length
full_utr<-merge(full_utr,ens_transcripts,all.x=TRUE)

#Calculate extensions (AA, NT, frame, and Prop UTR loss)
full_utr$nt_extension<-(as.numeric(full_utr$Stop_start_position)-1)+3 #+3 to add the reference stop as an amino acid that is otherwise not counted
full_utr$aa_extension<-(as.numeric(full_utr$nt_extension)/3)
full_utr$frame<-as.numeric(full_utr$nt_extension) %% 3
full_utr$frame[which(full_utr$flags == "No stop in any frame")]<-"No Stop in any frame"

#Include non-stop decay entries for transcripts where the only stops are in alternate frames (So there will always be at least one entry per transcript per frame)
no_in<-subset(ens_transcripts[1], !(ensembl_transcript_id %in% full_utr$ensembl_transcript_id[which(full_utr$frame == 0)]))
no_1<-subset(ens_transcripts[1], !(ensembl_transcript_id %in% full_utr$ensembl_transcript_id[which(full_utr$frame == 1)]))
no_2<-subset(ens_transcripts[1], !(ensembl_transcript_id %in% full_utr$ensembl_transcript_id[which(full_utr$frame == 2)]))

##FRAME
#inframe
no_in<-subset(full_utr[,c(1,2,5,8,9,15:20)], ensembl_transcript_id %in% no_in$ensembl_transcript_id)
no_in$nt_extension<-no_in$utr_length+3
no_in$aa_extension<-no_in$nt_extension/3
no_in$frame<-"0"
no_in$flags<-"No stop in this frame"
#plus1
no_1<-subset(full_utr[,c(1,2,5,8,9,15:20)], ensembl_transcript_id %in% no_1$ensembl_transcript_id)
no_1$nt_extension<-no_1$utr_length+3
no_1$aa_extension<-no_1$nt_extension/3
no_1$frame<-"1"
no_1$flags<-"No stop in this frame"
#plus2
no_2<-subset(full_utr[,c(1,2,5,8,9,15:20)], ensembl_transcript_id %in% no_2$ensembl_transcript_id)
no_2$nt_extension<-no_2$utr_length+3
no_2$aa_extension<-no_2$nt_extension/3
no_2$frame<-"2"
no_2$flags<-"No stop in this frame"
#limit to a single row each
no_in<-unique(no_in)
no_1<-unique(no_1)
no_2<-unique(no_2)
no_in[setdiff(names(full_utr), names(no_in))] <- NA
no_1[setdiff(names(full_utr), names(no_1))] <- NA
no_2[setdiff(names(full_utr), names(no_2))] <- NA
full_utr<-rbind(full_utr,no_in,no_1,no_2)
remove(no_in,no_1,no_2)
gc()

#Calculate Utr losses
full_utr$utr_loss_nt<-full_utr$nt_extension #Not inc ref stop, but does include new stop.
full_utr$prop_utr_loss<-(full_utr$utr_loss_nt / full_utr$utr_length)*100
full_utr$utr_loss_nt[which(full_utr$frame == "No Stop in any frame" | full_utr$frame == "No stop in this frame")]<-full_utr$utr_length[which(full_utr$frame == "No Stop in any frame" | full_utr$frame == "No stop in this frame")]
full_utr$prop_utr_loss[which(full_utr$frame == "No Stop in any frame" | full_utr$frame == "No stop in this frame")] <- "100"

#Get relative distance for each stop for each frame
full_utr$proximity<-""
full_utr$frame_id<-paste0(full_utr$ensembl_transcript_id,"_",full_utr$frame,sep="")
pos<-full_utr[which(full_utr$strand=="+"),]
neg<-full_utr[which(full_utr$strand=="-"),]
none<-full_utr[which(is.na(full_utr$strand)),] #no stop in any frame
pos<- pos[order(pos$frame_id, pos$nt_extension), ] 
neg<- neg[order(neg$frame_id,-neg$nt_extension),]
pos$proximity <- ave(pos$nt_extension, pos$frame_id, FUN = seq_along)
neg$proximity <- ave(neg$nt_extension, neg$frame_id, FUN = seq_along)
full_utr<-rbind(pos,neg,none)
remove(pos, neg, none)
gc()

#Calculate AA % gain
full_utr$percent_protein_gain<-(full_utr$nt_extension / as.numeric(full_utr$CDS_length_nt))*100

#ADD MANE FLAG #INCLUDING Clinical
mane<-read.gff(file = "/input/MANE.GRCh38.v1.2.ensembl_genomic.gff", na.strings = c(".", "?"), GFF3 = TRUE)[c(3,9)]
mane<-subset(mane,type=="transcript")
mane$ensembl_transcript_id<- str_match(mane$attributes, "transcript_id=*(.*?)\\.")[,2]
mane$mane<- str_match(mane$attributes, "tag=*(.*?)[\\,\\;]")[,2]
mane<-mane[c(3,4)]
full_utr<-merge(full_utr,mane,all.x=TRUE)
full_utr$mane[which(is.na(full_utr$mane))]<-"Not a MANE transcript" 


####################################################
################# OUTPUT CORE DATA #################
####################################################
fu2<-full_utr[,c(1:5,8,9,11,12,14:23,25,26)]
write.table(fu2, file='/data/All_stops.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
closest_stop<-subset(fu2, proximity == 1)
closest_stop$strand[which(closest_stop$strand == "-" )]<-"Negative"
closest_stop$strand[which(closest_stop$strand == "+" )]<-"Positive"
write.table(closest_stop, file='/data/Closest_stops.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
remove(full_utr, fu2)
gc()


####################################################
############ GENERATE PLOT DATA FOR APP ############
####################################################
full_utr <- read.table('/data/All_stops.txt', header = TRUE, sep = '\t')
closest_stop <- read.table('/data/Closest_stops.txt', header = TRUE, sep = '\t')
utr_tbl<-closest_stop[,c(1,2,6,7,9,12,14,16,21,20,11)]#ENS_Transcript_id	utr_length_3	strand	ENS_Gene_id	end_rna	CDS_length_nt	extension_NT	frame	MANE full_length	Original_cds
colnames(utr_tbl)<-c("ENS_Transcript_id","utr_length_3","strand","Gene_symbol","end_rna","CDS_length_nt","extension_NT","frame","MANE","Proportion_extended","flags")
utr_tbl$full_length<-utr_tbl$utr_length_3 + utr_tbl$CDS_length_nt
utr_tbl$New_cds<-utr_tbl$CDS_length_nt + utr_tbl$extension_NT
utr_tbl$flags[which(utr_tbl$extension_NT >= 150 & utr_tbl$flags != "No stop in this frame")]<-"Extenion more than 150 NT" #Flag any that are more than 150bp downstream of the original stop
write.table(utr_tbl, file='/data/plot_full_utr.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)


####################################################
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
