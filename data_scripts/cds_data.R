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

####################################################
############### GET TRANSCRIPT DATA ################
####################################################
#get a full set of havanna/havanna_ensembl protein coding transcripts
ens_data <- read.table('/Users/alexmg/shiny_apps/stop_start_app/input/Homo_sapiens.GRCh38.110.gtf', header = FALSE, sep = '\t')
ens_data$biotype<- str_match(ens_data$V9, "gene_biotype\\s*(.*?);\\s*transcript_")[,2]
ens_data<-subset(ens_data, biotype %in% c("protein_coding") & (V2 == "havana" | V2 == "ensembl_havana"))
ens_data<-subset(ens_data, V1 %in% c(1:22,"X","Y"))
ens_data$ensembl_transcript_id<- str_match(ens_data$V9, "transcript_id\\s*(.*?);\\s*transcript_version")[,2]
ens_data$gene_symbol<- str_match(ens_data$V9, "gene_name\\s*(.*?);\\s*gene_source")[,2]

#Some transcripts have unusual stop codons e.g. ACT in ENST00000419164, It looks like when this happens, no stop is recorded in the data, and instead the stop is part of the last CDS exon. For the puropses of the App we do not want to retain these transcripts. 
stops<-subset(ens_data, V3 == "stop_codon")
ens_data<-subset(ens_data, ensembl_transcript_id %in% stops$ensembl_transcript_id)
utr<-subset(ens_data, V3 == "three_prime_utr") #Get the CDS exons

#Get a list of transcripts without UTRs, and exclude them from core transcript list
no_utr<-subset(ens_data, !(ensembl_transcript_id %in% utr$ensembl_transcript_id) & V3=="transcript")[,c(1:3,7,11)]
ens_data<-subset(ens_data, !(ensembl_transcript_id %in% no_utr$ensembl_transcript_id))
cds<-subset(ens_data, V3 == "CDS") [c(1,4,5,7,11,12)]#Get the CDS exons

##########################################################################
############ GET REFERENCE GENOME SEQUENCES AND 'REBUILD' RNA ############
##########################################################################
#get NCBI reference genome sequence and limit to nuclear chromosomes
fasta<- readDNAStringSet('/Users/alexmg/shiny_apps/stop_start_app/input/GRCh38_latest_genomic.fna', format="fasta", skip=0L, seek.first.rec=FALSE, use.names=TRUE)
seqName <- as.data.frame(names(fasta))
colnames(seqName)[1]<-"header"
seqName <- subset(seqName, grepl("NC", header) & !(grepl("mito", header)))
fasta<-fasta[seqName$header]
seqName$chrom<- str_match(names(fasta), "chromosome\\s*(.*?),\\s*GRC")[,2]## assign new sequence names  by mapping fasta seq name to data frame names
names(fasta) <- seqName[match(names(fasta) , seqName$header) , "chrom"]

#Add CDS exon sequences to data
cds$seq<-""
cds$seq[which(cds$V7=="-")]<-reverseComplement(subseq(fasta[cds$V1[which(cds$V7=="-")]], start=cds$V4[which(cds$V7=="-")], end=cds$V5[which(cds$V7=="-")]))
cds$seq[which(cds$V7=="+")]<-subseq(fasta[cds$V1[which(cds$V7=="+")]], start=cds$V4[which(cds$V7=="+")], end=cds$V5[which(cds$V7=="+")])
remove(fasta,seqName,no_utr,utr)
gc()

#determine exon order per transcript
colnames(cds)<-c("chrom","dna_start", "dna_end","strand", "ensembl_transcript_id", "Gene", "Sequence")
pos<-cds[which(cds$strand=="+"),]
neg<-cds[which(cds$strand=="-"),]
pos<- pos[order(pos$ensembl_transcript_id, pos$dna_end), ] 
neg<- neg[order(neg$ensembl_transcript_id,-neg$dna_end),]
pos$rel_pos_frame <- ave(pos$dna_end, pos$ensembl_transcript_id, FUN = seq_along)
neg$rel_pos_frame <- ave(neg$dna_end, neg$ensembl_transcript_id, FUN = seq_along)
cds<-rbind(pos,neg)
remove(pos,neg,ens_data)
gc()

#Calculate position in RNA of each 3' exon 
cds$rna_start<-""
cds$rna_end<-""
cds$exon_sum<-(cds$dna_end-cds$dna_start)+1
cds$rna_end <- ave(cds$exon_sum, cds$ensembl_transcript_id, FUN=cumsum)
cds$rna_start<-(cds$rna_end - cds$exon_sum)+1
cds<-cds[c(1:10)]

#Combine all exon sequences to create a 'complete' RNA sequence
full_cds<-cds %>%
  group_by(ensembl_transcript_id) %>%
  summarise(full_cds_seq=paste(Sequence,collapse=''))

#Calculate the total length of cds (RNA)
full_cds$cds_length<-nchar(full_cds$full_cds_seq)
full_cds<-full_cds[c(1,3)]

write.table(cds, file='/Users/alexmg/shiny_apps/stop_start_app/data/cds_exon_sequences.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
write.table(full_cds, file='/Users/alexmg/shiny_apps/stop_start_app/data/cds_exon_length.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
