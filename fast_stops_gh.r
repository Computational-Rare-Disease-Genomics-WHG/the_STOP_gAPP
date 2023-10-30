rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(beepr)
library("Biostrings")
library(biomaRt)
options(scipen=1000000)

#Define stop codons
forward_stops<-PDict(c("TAA","TAG","TGA"))

#Pull in the Ensembl data
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="https://www.ensembl.org")

#Get transcript details for protein coding/protein coding LoF biotype proteins
ens_details <- getBM(attributes = c('ensembl_transcript_id','cds_length','transcript_biotype'), filters='transcript_biotype', values= c('protein_coding','protein_coding_LoF'), mart= ensembl)
colnames(ens_details)<-c("ensembl_transcript_id","cds_length_nt","transcript_biotype")

#Get transcript list and coordinates (for above transcripts)
ens_transcripts <- getBM(attributes = c('ensembl_transcript_id','strand','chromosome_name','exon_chrom_start','exon_chrom_end','3_utr_start','3_utr_end'), filters= 'ensembl_transcript_id', values = c(ens_details$ensembl_transcript_id), mart = ensembl)
ens_transcripts<-na.omit(ens_transcripts) #Drop rows with incomplete values
colnames(ens_transcripts)<-c('ensembl_transcript_id','strand','chromosome_name','exon_chrom_start','exon_chrom_end','utr_start_3','utr_end_3')
ens_transcripts<-merge(ens_transcripts, ens_details,all.x=TRUE)
ens_transcripts$rna_start<-""
ens_transcripts$rna_end<-""
ens_transcripts$exon_sum<-(ens_transcripts$utr_end_3-ens_transcripts$utr_start_3)+1

#Get a list of unique Transcripts
unique_transcripts<-as.data.frame(unique(ens_transcripts$ensembl_transcript_id))
colnames(unique_transcripts)<-"id"

#map each exon to position in RNA (Note - takes a long time)
df2<-ens_transcripts[0,]
pb = txtProgressBar(min = 0, max = nrow(unique_transcripts), initial = 0)
for(i in 1:nrow(unique_transcripts)){
  setTxtProgressBar(pb,i)
  temp<-subset(ens_transcripts, ensembl_transcript_id == unique_transcripts$id[i])
  if(unique(temp$strand=="1")){
    temp$rna_start[which(temp$utr_start_3 == min(temp$utr_start_3))]<-"1"
    temp$rna_end[which(temp$utr_start_3 == min(temp$utr_start_3))]<-temp$exon_sum[which(temp$utr_start_3 == min(temp$utr_start_3))]
    temp<-temp[order(temp$exon_chrom_start, decreasing = FALSE), ]
    if(nrow(temp) > 1){  
      for(r in 2:nrow(temp)){
        temp$rna_start[r]<-as.numeric(temp$rna_end[(r-1)])+1
        temp$rna_end[r]<-(as.numeric(temp$rna_start[r]) + as.numeric(temp$exon_sum[r]))-1
      }
    }
  }
  if(unique(temp$strand=="-1")){
    temp$rna_end[which(temp$utr_end_3 == max(temp$utr_end_3))]<-"1"
    temp$rna_start[which(temp$utr_end_3 == max(temp$utr_end_3))]<-temp$exon_sum[which(temp$utr_end_3 == max(temp$utr_end_3))]
    temp<-temp[order(temp$exon_chrom_start, decreasing = TRUE), ]
    if(nrow(temp) > 1){  
      for(r in 2:nrow(temp)){
        temp$rna_end[r]<-as.numeric(temp$rna_start[(r-1)])+1
        temp$rna_start[r]<-(as.numeric(temp$rna_end[r]) + as.numeric(temp$exon_sum[r]))-1
      }
    }
  }
  df2<-rbind(df2,temp)
}
close(pb)
beep(sound = 3)
ens_transcripts<-df2

#Declare functions
get_position<-function(x,y,z){
  for(x in 1:nrow(z)){
    if(between(y$start, as.numeric(z$rna_start[x]),as.numeric(z$rna_end[x]))){
      temp<-(y$start - as.numeric(z$rna_start[x]))+z$utr_start_3[x]
      return(temp)
    }
  }
}

#1. get the sequence
ens_sequences_3 <- getSequence(id = c(unique_transcripts$id),
                               type = "ensembl_transcript_id", 
                               seqType = "3utr", 
                               mart = ensembl, 
                               verbose = FALSE)

#Prepare output frame and frame to populate missing data with NAs 
stops<-data.frame(start=numeric(),end=numeric(),width=numeric(),motif=character(),Feature=character(),strand=character(),aa_extension=character(),nt_extension=character(),genomic_pos_start=character(),stringsAsFactors=FALSE)
missing_data<-data.frame(start=NA,end=NA,width=NA,motif=NA,aa_extension=NA,nt_extension=NA,genomic_pos_start=NA,stringsAsFactors=FALSE)

pb = txtProgressBar(min = 0, max = nrow(unique_transcripts), initial = 0)
for(i in 1:nrow(unique_transcripts)){
  setTxtProgressBar(pb,i)
  temp<-subset(ens_transcripts, ensembl_transcript_id == unique_transcripts$id[i])
  strand<-unique(temp$strand)
  t2<-stops[0,]#prep output frame
  #1. get the sequence
  s = DNAStringSet(ens_sequences_3[ens_sequences_3$ensembl_transcript_id==unique_transcripts$id[i],1])
  #2. check for stops (inc ones that overlap an exon/intron junction) 
    fs<-matchPDict(forward_stops,s[[1]])
    for(j in 1:3){ #For each stop motif
      out<-as.data.frame(fs[[j]]) #Get stops per motif
      if(nrow(out[which(as.numeric(out$end) %% 3==0),])!=0){#If one or more stop is in frame
        out<-out[which(out$start == min(out$start[which(as.numeric(out$end) %% 3==0)])),] #first in frame stop
        out$Feature<-unique_transcripts$id[i]
        out$aa_extension<-as.numeric(out$end)/3
        out$nt_extension<-as.numeric(out$end)
        out$motif<-as.character(forward_stops[[j]])
        if(unique(temp$strand)=="1"){
          for(p in 1:nrow(temp)){
            if(between(out$start, as.numeric(temp$rna_start[p]),as.numeric(temp$rna_end[p]))){
              out$genomic_pos_start<-(out$start - as.numeric(temp$rna_start[p]))+temp$utr_start_3[p]
              out$strand<-"+"
            }
          }
        }
        if(unique(temp$strand)=="-1"){
          for(p in 1:nrow(temp)){
            if(between(out$start, as.numeric(temp$rna_end[p]), as.numeric(temp$rna_start[p]))){
              out$genomic_pos_start<-temp$utr_end_3[p]-(out$start - as.numeric(temp$rna_end[p]))-2
              out$strand<-"-"
            }
          }
        }
        t2<-rbind(t2,out)
      }
      if(nrow(out[which(as.numeric(out$end) %% 3==0),])==0){
        out<-missing_data
        out$Feature<-unique_transcripts$id[i]
        out$strand<-strand
        t2<-rbind(t2,out)
      }
    }
  #3. return the nearest stop, NT and AA extensions (one row of output per transcript)
  if(nrow(na.omit(t2))>=1){
    t2<-na.omit(t2)
    t2<-t2[which(t2$aa_extension == min(t2$aa_extension)),]
  }
  t2<-unique(t2)
  stops<-rbind(stops,t2)
}
close(pb)
beep(sound = 2) #Notify of completion

#merge with additional transcript info
md<-ens_transcripts[,c(1,3,8)]
md<-unique(md)
colnames(md)[1]<-"Feature"
stops<-merge(stops,md,by="Feature",all.x=TRUE)
colnames(stops)<-c("ENS_Transcript_id","start_rna","end_rna","width","extension_AA","extension_NT","Motif_forward","start_genomic","strand","chrom","CDS_length_nt")

#add aditional info]
stops$end_genomic<-stops$start_genomic+2
stops$CDS_length_aa<-stops$CDS_length_nt/3 #Original AA length
stops$prop_AA_ext<-(stops$extension_AA/stops$CDS_length_aa)*100 #Additional AA %
stops$prop_NT_ext<-(stops$extension_NT/stops$CDS_length_nt)*100 #Additional NT %

stops<-stops[,c(1,7,10,2,3,8,12,9,11,13,5,6,14,15)]

ens_names <- getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id','external_gene_name'), filters="ensembl_transcript_id", values= c(stops$ENS_Transcript_id), mart= ensembl)
colnames(ens_names)<-c("ENS_Transcript_id","ENS_Gene_id","Gene_Symbol")
stops<-merge(stops,ens_names,all.x=TRUE)
rm(df2,ens_details, ens_names,ens_sequences_3,forward_stops,fs,md,missing_data,out,pb,s,t2,temp,i,j,p,r,strand)

#Make plot data 
#For each transcript, get the total mRNA length, the number of bases in the 5'UTR/ 3'UTR/ CDS
#Get transcript list and coordinates (for above transcripts)
#Then a second row with total mRNA length, the number of bases in the 5'UTR/ 3'UTR/ CDS
utr_tbl <- ens_transcripts %>% group_by(ensembl_transcript_id) %>% 
  summarise(utr_length_3=sum(exon_sum),.groups = 'drop')
colnames(utr_tbl)[1]<-"ENS_Transcript_id"

#For transcripts with no alternative stop - make the end of the 3'UTR the same ase the end of the MRNA and flag them as non-stop decay
#Flag any that are more than 150bp downstream of the original stop
s1<-stops[,c(1,8,16,5,9,12)]
utr_tbl<-merge(utr_tbl, s1)
utr_tbl$full_length<-utr_tbl$utr_length_3 + utr_tbl$CDS_length_nt
utr_tbl$Original_cds<-utr_tbl$CDS_length_nt / utr_tbl$full_length
utr_tbl$Proportion_extended<-utr_tbl$extension_NT / utr_tbl$full_length
write.table(utr_tbl, file='plot_stops.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

stops$prop_utr_loss<-stops$extension_NT/utr_tbl$utr_length_3[which(utr_tbl$ENS_Transcript_id==stops$ENS_Transcript_id)]*100
write.table(stops, file='All_ENS_rescue_stops.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
