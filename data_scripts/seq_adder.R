rm(list = ls(all.names = TRUE))
library(ggplot2)
library(dplyr)
library(Biostrings)

######################################
######## Additional sequences ########
######################################
rm(list = ls(all.names = TRUE))
#For ease (and to account for stops containing introns) the following will pull expanded sequence for each individual base in a stop codon
#get a full set of stop coordinates
stops<-read.csv(file = "/Users/alexmg/shiny_apps/stop_start_app/data/basic_reference_stops.txt", sep = '\t', header = TRUE, na.strings=c("","NA")) #transcript ID, chrom, start, end and gene symbol for all stops (stops with introns span 2 lines)

############################################
########## GET REFERENCE GENOME ############
############################################
fasta<- readDNAStringSet('/Users/alexmg/shiny_apps/stop_start_app/input/GRCh38_latest_genomic.fna', format="fasta", skip=0L, seek.first.rec=FALSE, use.names=TRUE)
seqName <- as.data.frame(names(fasta))
colnames(seqName)[1]<-"header"
seqName <- subset(seqName, grepl("NC", header) & !(grepl("mito", header)))
fasta<-fasta[seqName$header]
## assign new seq names  by mapping fasta seq name to data frame names
seqName$chrom<- str_match(names(fasta), "chromosome\\s*(.*?),\\s*GRC")[,2]
names(fasta) <- seqName[match(names(fasta) , seqName$header) , "chrom"]


##########################################################################
############ GET SEQUENCE FOR REFERENCE STOPS (forward strand) ###########
##########################################################################
stops$left<-stops$end-52 #Includes stop Motif as final 3 bases
stops$right<-stops$end+50
stops$left<-subseq(fasta[stops$chrom], start=stops$left, end=stops$end)
stops$right<-subseq(fasta[stops$chrom], start=(stops$end+1), end=stops$right)
colnames(stops)[2]<-"position"

write.table(stops, file='/Users/alexmg/shiny_apps/stop_start_app/data/expanded_stop_seq.txt', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
