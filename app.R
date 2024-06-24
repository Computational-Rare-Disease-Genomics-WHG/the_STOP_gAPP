rm(list = ls(all.names = TRUE))
library(shiny)
library(shinythemes)
library(ggplot2); theme_set(theme_linedraw())
library(utils)
library(datasets)
library(shinyjs)
library(png)
library(dplyr)
library(stringr)
library(Biostrings)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##### Test Variants #####
#x<-"cHr5:88722604_T/G" SNV (Disrupting) MEF2C -REV check transcript ENST00000504921
#x<-"cHr12:8941818_C/T" SNV (Preserving) M6PR -REV t<-"ENST00000000412"
#x<-"cHr5:88722604_TC/GA" MNV MEF2C -REV t<-"ENST00000504921"
#x<-"cHr5:88722604_T/TA" Insertion MEF2C -REV
#x<-"cHr5:88722604_TCA/T" Deletion MEF2C -REV #Double check transcript ENST00000503544 
#x<-"cHr12:106992782_TTTACCTAA/T" Deletion overlapping intron CRY1 -REV t<-"ENST00000549356" <-Does not loadfor T8527(Mane Select)
#x<-"chr8:24508544_AGAG/TTGC" MNV starting in an intron (in a split stop)
#x<-"chr13:97347004_GCGCAGCACCAAGCCAACCAAGCTGCGGTGGCCGCCCAGGCAGCCGCGGCCGCGGCCACAGTCA/GAGC" t<-"ENST00000469707" MNV overlapping an upstream intron
#x<-"chR19:41331050_CCTCA/CGTTA" t<-"ENST00000221930" Stop preserving MNV reverse TGFB1 #Not detecting the replacement stop at the wt position. check the reverse counting is adjusting (adding 2 to location) and figure out why it cant see the first one. 
#x<-"chR16:67637870_TGATGG/TAAGCCATT" t<-"ENST00000264010" Stop preserving MNV forward CTCF
#x<-"chR16:67637870_CGGTGA/TCCGGGTAA" t<-"ENST00000264010" Stop Creating MNV forward CTCF +1AA
#x<-"chR16:67637867_CGGTGA/TAAGCC" t<-"ENST00000264010" Stop Creating MNV forward CTCF -1AA

############################################################################################################################################

### COLLECT DATA ###
stops<-read.csv(file = "data/Closest_stops_loeuf.txt", sep = '\t', header = TRUE, na.strings=c("","NA")) #The closest alternative stop in each of the three frames
utr_tbl<-read.csv(file = "data/plot_full_utr.txt", sep = '\t', header = TRUE, na.strings=c("","NA")) #pre-formatted data specifically for plotting
pres_var<-read.csv(file = "data/stop_preserving_variants.txt", sep = '\t', header = TRUE, na.strings=c("","NA")) #All SNVs that 'preserve' a stop sequence
pres_var$match<-paste0(pres_var$chrom,":",pres_var$end,"_",pres_var$reference_allele,"/",pres_var$Pres,sep="")
ref_stops<-read.csv(file = "data/basic_reference_stops.txt", sep = '\t', header = TRUE, na.strings=c("","NA")) #transcript ID, chrom, start, end and gene symbol for all stops (stops with introns span 2 lines)
introns<-read.csv(file = "data/Reference_introns.txt", sep = '\t', header = FALSE, na.strings=c("","NA")) #Coordinates of introns within 50bp of a stop
stop_seq<-read.csv(file = "data/expanded_stop_seq.txt", sep = '\t', header = FALSE, na.strings=c("","NA")) #Stop codons with details of DNA sequence either side (by base)

#Add Mane Flag to transcript dropdown
mane<-unique(stops[which(stops$mane %in% c("MANE_Select","MANE_Plus_Clinical")),c(1,21)])
ref_stops<-merge(ref_stops,mane,all.x=TRUE)
ref_stops$mane[which(ref_stops$mane == "MANE_Select")]<- "MANE (Select)"
ref_stops$mane[which(ref_stops$mane == "MANE_Plus_Clinical")]<- "MANE (+Clinical)"
ref_stops$mane_transcript<-ref_stops$ensembl_transcript_id
ref_stops$mane_transcript[which(ref_stops$mane %in% c("MANE (Select)","MANE (+Clinical)"))]<-paste0(ref_stops$ensembl_transcript_id[which(ref_stops$mane %in% c("MANE (Select)","MANE (+Clinical)"))],"<-",ref_stops$mane[which(ref_stops$mane %in% c("MANE (Select)","MANE (+Clinical)"))])
ref_stops<-subset(ref_stops, ensembl_transcript_id %in% utr_tbl$ENS_Transcript_id) #drop transcripts without a 3'UTR

stops$polyA_flag[which(stops$polyA_flag == "stop_prior_to_primary_polyA")]<-"New stop before primary polyA site"
############################################################################################################################################

### DEFINE FUNCTIONS ###

#Define variant annotation functions
var_split<-function(x){ #Separates input variant into constituent parts
  x<-gsub("chr", "", x, ignore.case = T)
  strarray<-as.data.frame(strsplit(x, split = "[:|,\\_|,\\/]")) #array with entries - chrom;pos;ref;alt
  return(strarray)
}

tran_strip<-function(t){ #Format Transcript to app form
  t2<-as.data.frame(strsplit(t, split = "[\\<]"))[1,1]
  return(t2)
}

#Determine if the variant is an SNV/balanced-MNV/Insertion/Deletion
get_var_type<-function(x){
  strarray<-var_split(x)
  var_type_f<-""
  if(nchar(strarray[3,1])==1 & nchar(strarray[4,1])==1){ #If ref and alt are a single NT
    var_type_f<-"SNV"
  }
  if(nchar(strarray[3,1])>1 & nchar(strarray[4,1])==nchar(strarray[3,1])){ #If it is a balanced variant greater than 1NT 
    var_type_f<-"MNV"
  }
  if((nchar(strarray[4,1])) > (nchar(strarray[3,1]))){ #If alt is larger than ref
    var_type_f<-"Insertion"
  }
  if((nchar(strarray[4,1])) < (nchar(strarray[3,1]))){ #If alt is shorter than ref
    var_type_f<-"Deletion"
  }
  return(var_type_f)
}


# match stop codons 
get_forward_stops<-function(x,y){
  forward_stops<-PDict(c("TAA","TAG","TGA"))#Define stop codons
  m<-paste(as.data.frame(matchPDict(forward_stops,DNAStringSet(x)[[1]])[[y]])[,c(1)],collapse="|")
  return(m)
}
get_reverse_stops<-function(x,y){
  reverse_stops<-PDict(c("TTA","CTA","TCA"))#Define reversed and complemented stop codons
  m<-paste(as.data.frame(matchPDict(reverse_stops,DNAStringSet(x)[[1]])[[y]])[,c(1)],collapse="|")
  return(m)
}

#search a specific string for stops
get_stops_in_string<-function(x,s){
  stop_frame<-data.frame(TAA=character(), TAG=character(), TGA=character())#Get all stops (in any frame)
  stop_frame[1,]<-c("","","")
  if(s=="Positive"){  
    stop_frame$TAA<-lapply(x, get_forward_stops, y=1) 
    stop_frame$TAG<-lapply(x, get_forward_stops, y=2)
    stop_frame$TGA<-lapply(x, get_forward_stops, y=3)
    return(stop_frame)
  }
  if(s=="Negative"){  
    stop_frame$TAA<-lapply(x, get_reverse_stops, y=1) 
    stop_frame$TAG<-lapply(x, get_reverse_stops, y=2)
    stop_frame$TGA<-lapply(x, get_reverse_stops, y=3)
    return(stop_frame)
  }
}

#Check if MNV preserves or creates a stop
get_mnv_impact<-function(x,t){
  if(nrow(stop_seq[which(stop_seq$V1 == t),]) == 1){ #This will only look at stops that do not contain an intron 
    #Get 50bp either side of the stop, and the stop motif & Find out the variant position relative to the stop
    seq_long<-paste0(stop_seq$V6[which(stop_seq$V1 == t)],stop_seq$V7[which(stop_seq$V1 == t)])
    sp = as.numeric(unique(ref_stops$start[which(ref_stops$mane_transcript == t | ref_stops$mane_transcript == paste0(t,"<-MANE (Select)"))])) #Start of the stop codon
    vp<-as.numeric(strsplit(x, split = "[:|,\\_|,\\/]")[[1]][2]) #Start of the variant
    #Replace ref with alt - long winded due to issues replacing string with longer string
    leader_seq<-substr(seq_long, (1), ((sp-vp)+(50))) #All bases up to the variant
    flexi_string<-substr(seq_long, ((sp-vp)+51), ((sp-vp)+(51+(nchar(strsplit(x, split = "[:|,\\_|,\\/]")[[1]][3])-1)))) #Reference sequence
    flexi_string<- strsplit(x, split = "[:|,\\_|,\\/]")[[1]][4]
    trailing_seq<-substr(seq_long, ((sp-vp)+(51+(nchar(strsplit(x, split = "[:|,\\_|,\\/]")[[1]][3])-1)))+1, (nchar(seq_long))) #All bases after the variant
    seq_new<-paste0(leader_seq,flexi_string,trailing_seq,sep="")
    #Find new stops
    d1<-get_stops_in_string(seq_new,unique(stops$strand[which(stops$ensembl_transcript_id == t)])) #Get stop positions
    
    #Unnest stop positions, and combine all stops into one column with an additional column as motif
    taa<-as.data.frame(d1[,1])
    colnames(taa)<-"TAA"
    tag<-as.data.frame(d1[,2])
    colnames(tag)<-"TAG"
    tga<-as.data.frame(d1[,3])
    colnames(tga)<-"TGA"
    
    taa<-taa %>% mutate(TAA = strsplit(as.character(TAA), "\\|")) %>% tidyr::unnest(TAA)
    taa<-subset(taa, TAA != "|")
    taa$motif<-"TAA"
    colnames(taa)[1]<-"stop_pos"
    taa<-unique(taa)
    
    tag<-tag %>% mutate(TAG = strsplit(as.character(TAG), "\\|")) %>% tidyr::unnest(TAG)
    tag<-subset(tag, TAG != "|")
    tag$motif<-"TAG"
    colnames(tag)[1]<-"stop_pos"
    tag<-unique(tag)
    
    tga<-tga %>% mutate(TGA = strsplit(as.character(TGA), "\\|")) %>% tidyr::unnest(TGA)
    tga<-subset(tga, TGA != "|")
    tga$motif<-"TGA"
    colnames(tga)[1]<-"stop_pos"
    tga<-unique(tga)
    d1<-rbind(taa,tag,tga)
    d1$frame<-""
    
    #If the first position in the original sequence was base 2 in the codon (where the positions are 1,2,3) calculate the new positions
    if(unique(stops$strand[which(stops$ensembl_transcript_id == t)])=="Positive"){  
      d1$frame<-as.numeric(as.numeric(d1$stop_pos)-51) %% 3
      d1<-subset(d1, frame == 0)
      if(nrow(d1)>0){
        d1<-subset(d1, stop_pos == min(d1$stop_pos)) #Get the earliest in-frame stop
        d1$aa_difference<-""
        d1$aa_difference<-((as.numeric(d1$stop_pos)-51)/3)-1
        d1<-d1$aa_difference
      }
      if(is.data.frame(d1)){
        if(nrow(d1)==0){
          d1<-"Stop disrupting "
        }
      }
      return(d1)
    }
    if(unique(stops$strand[which(stops$ensembl_transcript_id == t)])=="Negative"){ 
      d1$stop_pos<-as.numeric(d1$stop_pos)+2 #Add 2 to each to get the start rather than the end, and count from the right
      d1$frame<-as.numeric(53-as.numeric(d1$stop_pos)) %% 3
      d1<-subset(d1, frame==0)
      if(nrow(d1)>0){
        d1<-subset(d1, stop_pos == max(d1$stop_pos)) #Get the earliest in-frame stop
        d1$aa_difference<-""
        d1$aa_difference<-(53-as.numeric(d1$stop_pos))/3
        d1<-d1$aa_difference
      }
      if(is.data.frame(d1)){
        if(nrow(d1)==0){
          d1<-"Stop disrupting "
        }
      }
      return(d1)
      }
    }
  if(nrow(stop_seq[which(stop_seq$V1 == t),]) != 1){
    d1 <- "intron"
    return(d1)
  }
}
  
#Check if SNV preserves a stop (precomputed)
get_pres<-function(x){ #Checks variant against a precomputed list of stop creating SNVs
  x<-gsub("chr", "", x, ignore.case = T)
  if(x %in% pres_var$match){
    dp<-"Stop preserving "
  }
  if(!x %in% pres_var$match){
    dp<-"Stop disrupting "
  }
  return(dp)
}

#Check if variant intersects or contains an intron - these can have more complex interpretations including splice disruption
get_intron_vars<-function(x,t){ 
  strarray<-var_split(x)
  t<-tran_strip(t)
  intR<-introns[which(introns$V1 == t & ( (between(introns$V4, as.numeric(strarray[2,1]), (as.numeric(strarray[2,1]))+(nchar(strarray[3,1])-1))) | (between(introns$V3, as.numeric(strarray[2,1]), (as.numeric(strarray[2,1])+(nchar(strarray[3,1])-1)))) | ((as.numeric(strarray[2,1]) <= introns$V3) & (as.numeric(strarray[2,1])+(nchar(strarray[3,1])-1)) >= introns$V4))),] #If intron starts inside the variant or If intron ends within variant or if intron is wholly contained within the variant
  if(nrow(intR)>0){
    return("intron")
  }
  if(nrow(intR)==0){
    return("not_intron")
  }
}

#Identify frame shifts
get_frame<-function(x){ 
  strarray<-var_split(x)
  frame_f<-""
  val_1<-""
  type<-""
  if(nchar(strarray[3,1])==1 & nchar(strarray[4,1])==1){ #SNV #If ref and alt are a single NT
    frame_f<-0
  }
  if(nchar(strarray[3,1])>1 & nchar(strarray[4,1])==nchar(strarray[3,1])){ #MNV #If it is a balanced variant greater than 1NT
    frame_f<-0
  }
  if((nchar(strarray[4,1])) > (nchar(strarray[3,1]))){ #Insertion #If alt is larger than ref
    frame_f<-(nchar(strarray[4,1]) - (nchar(strarray[3,1]))) %% 3 #as we are just looking for the difference, we don't need to deduct 1 in order to account for the "anchoring" base at the start
  }
  if((nchar(strarray[4,1])) < (nchar(strarray[3,1]))){ #Deletion #If alt is shorter than ref
    frame_f<-(nchar(strarray[3,1]) - (nchar(strarray[4,1]))) %% 3 #as we are just looking for the difference, we don't need to deduct 1 in order to account for the "anchoring" base at the start
    val_1<-(nchar(strarray[3,1]) - (nchar(strarray[4,1])))
    type<-"D"
  }
  return(frame_f)
}

#Get a list of affected transcripts
get_transcripts<-function(x){ #Retrieves a list of transcript IDs where the stop codon intersects with the variant
  strarray<-var_split(x)
  if(nchar(strarray[3,1])==1 & nchar(strarray[4,1])==1){   #SNV #If ref and alt are a single NT
    choices_f = unique(ref_stops$mane_transcript[which(ref_stops$chrom == strarray[1,1] & (ref_stops$start <= as.numeric(strarray[2,1])) & (ref_stops$end >= as.numeric(strarray[2,1])))]) #If the SNV position is within an identified reference stop
  }
  if(nchar(strarray[3,1])>1 & nchar(strarray[4,1])==nchar(strarray[3,1])){ #MNV #If it is a balanced variant greater than 1NT
    choices_f = unique(ref_stops$mane_transcript[which(ref_stops$chrom == strarray[1,1] & ((between(ref_stops$end, as.numeric(strarray[2,1]), (as.numeric(strarray[2,1]))+(nchar(strarray[3,1])-1))) | (between(ref_stops$start, as.numeric(strarray[2,1]), (as.numeric(strarray[2,1])+(nchar(strarray[3,1])-1)))) | ((as.numeric(strarray[2,1]) <= ref_stops$start) & (as.numeric(strarray[2,1])+(nchar(strarray[3,1])-1)) >= ref_stops$end)))]) #If the MNV starts or ends within an identified reference stop, or surrounds it
  }
  if((nchar(strarray[4,1])) > (nchar(strarray[3,1]))){ #Insertion #If alt is larger than ref
    choices_f = unique(ref_stops$mane_transcript[which(ref_stops$chrom == strarray[1,1] & ((between(ref_stops$end, as.numeric(strarray[2,1]), (as.numeric(strarray[2,1]))+(nchar(strarray[3,1])-1))) | (between(ref_stops$start, as.numeric(strarray[2,1]), (as.numeric(strarray[2,1])+(nchar(strarray[3,1])-1)))) | ((as.numeric(strarray[2,1]) <= ref_stops$start) & (as.numeric(strarray[2,1])+(nchar(strarray[3,1])-1)) >= ref_stops$end)))]) #Follows the same logic as MNV
  }
  if((nchar(strarray[4,1])) < (nchar(strarray[3,1]))){ #Deletion #If alt is shorter than ref
    choices_f = unique(ref_stops$mane_transcript[which(ref_stops$chrom == strarray[1,1] & ((between(ref_stops$end, (as.numeric(strarray[2,1])+1), ((as.numeric(strarray[2,1])+1))+(nchar(strarray[3,1])-2))) | (between(ref_stops$start, (as.numeric(strarray[2,1])+1), ((as.numeric(strarray[2,1])+1)+(nchar(strarray[3,1])-2)))) | (((as.numeric(strarray[2,1])+1) <= ref_stops$start) & (((as.numeric(strarray[2,1])+1)+(nchar(strarray[3,1])-2)) >= ref_stops$end))))]) 
  }
  return(choices_f)
}

## Generate the plot for GUI ##
get_plot<-function(x,y,t,v){ #Transcript, frame,  type, variant
  #v<-"cHr5:88722604_TCA/T"
  #t<-"Deletion"
  #z<-""
  #y<-0
  #x<-tran_strip("ENST00000503554")
  if(t == "SNV"){
    z<-get_pres(v) #Determine if the variant is stop disrupting
  }
  if(t != "SNV"){
    z<-get_mnv_impact(v,x)
  }
  z2<-z
  if(grepl("[[:digit:]]", z)){
    z2<-"Stop creating "
  }
  if(nchar(x) > 15){
    lab<-paste0("Transcript ",x," has an unexpectedly long ID!")
    p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
  }
  if(nchar(x) == 15 & !(x %in% utr_tbl$ENS_Transcript_id)){
    lab<-paste0("Transcript ",x," Is not a known transcript.")
    p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
  }
  if(x %in% utr_tbl$ENS_Transcript_id){
    plotframe<-data.frame(prop_pos=numeric(0),rowheight=numeric(0),weight=numeric(0),category=character(0),region=character(0),stringsAsFactors = FALSE)
    plotframe[nrow(plotframe) + 1,] = c(0,0.7,9,"Wild Type","CDS")
    plotframe[nrow(plotframe) + 1,] = c(utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id== x & utr_tbl$frame == y)],0.7,9,"Wild Type","CDS")
    plotframe[nrow(plotframe) + 1,] = c(utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)],0.7,4,"Wild Type","UTR")
    plotframe[nrow(plotframe) + 1,] = c((utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]+utr_tbl$utr_length_3[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]),0.7,4,"Wild Type","UTR")
    plotframe[nrow(plotframe) + 1,] = c(0,0.4,9,"Stop lost","CDS")
    plotframe[nrow(plotframe) + 1,] = c(utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y )],0.4,9,"Stop lost","CDS")
    plotframe[nrow(plotframe) + 1,] = c(utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y )],0.4,9,"Stop lost","Extension")
    plotframe[nrow(plotframe) + 1,] = c((utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]+utr_tbl$extension_NT[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]),0.4,9,"Stop lost","Extension")
    plotframe[nrow(plotframe) + 1,] = c((utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]+utr_tbl$extension_NT[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]),0.4,4,"Stop lost","UTR")
    plotframe[nrow(plotframe) + 1,] = c((utr_tbl$CDS_length_nt[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]+utr_tbl$utr_length_3[which(utr_tbl$ENS_Transcript_id == x & utr_tbl$frame == y)]),0.4,4,"Stop lost","UTR")
    frametext<-""
    if(y==0){ frametext<-" - In Frame" }
    if(y!=0){ frametext<-paste0(" - Out of Frame + ",y) }
    stops$LOEUF_decile[which(is.na(stops$LOEUF_decile))]<-"Not recorded"
    stops$polyA_flag[which(stops$polyA_flag == "New stop before primary polyA site" & stops$frame == y)]<-paste0(" - ","New stop before primary polyA site",sep="")
    stops$polyA_flag[which(stops$polyA_flag == "No stop in any frame" & stops$frame == y)]<-paste0(" - ","No stop in any frame",sep="")
    stops$polyA_flag[which(stops$polyA_flag == "No stop in this frame" & stops$frame == y)]<-paste0(" - ","No stop in this frame",sep="")
    stops$polyA_flag[which(stops$polyA_flag == "Primary polyadenylation site unknown" & stops$frame == y)]<-paste0(" - ","Primary polyadenylation site unknown",sep="")
    group.colors <- c(CDS = "darkcyan", Extension = "darkgoldenrod2", UTR ="black", POLYA ="white")
    if(!is.na(stops$rna_start[which(stops$ensembl_transcript_id==x & stops$frame==y)])){
      p<- ggplot(plotframe,aes(as.numeric(prop_pos),as.numeric(rowheight),linewidth=as.numeric(weight),group=category,colour=region))+
        geom_line(size=8)+
        scale_y_continuous(limits = c(0, 1))+
        theme_void()+
        scale_colour_manual(values=group.colors,guide="legend")+
        geom_text(aes( x=(as.numeric(plotframe[10,1])/2), y=0, label= paste0(z2,t,sep="")),family="Arial",color="deeppink4", size=8, fontface="bold")+
        geom_text(aes( x=(as.numeric(plotframe[10,1])/2), y=0.1, label= paste0(stops$strand[which(stops$ensembl_transcript_id==x & stops$frame == y)]," strand gene. LOEUF Decile = ",stops$LOEUF_decile[which(stops$ensembl_transcript_id==x & stops$frame==y)],sep="")),family="Arial",color="grey10", size=6)+
        geom_text(aes( x=(as.numeric(plotframe[10,1])/2), y=0.2, label= paste0("CDS Extension = ",format(round(stops$percent_protein_gain[which(stops$ensembl_transcript_id==x & stops$frame == y)], 0), nsmall = 1), "% / ",trunc(stops$aa_extension[which(stops$ensembl_transcript_id==x & stops$frame == y)]),"aa",sep="")),family="Arial",color="black", size=6)+ 
        geom_text(aes( x=(as.numeric(plotframe[10,1])/2), y=0.3, label= paste0(utr_tbl$MANE[which(utr_tbl$ENS_Transcript_id==x & utr_tbl$frame == y)]," ",stops$polyA_flag[which(stops$ensembl_transcript_id==x & stops$frame == y)],sep="")),family="Arial",color="grey10", size=6)+
        geom_text(aes( x=0.05, y=0.48, label="Stop Lost"),family="Arial",color="black", size=6)+
        geom_text(aes( x=0.05, y=0.78, label="Wild Type"),family="Arial",color="black", size=6)+
        geom_text(aes( x=(as.numeric(plotframe[10,1])/2), y=0.8, label= x),family="Arial",color="black", size=8, fontface="bold")+
        geom_text(aes( x=(as.numeric(plotframe[10,1])/2), y=0.9, label= paste0(stops$gene_symbol[which(stops$ensembl_transcript_id==x & stops$frame==y)],frametext,sep="")),family="Arial",color="black", size=8, fontface="bold")+
        theme(legend.key.height=unit(5,"lines"),legend.position = "bottom",legend.text=element_text(size=16),legend.key.width =unit(4,"cm"))+
        guides(linewidth=FALSE, fill = guide_legend(override.aes = list(shape = 15, size = 5)))+
        guides(color = guide_legend(reverse = TRUE))+
        labs(colour = "")
    }
    if(is.na(stops$rna_start[which(stops$ensembl_transcript_id==x & stops$frame==y)])){
      lab<-paste0("No Stop has been identified.", 
                  "We have not been able to identify a downstream stop in this frame.
      It is possible that this variant would lead to non-stop decay")
      p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
    }
  }
  if(t=="SNV"){
    if(z=="Stop preserving "){
      lab<-paste0("This SNV results in the formation of an alternative stop motif.
      It is unlikely that this variant would disrupt the stop codon.
                  If you feel this is an error please check your entry format, or contact us if the issue persists")
      p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
    }
    if(z != "Stop preserving "){
      if(is.na(stops$rna_start[which(stops$ensembl_transcript_id==x & stops$frame==y)])){
        lab<-paste0("No alternative stops in this frame were found within the mRNA in this transcript. 
                  It is likely that the consequence of a distrupted stop in transcript "
                    ,x," is non-stop decay")
        p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
      }
      if(stops$polyA_flag[which(stops$ensembl_transcript_id==x & stops$frame==y)]=="Likely Polyadenylation related NSD"){ #12:8941818_C/T
        lab<-paste0("The Primary Polyadenylation signal is prior to the first alternative stop codon in this frame in this transcript.",'\n',"It is likely that the consequence of a distrupted stop in transcript "
                    ,x," is non-stop decay")
        p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
      }
    }
  }
  if(t!="SNV"){
    if(z == "0"){
      lab<-"This MNV results in the formation of an alternative stop motif. 
      It is unlikely that this variant would disrupt the stop codon.
                If you feel this is an error please check your entry format, or contact us if the issue persists"
      p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
    }
    if(z != "0" & z != "Stop disrupting "){
      lab<-paste0("This MNV creates a new in-frame stop motif.",sep="")
      p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
    } 
    if(z=="intron"){
      lab<-paste0("WARNING Intron:Exon junction variant",
                  "This variant spans one or more intron-exon junction and may affect splicing in addition to stop disruption. 
                We would recommend running the variant through spliceAI Lookup to determine if splicing impact could supplement interpretation")
      p<- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()
    }
  }
  return(p)
}


############################################################################################################################################

### LAYOUT ###

#Populate URL text
url1 <- a("Recommendations for clinical interpretation of variants found in non-coding regions of the genome", href="https://link.springer.com/article/10.1186/s13073-022-01073-3")
url2 <- a("Systematic identification of disease-causing promoter and untranslated region variants in 8,040 undiagnosed individuals with rare disease", href="https://www.medrxiv.org/content/10.1101/2023.09.12.23295416v1")
url3 <- a("Loss of function SMAD4 nonstop mutations in human cancer", href="https://onlinelibrary.wiley.com/doi/full/10.1111/his.14880")
url4 <- a("website", href="https://rarediseasegenomics.org")
url5 <- a("GitHub repo", href="https://github.com/Computational-Rare-Disease-Genomics-WHG/inframe_stop_app")
url6 <- a("PolyA DB", href="https://exon.apps.wistar.org/PolyA_DB/")
url7 <- a("MANE", href="https://www.ncbi.nlm.nih.gov/refseq/MANE/")
url8 <- a("LOEUF", href="https://www.nature.com/articles/s41586-020-2308-7")

#Define user interface
ui <- fluidPage(
  theme = shinytheme("lumen"), # https://rstudio.github.io/shinythemes/
  titlePanel(title=div(h1("the STOP g.APP [BETA]", style = "font-size:40px;color:navy;", img(src="dna.png", height = 100, width = 110, style="position:absolute;right:40px;z-index:1000000;")))),
  sidebarLayout(
    sidebarPanel(
      h2("Data Selection Options"),
      #h5("If output doesn't update automatically - hit submit!"),
      h3("Search by variant position"),
      textInput(inputId = "free_var",label = "Input Variant position (i.e. cHr5:88722604_T/G)",placeholder = "e.g cHr5:88722604_T/G"),
      selectInput("drop_var","Select transcript",choices = c("First enter variant position"='')),
      actionButton("goVarButton", "   Submit variant", class = "btn-warning", width = '100%', class="btn-xs",style="color: #fff; background-color: #355BD6; border-color: #FAFAFA",icon=icon("dna")),
      h3("Search by gene"),
      textInput(inputId = "free_gene",label = "Input gene ID (ENS/Stable) - Note: Greek letters should be converted to English form e.g. β = B",placeholder = "e.g MEF2C"),
      selectInput("drop_transcript","Select transcript",choices = c("First enter gene ID"='')),
      actionButton("goGeneButton", "   Submit Gene", class = "btn-primary", width = '100%', class="btn-xs",style="color: #fff; background-color: #355BD6; border-color: #FAFAFA",icon=icon("dna")),
      h3("OR Search directly by Transcript"),
      textInput(inputId = "free_transcript",label = "Input transcript ID (ENS)",placeholder = "e.g ENST00000340208"),
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(strong("Consequence plot"), plotOutput("UTR")),
                  tabPanel(strong("Consequence information"), tableOutput("tabledata"),),
                  tabPanel(strong("About"), h3("About the Rescue Stop Site"),
                           tagList("The disruption of stop codons has been found to cause disease, but the interpretation of possible consequences of such disruptions is often difficult.One way of identifying the potential deleteriousness of a disrupted stop codon is to determine the extent to which the CDS has been extended, and the 3’UTR expunged."),
                           h5(""),
                           tagList("The effects of ntroducing additional amino acids may vary between transcripts, however, some key outcomes include protein misfolding and the potential for changes in phobicity of the protein product"),
                           h5(""),
                           tagList("The 3’ UTR of a given gene includes important regulatory elements, including polyA sies and signals, downstream open reading frames (dORFs), and RNA binding protein binding sites. When these elements are disrupted it has been shown to lead to disease phenotypes."),
                           h5(""),
                           tagList("The Rescue Stop Site provides the location of the next available in-frame stop codon, and both the extension of CDS that would be created if it were to be used, and the proportion of 3’ UTR that would be lost. This is currently limited to Ensembl transcripts (release 110). Read more about the role of 3' UTRs, regulatory variants in disease, and the consequences of stop loss in the following papers:"),
                           h5(""),
                           tagList("1 : Ellingford et al, ", url1, ", Genome Medicine (2023)."),
                           h5(""),
                           tagList("2 : Martin-Geary et al, ", url2, ", medRxiv (2023)."),
                           h5(""),
                           tagList("3 : Bauer et al (2023),", url3, ", Histopathology (2023)."),
                           h3("About the Team"),
                           h5(""),
                           tagList("The Rescue Stop Site was developed by the Computational Rare Disease Genomics Lab (CRDG). CRDG specialises in identifying the role of non-protein-coding variants in rare genetic diseases through the use of bioinformatics on large exome and genome sequencing datasets. We are based out of the Big Data Institute at the University of Oxford. You can learn more about our lab at our ",url4,"."),
                           h3("Tools and Data"), 
                           h5(""),
                           tagList("The Rescue Stop Site was built using R version 4.3.1 (2023-06-16), R studio version 2023.09.0+463, and the packages ‘base’, ‘utils’, ‘ggplot2’, ’shiny’, 'shinyjs', ’shinytheme’, ’stats’, ’dplyr’, ’tidyr’, ’stringr’, ’beepr’ and ‘Biostrings’. Stop codon and sequence data was collected using the ‘biomaRt’ R package and Ensembl release 110. This was further supplemented using polyadenylation data from ",url6,", alongside ",url7," transcript data, and ",url8," scores."),
                           h5(""),
                           tagList("All code and data is available through our GitHub repo <insert here> If you face any issues or wish to report a bug whilst using the Rescue Stop Site please get in touch with us by opening up an issue in our repository's ", url5),
                           h3("Acknowledgements"),
                           h5(""),
                           tagList("This work was funded through Wellcome Trust and the Royal Society (Sir Henry Dale Fellowship awarded to Nicola Whiffin; 220134/Z/20/Z) and the Rosetrees Trust (H5R01320). The ShinyTheme \"lumen\" was used for generating the template for this application. Our DNA Logo icon was created by  Freepik - Flaticon. All rights reserved for their respective owners."),),
                  tabPanel(strong("Usage"), h3("How to use the Rescue Stop Site"),
                           h3("*BETA PHASE*"),
                           tagList("the STOP g.APP is currently open for beta testing while we implement some additional features. We would appreciate any and all feedback, comments or suggestions you may have during this time!"),
                           h3("Variant queries"),
                           tagList("Variants should be input in the format chr[N]:[POSITION]_REF/ALT, i.e. chr19:45163192_T/C or 19:45163192_T/CGG. For single nucleotide polymorphisms the STOP gAPP can identify if the variant would disrupt, or preserve (i.e. result in an alternate stop motif) the stop codon. For multi nucleotide variants, alternate stops are calculated with reference to the frame change (if any) resulting from a frameshift due to the variant."),
                           h5(""),
                           h3("Gene level queries"),
                           tagList("At a gene level, the STOP gAPP will take in a gene symbol, and show the next in-frame stop downstream from the canonical stop for each transcript associated with that gene. Simply input a gene symbol i.e. MEF2c, TGFB1, tp53, and select a transcript of interest from the drop down menu."),
                           h5(""),
                           h3("Transcript specific queries"),
                           tagList("For queries relating to a particular transcript, type or paste the Ensembl transcript I.D of the transcript of interest into the box, and the next in-frame stop will be shown."),
                           h3("The visualisation tab"),
                           tagList("The visualisation tab will give you an overview of the impact of a stop-loss variant, including an illustration of the wild-type CDS length, and any additions to the CDS that would occur as a result of a downstream stop being used"),
                           h3("The consequence information tab"),
                           tagList("The consequence information tab provides a table view of the impact of a stop-loss variant. This includes general information about the transcript in question, alongside specific information about the predicted stop, including the downstream stop motif, genomic coordinates, and extension statistics."),)
      )
    )
  )
)


############################################################################################################################################

### MAIN ###
server <- function(input, output, session) {
  myframe<-NULL
  dis_pres<-""
  var_type<-""
  var2<-""
  ##### VARIANT CAPTURE:
  observe({
    var <- input$free_var
    var <- str_trim(var, side = "both")
    var2<<-var
    #If User has entered a variant
    if(!var=="" ){ 
      if(!grepl("^[a-zA-Z0-9]+[:][0-9]+[_][A-Za-z]+[/][A-Za-z]+",var)){
        lab<-"Incorrect variant format"
        output$UTR <- renderPlot({ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()})
      }
    }
    #If User has not entered a variant
    if(var==""){ 
      lab<-"Hello you ;)
      <- Please choose a search option."
      output$UTR <- renderPlot({ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = lab) + theme_void()})
    }
    #If the variant is in the correct format 'chrom:N_pos_ref/alt'
    if(grepl("^[a-zA-Z0-9]+[:][0-9]+[_][A-Za-z]+[/][A-Za-z]+",var)){
      choice_list<-get_transcripts(var)#Get the transcript variant is in, or warn "Not a known stop location" 
      print(choice_list)
      var_type<<-get_var_type(var)
      myframe<<-get_frame(var)
      dis_pres<-""
      
      #If transcripts corresponding with the variant have been found, get the desired transcript
      if(length(choice_list)!=0){ 
        lab<-paste0("Please select a transcript from the dropdown.")
        output$UTR <- renderPlot({
          ggplot() + annotate("text", x = 10,  y = 10,
                              size = 6,
                              label = lab) + theme_void()
        })
        updateSelectInput(session, "drop_var",
                          label = paste("Transcripts with stops intersecting ", var),
                          choices = choice_list)
      }
      #If transcripts corresponding with the variant have NOT been found Warn that this is not a stop intersecting variant
      if(length(choice_list)==0){
        updateSelectInput(session, "drop_var",
                          label = paste("Transcripts with stops intersecting ", var),
                          choices = choice_list)
        showModal(modalDialog(
          title = "Not a known stop location!",
          "This position is not recognised as being in a known stop codon.
              If you feel this is an error please check your entry format, or contact us if the issue persists",
          easyClose = TRUE,
          footer = modalButton("Dismiss")
        ))}
      }
  })
  
  #Plot the results
  observeEvent(input$goVarButton,{
    output$tabledata <- renderTable({
      tabledata<-stops
      x1<-c("Ensembl transcript ID","Gene symbol","Rescue stop motif","Chromosome","Position in DNA: start","Strand","Wild-type CDS length (amino acids)","Wild-type CDS length (nucleotides)","Extension length (amino acids)","Extension length (nucleotides)","Protein gain (%)","Frame","MANE status","Polyadenylation flags","LOEUF decile","Additional stop details")
      x2<-tabledata[which(tabledata$ensembl_transcript_id == tran_strip(input$drop_var) & tabledata$frame == myframe),c(1,7,4,5,10,6,13,12,15,14,20,16,21,23,24,11)]
      tabledata<-data.frame("Features"= as.character(x1), "details" = as.character(x2))
    })
    int_var<<-get_intron_vars(var2,tran_strip(input$drop_var))
    dis_pres<<-get_mnv_impact(var2,tran_strip(input$drop_var))
    if(int_var=="intron" & var_type %in% c("Deletion","MNV") & dis_pres != 0){ 
      strarray<-var_split(var2)
      sAI_url<-paste0("https://spliceailookup.broadinstitute.org/#variant=",strarray[1,1],"-",strarray[2,1],"-",strarray[3,1],"-",strarray[4,1],"&hg=38&distance=500&mask=0&ra=0",sep="")
      showModal(modalDialog(
        title = "WARNING Intron:Exon junction variant",
        tags$div("This variant spans one or more intron-exon junction and may affect splicing in addition to stop disruption. Due to this, frame shifting may not be accurate - plot view is therefore a guide only. We would recommend running the variant through ", tags$a(href = sAI_url, "spliceAI")," to determine if splicing impact could supplement interpretation"),
        easyClose = TRUE,
        footer = modalButton("Dismiss")
      ))
      output$UTR <- renderPlot({get_plot(tran_strip(input$drop_var),myframe,var_type,var2)})
    }
    if(int_var=="not_intron"){
      if(get_var_type(var2)=="SNV"){
        if(get_pres(var2)=="Stop preserving "){
          showModal(modalDialog(
            title = "STOP PRESERVING VARIANT",
            tags$div("This SNV results in the formation of an alternative stop motif. It is unlikely that this variant would disrupt the stop codon. If you feel this is an error please check your entry format, or contact us if the issue persists"),
            easyClose = TRUE,
            footer = modalButton("Dismiss")
          ))
          output$UTR <- renderPlot({get_plot(tran_strip(input$drop_var),myframe,var_type,var2)})
        }
        if(get_pres(var2) != "Stop preserving "){
          if(is.na(stops$rna_start[which(stops$ensembl_transcript_id==tran_strip(input$drop_var) & stops$frame==get_frame(var2))])){
            showModal(modalDialog(
              title = "STOP DISRUPTING VARIANT",
              tags$div(paste0("No alternative stops in this frame were found within the mRNA in this transcript.It is likely that the consequence of a distrupted stop in transcript ",tran_strip(input$drop_var)," is non-stop decay")),
              easyClose = TRUE,
              footer = modalButton("Dismiss")
            ))
            output$UTR <- renderPlot({get_plot(tran_strip(input$drop_var),myframe,var_type,var2)})
            }
          if(stops$polyA_flag[which(stops$ensembl_transcript_id==tran_strip(input$drop_var) & stops$frame==get_frame(var2))]=="Likely Polyadenylation related NSD"){ 
            showModal(modalDialog(
              title = "POLYA TRIGGERED NSD",
              tags$div(paste0("The Primary Polyadenylation signal is prior to the first alternative stop codon in this frame in this transcript.",'\n',"It is likely that the consequence of a distrupted stop in transcript ",tran_strip(input$drop_var)," is non-stop decay")),
              easyClose = TRUE,
              footer = modalButton("Dismiss")
            ))
            output$UTR <- renderPlot({get_plot(tran_strip(input$drop_var),myframe,var_type,var2)})
          }
        }
      }
      if(get_var_type(var2)!="SNV"){
        if(dis_pres == "0"){
          showModal(modalDialog(
            title = "STOP PRESERVING MNV",
            tags$div(paste0("This MNV results in the formation of an alternative stop motif in the same location as the wild-type stop. It is unlikely that this variant would disrupt the stop codon. If you feel this is an error please check your entry format, or contact us if the issue persists")),
            easyClose = TRUE,
            footer = modalButton("Dismiss")
          ))
          output$UTR <- renderPlot({get_plot(tran_strip(input$drop_var),myframe,var_type,var2)})
        }
        if(dis_pres != "Stop disrupting " & dis_pres != "0"){
          showModal(modalDialog(
            title = "STOP PRESERVING MNV",
            tags$div(paste0("This MNV creates a new in-frame stop motif. It alters the length of the peptide by ",dis_pres," amino acids",sep="")),
            easyClose = TRUE,
            footer = modalButton("Dismiss")
          ))
          output$UTR <- renderPlot({get_plot(tran_strip(input$drop_var),myframe,var_type,var2)})
        } 
      }
      output$UTR <- renderPlot({get_plot(tran_strip(input$drop_var),myframe,var_type,var2)})
    }
  })
  
  ##### GENE CAPTURE:
  observe({
    gene <- toupper(input$free_gene)
    gene<-str_trim(gene, side = "both")
    var<-""
    if (gene==""){ # Can use character(0) to remove all choices
      gene <- character(0)
    }
    if(!(identical(gene, character(0)))){
      if(!gene %in% stops$gene_symbol){
        lab<-paste0("Gene ",gene," not found.")
        output$UTR <- renderPlot({
          ggplot() + annotate("text", x = 10,  y = 10,
                              size = 6,
                              label = lab) + theme_void()
        })
      }
      if(gene %in% stops$gene_symbol){
        lab<-paste0("Please select a transcript from the dropdown.")
        output$UTR <- renderPlot({
          ggplot() + annotate("text", x = 10,  y = 10,
                              size = 6,
                              label = lab) + theme_void()
        })
        choice_list<-unique(ref_stops$mane_transcript[which(ref_stops$gene_symbol == gene)])
        myframe<<-0
        updateSelectInput(session, "drop_transcript",
                          label = paste("Transcripts by gene ", gene),
                          choices = choice_list
        )
      }
    }
  })
  observeEvent(input$goGeneButton,{
    output$tabledata <- renderTable({
      tabledata<-stops
      x1<-c("Ensembl transcript ID","Gene symbol","Rescue stop motif","Chromosome","Position in DNA: start","Strand","Wild-type CDS length (amino acids)","Wild-type CDS length (nucleotides)","Extension length (amino acids)","Extension length (nucleotides)","Protein gain (%)","Frame","MANE status","Polyadenylation flags","LOEUF decile","Additional stop details")
      x2<-tabledata[which(tabledata$ensembl_transcript_id == tran_strip(input$drop_transcript) & tabledata$frame == myframe),c(1,7,4,5,10,6,13,12,15,14,20,16,21,23,24,11)]
      tabledata<-data.frame("Features"= as.character(x1), "details" = as.character(x2))
    })
    print(myframe)
    var_type<-""
    var<-""
    output$UTR <- renderPlot({get_plot(input$drop_transcript,myframe,var_type,var)})
  })
  
  ##### TRANSCRIPT CAPTURE
  observe({
    var<-""
    transcript <- toupper(input$free_transcript)
    if(nchar(transcript)==15){
      myframe<<-0
      output$tabledata <- renderTable({
        tabledata<-stops
        x1<-c("Ensembl transcript ID","Gene symbol","Rescue stop motif","Chromosome","Position in DNA: start","Strand","Wild-type CDS length (amino acids)","Wild-type CDS length (nucleotides)","Extension length (amino acids)","Extension length (nucleotides)","Protein gain (%)","Frame","MANE status","Polyadenylation flags","LOEUF decile","Additional stop details")
        x2<-tabledata[which(tabledata$ensembl_transcript_id == transcript & tabledata$frame == myframe),c(1,7,4,5,10,6,13,12,15,14,20,16,21,23,24,11)]
        tabledata<-data.frame("Features"= as.character(x1), "details" = as.character(x2))
      })
      output$UTR <- renderPlot({get_plot(input$free_transcript,myframe,var_type,var)})
    }
  })
}

shinyApp(ui = ui, server)

#END#