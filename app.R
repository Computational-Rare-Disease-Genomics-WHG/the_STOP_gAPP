library(shiny)
library(shinythemes)
library(ggplot2); theme_set(theme_linedraw())
library(utils)
library(datasets)
library(shinyjs)
library(png)


stops<-read.csv(file = "All_ENS_rescue_stops.txt", sep = '\t', header = TRUE, na.strings=c("","NA"))#[c(2:18,1)]
stops$strand[which(stops$strand == "-")]<-"Negative"
stops$strand[which(stops$strand == "+")]<-"Positive"
utr_tbl<-read.csv(file = "plot_stops.txt", sep = '\t', header = TRUE, na.strings=c("","NA"))#[c(2:18,1)]
url1 <- a("Recommendations for clinical interpretation of variants found in non-coding regions of the genome", href="https://link.springer.com/article/10.1186/s13073-022-01073-3")
url2 <- a("Systematic identification of disease-causing promoter and untranslated region variants in 8,040 undiagnosed individuals with rare disease", href="https://www.medrxiv.org/content/10.1101/2023.09.12.23295416v1")
url3 <- a("Loss of function SMAD4 nonstop mutations in human cancer", href="https://onlinelibrary.wiley.com/doi/full/10.1111/his.14880")
url4 <- a("website", href="https://rarediseasegenomics.org")
url5 <- a("GitHub repo", href="https://github.com/Computational-Rare-Disease-Genomics-WHG/inframe_stop_app")


ui <- fluidPage(
  theme = shinytheme("lumen"), # https://rstudio.github.io/shinythemes/
  titlePanel(title=div(h1("the Rescue Stop Site", style = "font-size:40px;color:navy;", 
                          img(src="dna.png", height = 100, width = 110,
                              style="position:absolute;right:40px;z-index:1000000;")))),
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      h2("Data Selection Options"),
      radioButtons("radio",
                   h6(""), 
                   choices = list( 
                                  "Disrupted stop codons" = 2),
                   selected = 2),
      h3("Search by gene"),
      textInput(inputId = "free_gene",label = "Input gene ID (ENS/Stable) - Note Greek letters should be converted to English form e.g. β = B",placeholder = "e.g MEF2C"),
      verbatimTextOutput("placeholder"),
      selectInput("drop_transcript","Select transcript",choices = c("First enter gene ID"='')),
      h3("OR Search directly by Transcript"),
      textInput(inputId = "free_transcript",label = "Input transcript ID (ENS)",placeholder = "e.g ENST00000340208"),
      h5("Note Some Ensembl transcripts are incomplete - due to this, a small number of entries will show as having amino acid counts that are not whole numbers"),
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(strong("Consequence information"), tableOutput("tabledata"),),
                  tabPanel(strong("Consequence plots"), plotOutput("UTR")),
                  tabPanel(strong("About"), h3("About the Rescue Stop Site"),
                           tagList("The disruption of stop codons has been found to cause disease, but the interpretation of possible consequences of such disruptions is often difficult.One way of identifying the potential deleteriousness of a disrupted stop codon is to determine the extent to which the CDS has been extended, and the 3’UTR expunged."),
                           h5(""),
                           tagList("<Note on effects of CDS extension>"),
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
                           tagList("The Rescue Stop Site was built using R version 4.3.1 (2023-06-16), R studio version 2023.09.0+463, and the packages ‘base’, ‘utils’, ‘ggplot2’, ’shiny’, 'shinyjs', ’shinytheme’, ’stats’, ’dplyr’, ’tidyr’, ’stringr’, ’beepr’ and ‘Biostrings’. Data was collected using the ‘biomaRt’ R package and Ensembl release 110."),
                           h5(""),
                           tagList("All code and data is available through our GitHub repo <insert here> If you face any issues or wish to report a bug whilst using the Rescue Stop Site please get in touch with us by opening up an issue in our repository's ", url5),
                           h3("Acknowledgements"),
                           h5(""),
                           tagList("This work was funded through Wellcome Trust and the Royal Society (Sir Henry Dale Fellowship awarded to Nicola Whiffin; 220134/Z/20/Z) and the Rosetrees Trust (H5R01320). The ShinyTheme \"lumen\" was used for generating the template for this application. Our DNA Logo icon was created by  Freepik - Flaticon. All rights reserved for their respective owners."),)
      )
    )
  )
)


server <- function(input, output, session) {
  dataset <- reactive({
    if(input$radio == 1){
      dataset <- starts
    }else if(input$radio == 2){
      dataset <- stops
    }
  })  
  
  observe({
    x <- input$free_gene
    
    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)
    
    # Can also set the label and select items
    updateSelectInput(session, "drop_transcript",
                      label = paste("Transcripts by gene ", x),
                      choices = dataset()$ENS_Transcript_id[which(dataset()$Gene_Symbol == x)]
    )
    
    
    output$tabledata <- renderTable({
      tabledata<-dataset()
      x1<-c("Ensembl transcript ID","Rescue stop motif","Chromosome","Position in RNA: start","Position in RNA: end","Position in DNA: start","Position in DNA: end","Strand","Wild-type CDS length (nucleotides)","Wild-type CDS length (amino acids)","Extension length (amino acids)","Extension length (nucleotides)","Percentage extended (amino acids)","Percentage extended (nucleotides)","Ensembl gene ID","Gene symbol","Percentage of UTR lost")
      x2<-tabledata[which(tabledata$ENS_Transcript_id == input$drop_transcript),]
      if(nchar(input$drop_transcript)>5){
        tabledata<-data.frame("Features"= as.character(x1), "details" = as.character(x2))
      }
    })
  })
  observe({
    x <- input$free_transcript
    if(nchar(x)>5){
      # Can use character(0) to remove all choices
      if (is.null(x))
        x <- character(0)
      
      output$tabledata <- renderTable({
        tabledata<-dataset()
        x1<-c("Ensembl transcript ID","Rescue stop motif","Chromosome","Position in RNA: start","Position in RNA: end","Position in DNA: start","Position in DNA: end","Strand","Wild-type CDS length (nucleotides)","Wild-type CDS length (amino acids)","Extension length (amino acids)","Extension length (nucleotides)","Percentage extended (amino acids)","Percentage extended (nucleotides)","Ensembl gene ID","Gene symbol","Percentage of UTR lost")
        x2<-tabledata[which(tabledata$ENS_Transcript_id == input$free_transcript),]
        if(nchar(input$free_transcript)>5){
          tabledata<-data.frame("Features"= as.character(x1), "details" = as.character(x2))
        }
      })
    }
  })
  observe({
    x <- input$drop_transcript
    y <- input$free_transcript
    if(nchar(x)>5){
      if(is.na(dataset()$start_rna[which(dataset()$ENS_Transcript_id == x)])){
        showModal(modalDialog(
          title = "NON-STOP DECAY!",
          "No inframe stops following the cannonical stop were found within the mRNA in this transcript. 
                It is likely that the consequence of a distrupted stop in this gene is non-stop decay",
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }
    if(nchar(y)>13){
      if(y %in% dataset()$ENS_Transcript_id){
        if(is.na(dataset()$start_rna[which(dataset()$ENS_Transcript_id == y)])){
          showModal(modalDialog(
            title = "NON-STOP DECAY!",
            "No inframe stops following the cannonical stop were found within the mRNA in this transcript. 
                  It is likely that the consequence of a distrupted stop in this gene is non-stop decay",
            easyClose = TRUE,
            footer = NULL
          ))
        }
      }
    }
  })
  observe({
    if(input$drop_transcript %in% utr_tbl$ENS_Transcript_id | input$free_transcript %in% utr_tbl$ENS_Transcript_id){
      output$UTR <- renderPlot({
        #Add plot with different weights of bar representing the proportional lengths of 5'UTR, 3'UTR and CDS limits for wild type and mutant transcript
        df2<-data.frame(x=1:10,y=1:10,z=c(1,1,1,2,2,2,3,3,4,4))
        ggplot(df2,aes(x,y,linewidth=as.numeric(z)))+geom_line()
        plotframe<-data.frame(prop_pos=numeric(0),
                              rowheight=numeric(0),
                              weight=numeric(0),
                              category=character(0),
                              region=character(0),
                              stringsAsFactors = FALSE)
        #Wild type
        plotframe[nrow(plotframe) + 1,] = c(0,0.7,8,"Wild Type","CDS")
        plotframe[nrow(plotframe) + 1,] = c(utr_tbl$Original_cds[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)],0.7,8,"Wild Type","CDS")
        plotframe[nrow(plotframe) + 1,] = c(utr_tbl$Original_cds[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)],0.7,2,"Wild Type","UTR")
        plotframe[nrow(plotframe) + 1,] = c(1,0.7,2,"Wild Type","UTR")
        #Stop lost
        plotframe[nrow(plotframe) + 1,] = c(0,0.4,8,"Stop lost","CDS")
        plotframe[nrow(plotframe) + 1,] = c(utr_tbl$Original_cds[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)],0.4,8,"Stop lost","CDS")
        plotframe[nrow(plotframe) + 1,] = c(utr_tbl$Original_cds[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)],0.4,8,"Stop lost","Extension")
        plotframe[nrow(plotframe) + 1,] = c((utr_tbl$Original_cds[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)]+utr_tbl$Proportion_extended[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)]),0.4,8,"Stop lost","Extension")
        plotframe[nrow(plotframe) + 1,] = c((utr_tbl$Original_cds[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)]+utr_tbl$Proportion_extended[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)]),0.4,2,"Stop lost","UTR")
        plotframe[nrow(plotframe) + 1,] = c(1,0.4,2,"Stop lost","UTR")
        
        group.colors <- c(CDS = "#ABABD4", Extension = "#FFBA0A", UTR ="#383971")
        ggplot(plotframe,aes(as.numeric(prop_pos),as.numeric(rowheight),linewidth=as.numeric(weight),group=category,colour=region))+
          geom_line()+
          scale_y_continuous(limits = c(0, 1))+
          theme_void()+
          scale_colour_manual(values=group.colors)+
          geom_text(aes( x=0.5, y=0.9, label= stops$Gene_Symbol[which(stops$ENS_Transcript_id==utr_tbl$ENS_Transcript_id[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)])]),color="#00665E", size=7, fontface="bold")+
          geom_text(aes( x=0.05, y=0.78, label="Wild Type"),color="black", size=5)+
          geom_text(aes( x=0.05, y=0.48, label="Stop Lost"),color="black", size=5)+
          geom_text(aes( x=0.05, y=0.05, label= paste0(stops$strand[which(stops$ENS_Transcript_id==utr_tbl$ENS_Transcript_id[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)])]," strand gene",sep="")),color="grey10", size=5,fontface="bold")+
          geom_text(aes( x=0.5, y=0.2, label= paste0("CDS Extension = ",format(round(stops$prop_AA_ext[which(stops$ENS_Transcript_id==utr_tbl$ENS_Transcript_id[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)])], 2), nsmall = 2), "% / ",stops$extension_AA[which(stops$ENS_Transcript_id==utr_tbl$ENS_Transcript_id[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)])],"aa",sep="")),color="black", size=5)+
          geom_text(aes( x=0.5, y=0.15, label= paste0("UTR loss = ",format(round((utr_tbl$extension_NT[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)]/utr_tbl$utr_length_3[which(utr_tbl$ENS_Transcript_id == input$free_transcript | utr_tbl$ENS_Transcript_id == input$drop_transcript)])*100, 2), nsmall = 2), "%",sep="")),color="black", size=5)+
          theme(legend.position = "none")
      })
    }
  })
}

shinyApp(ui = ui, server)

