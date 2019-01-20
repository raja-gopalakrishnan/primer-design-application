library(shiny)
library(XML)
library(Biostrings)
library(stringi)
library(stringr)

path="path/to/primer_design_app"

shinyui = fluidPage(
  titlePanel("Design Primers"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene","Gene name",""),
      selectInput("expt", "Experiment:",c("Deletion"="deletion", "Tagging"="tagging")),
      numericInput("homology","Number of homologous nucleotides:",40),
      conditionalPanel(
        condition = "input.expt == \"tagging\"",
        selectInput("tag","Tag:",c("HA","GFP","myc","Flag"))
      ),
      uiOutput("selection"),
      actionButton("submit","Submit")
    ),
    mainPanel(
      fluidRow(
        column(width=2,("FB plasmid: ")),
        column(width=3,textOutput("fb_number"))
      ),
      br(),
      tableOutput(outputId = "primers"), br(),
      "Checking primer:", tableOutput("checking"), br(),
      "Gene sequence:", 
      tags$div(HTML(paste("<pre><span>",htmlOutput("seq1",inline=T),"</span><span style=\"color:red\">",
                          htmlOutput("seq2",inline=T),"</span><span>",
                          htmlOutput("seq3",inline=T),"</span><span style=\"color:blue\">",
                          htmlOutput("seq4",inline=T),"</span><span style=\"color:#E1E100\">",
                          htmlOutput("seq5",inline=T),"</span><span style=\"color:blue\">",
                          htmlOutput("seq6",inline=T),"</span><span>",htmlOutput("seq7",inline=T),"</span></pre>",sep="")))
    )
  )
)


shinyServer = function(input,output)
{
  primer_5_HA="CGGATCCCCGGGTTAATTAA"
  primer_3_HA="GAATTCGAGCTCGTTTAAAC"
  sgd_ids=read.table(paste(path,"sgd_ids.tsv",sep=""),sep="\t",stringsAsFactors = F,head=F)
  plasmids=read.table(paste(path,"plasmids.txt",sep=""),stringsAsFactors = F, head=F)
  primer_table=NULL
  
  #Dynamic UI based on selected tag
  output$selection=renderUI({
    if(input$expt=="deletion"){selectInput("sel_cassette","Selection cassette:",c("kanMX","hygMX","natMX"))}
    else if(input$expt=="tagging" && input$tag=="HA"){selectInput("sel_cassette","Selection cassette:",c("kanMX","natMX"))}
    else if(input$expt=="tagging" && input$tag=="GFP"){selectInput("sel_cassette","Selection cassette:",c("kanMX"))}
    else if(input$expt=="tagging" && input$tag=="myc"){selectInput("sel_cassette","Selection cassette:",c("kanMX","natMX"))}
    else if(input$expt=="tagging" && input$tag=="Flag"){selectInput("sel_cassette","Selection cassette:",c("kanMX"))}
  })
  
  #Which plasmid to use based on expt, tag and selection cassette
  txt=eventReactive(input$submit,{
    if(input$expt=="deletion")
    {
      var1 = which(plasmids[,1]==input$expt)
      var2 = which(plasmids[var1,2]==input$sel_cassette)
      plasmids[var1,][var2,3]
    }
    else if(input$expt=="tagging")
    {
      var1 = which(plasmids[,1]==input$tag)
      var2 = which(plasmids[var1,2]==input$sel_cassette)
      plasmids[var1,][var2,3]
    }
  })
  output$fb_number = renderText(txt())
  
  #Extract primer sequence
  tab=eventReactive(input$submit,{
    i = toupper(input$gene)
    
    #Extract the SGD ID for the gene. This is used in the link
    sgd = sgd_ids[which(match(sgd_ids[,3],i)==T),1]
    if(identical(sgd,character(0)))
    {
      sgd = sgd_ids[which(match(sgd_ids[,2],i)==T),1]
    }
    
    #Download the HTML file containing the FASTA sequence of the gene +- 1kb
    link = paste("http://yeastgenome.org/cgi-bin/getSeq?map=a3map&seq=",sgd,"&flankl=1000&flankr=1000&rev=",sep="")
    file_name = paste(path,"webpages/",i,"_flanking.html",sep="")
    download.file(link,file_name)
    
    #Read the HTML file and extract the sequence from the file
    doc.html = htmlParse(file_name)
    gene_sequence = unlist(xpathApply(doc.html, '//pre', xmlValue))
    
    #Get rid of extraneous information
    cut=regexpr("\r",gene_sequence)[1]
    gene_sequence=substr(gene_sequence,cut,stri_length(gene_sequence))
    gene_sequence=gsub("\r\n","",gene_sequence)
    
    #Set start and stop codon coordinates
    length=stri_length(gene_sequence)
    start=c(1001,1003)
    stop=c(length-1002,length-1000)
    
    #Checking primers
    if(input$expt=="deletion"){coord=start}
    else if(input$expt=="tagging"){coord=stop}
    check_table_5prime=NULL
    for(j in (coord[1]-310):(coord[1]-65))
    {
      test_seq=substr(gene_sequence,j,j+19)
      G_count=str_count(test_seq,"G")
      C_count=str_count(test_seq,"C")
      A_count=str_count(test_seq,"A")
      T_count=str_count(test_seq,"T")
      gc_content=(G_count+C_count)/20*100
      tm=round(4*(G_count+C_count) + 2*(20-G_count-C_count),0)
      
      check_table_5prime=rbind(check_table_5prime,as.data.frame(cbind.data.frame(test_seq,gc_content,tm)))
    }
    closest_5prime = tail(which(check_table_5prime[,3]==max(check_table_5prime[,3])),1)
    check_5primer = cbind(i,check_table_5prime[closest_5prime,])
    colnames(check_5primer)=c("Gene","Sequence","GC content (%)","Melting temp")
    
    #Extract 5' and 3' flanking sequences (40 bp)
    hom=input$homology
    flank_5_del=substr(gene_sequence,start[1]-hom,start[1]-1)
    flank_5_tag=substr(gene_sequence,stop[1]-hom,stop[1]-1)
    flank_3=as.character(reverseComplement(DNAString(substr(gene_sequence,stop[2]+1,stop[2]+hom))))
    
    #Make the final oligo by pasting the flanking sequence with the priming sequence
    if(input$expt=="tagging")
    {
      final_oligo_5=(paste0(flank_5_tag,primer_5_HA))
    }
    else if(input$expt=="deletion")
    {
      final_oligo_5=(paste0(flank_5_del,primer_5_HA))
    }
    final_oligo_3=(paste0(flank_3,primer_3_HA))
    
    #Make the table with Gene name, primer (5' or 3'), and sequence
    primer_table=rbind(primer_table,c(i,paste("5' ",input$expt ," primer",sep=""),final_oligo_5))
    primer_table=rbind(primer_table,c(i,paste("3' ",input$expt ," primer",sep=""),final_oligo_3))
    colnames(primer_table)=c("Gene","Primer","Sequence")
    
    #Return the different parameters, a is the primer table, b is the checking primer, 
    #c is the gene sequence before checking primer, d is the checking primer sequence,
    #e is the gene sequence after checking primer, f is the 5' homology sequence, g is the deleted sequence,
    #h is the 3' homology sequence, i is the gene sequence after the 3' homology sequence
    list("a"=(primer_table),"b"=check_5primer,
         "c"=substring(gene_sequence,start[1]-500,coord[1]-312+closest_5prime),
         "d"=substring(gene_sequence,coord[1]-311+closest_5prime,coord[1]-292+closest_5prime),
         "e"=substring(gene_sequence,coord[1]-291+closest_5prime,coord[1]-(hom+1)),
         "f"=substring(gene_sequence,coord[1]-hom,coord[1]-1),
         "g"=substring(gene_sequence,coord[1],stop[2]),
         "h"=substring(gene_sequence,stop[2]+1,stop[2]+hom),
         "i"=substring(gene_sequence,stop[2]+(hom+1),stop[2]+500)
    )
  })
  output$primers = renderTable(tab()$a)
  output$checking = renderTable(tab()$b)
  
  letters=c("c","d","e","f","g","h","i")
  for(k in 1:length(letters))
  {
    eval(parse(text=paste("output$seq",k," = renderText(tab()$",letters[k],")",sep="")))
  }
}

shinyApp(ui=shinyui,server=shinyServer)
