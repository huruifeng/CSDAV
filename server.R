#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
options(shiny.maxRequestSize=1000*1024^2) 

# Define server logic required to draw a histogram
server = function(input, output, session) {
  ready <- reactiveValues(fastq=0, library=0,counts_path=0,library_path=0, annotated=0) # to control the appearance/dissapearance of widgets
  output$fastq <- reactive({ return(ready$fastq) }) # an output to detect the value of the flag in the UI
  output$counts_path <- reactive({ return(ready$counts_path) })
  output$library_path <- reactive({ return(ready$library_path) })
  output$annotated <- reactive({ return(ready$annotated) })
  
  outputOptions(output, name="fastq", suspendWhenHidden=FALSE) # a condition to use in the conditional panels of the UI
  outputOptions(output, name="counts_path", suspendWhenHidden=FALSE)
  outputOptions(output, name="library_path", suspendWhenHidden=FALSE)
  outputOptions(output, name="annotated", suspendWhenHidden=FALSE)
  
  
  source("R/readfastq.R")
  source("R/readcounts.R")
  source("R/quality.R")
  source("R/annotate.R")
  source("R/rra2tail.R")
  source("R/differential.R")
  
  source("R/mageckrra.R")
  
  
  ### PLOTTING FUNCTIONS ### ---------------------------------------------------
  
  # FASTQ level:
  
  plot.mapped1 <- function(anno){
    p1 <- ggplot2::ggplot(anno, ggplot2::aes_string("File","`Mapped%`")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::xlab("") +
      ggplot2::ylab("% mapped")
    return(p1)
  }
  
  plot.mapped2 <- function(anno){
    anno$Unmapped <- anno$Total_reads-anno$Mapped
    my_data <- reshape2::melt(anno[,c("File","Mapped","Unmapped")], id.vars="File", variable.name="Mapping", value.name="Depth")
    p2 <- ggplot2::ggplot(my_data, ggplot2::aes_string(x="File", y="Depth", fill="Mapping")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(
        legend.title=ggplot2::element_blank(),
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::xlab("") +
      ggplot2::ylab("No. of reads in the fastq file")
    return(p2)
  }
  
  # Read-count level:
  
  plot.classification <- function(id,counts) {
    if(!"Category"%in%colnames(id)) id$Category <- "All"
    id1 <- id[id$sgRNA%in%rownames(counts),c("sgRNA","Gene","Category")]
    tmp1 <- data.frame(table(id1$Category))
    tmp2 <- data.frame(table(unique(id1[,c("Gene","Category")])$Category))
    names(tmp1)[1] <- names(tmp2)[1] <- "Category"
    p1 <- ggplot2::ggplot(tmp1, ggplot2::aes_string(x="factor(1)", y="Freq", fill="Category")) +
      ggplot2::ggtitle("sgRNA level") +
      ggplot2::geom_bar(stat="identity") +
      ggplot2:: coord_polar(theta = "y") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(
        axis.text.y=ggplot2::element_blank(),
        axis.ticks.y=ggplot2::element_blank(),
        legend.title=ggplot2::element_blank(),
        text=ggplot2::element_text(face="bold",color="black",size=12)) +
      ggplot2::xlab("") +
      ggplot2::ylab("")
    p2 <- ggplot2::ggplot(tmp2, ggplot2::aes_string(x="factor(1)", y="Freq", fill="Category")) +
      ggplot2::ggtitle("Gene level") +
      ggplot2::geom_bar(stat="identity") +
      ggplot2:: coord_polar(theta="y") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(
        axis.text.y=ggplot2::element_blank(),
        axis.ticks.y=ggplot2::element_blank(),
        legend.title=ggplot2::element_blank(),
        text=ggplot2::element_text(face="bold",color="black",size=12)) +
      ggplot2::xlab("") +
      ggplot2::ylab("")
    tg1 = gridExtra::tableGrob(tmp1,rows=NULL)
    tg2 = gridExtra::tableGrob(tmp2,rows=NULL)
    return(gridExtra::grid.arrange(p1,tg1,p2,tg2,nrow=2))
  }
  
  plot.gini <- function(anno){
    colnames(anno)[1] <- "Sample"
    p1 <- ggplot2::ggplot(anno,ggplot2::aes_string("Sample","Gini")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
                     text=ggplot2::element_text(face="bold",color="black",size=15),
                     axis.title.x=ggplot2::element_blank()) +
      ggplot2::ylab("Gini index")
    return(p1)
  }
  
  plot.sizes <- function(library, counts_mat){
    if(!"Category"%in%colnames(library)) library$Category <- "All"
    lib1 <- library[library$sgRNA%in%rownames(counts_mat),c("sgRNA","Gene","Category")]
    tmp1 <- split(lib1,lib1$Category)
    tmp2 <- do.call(rbind,lapply(1:length(tmp1),function(x){
      out <- data.frame(table(table(tmp1[[x]]$Gene)))
      names(out) <- c("Content","Frequency")
      out$Category <- names(tmp1)[x]
      return(out)
    }))
    tmp2$Content <- as.character(tmp2$Content)
    p1 <- ggplot2::ggplot(tmp2, ggplot2::aes_string(x="Content", y="Frequency", fill="Category")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        legend.position="top",
        legend.title=ggplot2::element_blank(),
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::xlab("sgRNAs per gene") +
      ggplot2::ylab("No.of genes")
    # ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
    return(p1)
  }
  
  plot.depths <- function(counts){
    if(!"Category"%in%colnames(counts)) counts$Category <- "All"
    tmp3 <- do.call(rbind,lapply(split(counts[,-which(colnames(counts)%in%c("sgRNA","Gene","Category")),drop=F],counts$Category),colSums))
    tmp3 <- reshape2::melt(tmp3,value.name="Depth",varnames=c("Category","Sample"))
    p1 <- ggplot2::ggplot(tmp3, ggplot2::aes_string(x="Sample", y="Depth", fill="Category")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        legend.title=ggplot2::element_blank(),
        axis.title.x=ggplot2::element_blank(),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::ylab("No. of counts")
    return(p1)
  }
  
  plot.zeros <- function(counts){
    if(!"Category"%in%colnames(counts)) counts$Category <- "All"
    tmp3 <- do.call(rbind,lapply(split(counts[,-which(colnames(counts)%in%c("sgRNA","Gene","Category")),drop=F],counts$Category), function(x) apply(x,2,function(y) sum(y==0))))
    tmp3 <- reshape2::melt(tmp3,value.name="Zeroes",varnames=c("Category","Sample"))
    p1 <- ggplot2::ggplot(tmp3, ggplot2::aes_string(x="Sample", y="Zeroes", fill="Category")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme_linedraw() +
      ggplot2::theme(
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        legend.title=ggplot2::element_blank(),
        axis.title.x=ggplot2::element_blank(),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::ylab("No. of zeros")
    return(p1)
  }
  
  plot.correlation <- function(counts,method){
    tmp1 <- stats::cor(counts,method=method)
    return(gplots::heatmap.2(tmp1, cexRow=1, cexCol=1, margins=c(15,15), #cexRow and cexCol change the font size
                             key = FALSE, dendrogram = "row", lhei = c(0.1,1), lwid=c(0.1,1),
                             trace="none", density.info="none", #show trace lines?
                             cellnote=round(tmp1,2), notecol="white", notecex=1.2,
                             col=grDevices::colorRampPalette(c("black","blue"))))
  }
  
  plot.pca <- function(pca,Y) {
    p1 <- ggplot2::ggplot(pca, ggplot2::aes_string("Z1",Y,label="Sample")) +
      ggrepel::geom_text_repel() +
      ggplot2::geom_point() +
      ggplot2::ggtitle("2nd vs 1st") +
      ggplot2::xlab("First Principal Component for Samples") +
      ggplot2::ylab("Second Principal Component") +
      ggplot2::theme_linedraw()
    p2 <- ggplot2::ggplot(pca, ggplot2::aes_string("Z1",Y,label="Sample")) +
      ggrepel::geom_text_repel() +
      ggplot2::geom_point() +
      ggplot2::ggtitle("3rd vs 1st") +
      ggplot2::xlab("First Principal Component for Samples") +
      ggplot2::ylab("Third Principal Component") +
      ggplot2::theme_linedraw()
    p3 <- gridExtra::grid.arrange(p1,p2,nrow=1)
    return(p3)
  }
  
  plot.pve <- function(pve_sample) {
    pve <- data.frame(PC=seq(1:length(pve_sample)),pve_sample=pve_sample)
    cpve <- data.frame(PC=seq(1:length(pve_sample)),Cpve_sample=cumsum(pve_sample))
    p1 <- ggplot2::ggplot(pve, ggplot2::aes_string("PC","pve_sample")) +
      ggplot2::geom_point(col="darkblue",shape="O",size=5) +
      ggplot2::geom_line(col="blue") +
      ggplot2::xlab("Principal Component") +
      ggplot2::ylab("Proportion of Variance Explained") +
      ggplot2::theme_linedraw()
    p2 <- ggplot2::ggplot(cpve, ggplot2::aes_string("PC","Cpve_sample")) +
      ggplot2::geom_line(col="blue") +
      ggplot2::geom_point(col="darkblue",shape="O",size=5) +
      ggplot2::xlab("Principal Component") +
      ggplot2::ylab("Cumulative PVE") +
      ggplot2::theme_linedraw()
    p3 <- gridExtra::grid.arrange(p1,p2,nrow=1)
    return(p3)
  }
  
  plot.folds <- function(folded, loglevel) {
    if(input$lfc.style=="Box plot") {
      p1 <- ggplot2::ggplot(folded,ggplot2::aes_string("Category","Value",col="Category")) +
        ggplot2::xlab("") +
        ggplot2::geom_boxplot(lwd=2,fatten=1/4) +
        ggplot2::stat_summary(fun.y=mean,geom="point",shape=18,size=2) +
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
          legend.position="top",
          text=ggplot2::element_text(face="bold",color="black",size=15),
          panel.background=ggplot2::element_blank(),
          panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
          legend.title=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          strip.background=ggplot2::element_rect(fill="#0033cc"),
          strip.text = ggplot2::element_text(face="bold",color="white"),
          axis.text.x=ggplot2::element_blank())
      if(loglevel=="Fold"){ p1 <- p1 + ggplot2::ylab("Log-fold-change") } else p1 <- p1 + ggplot2::ylab("Log-read-count")
    } else {
      p1 <- ggplot2::ggplot(folded,ggplot2::aes_string("Value",color="Category")) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
          legend.position="top",
          legend.title=ggplot2::element_blank(),
          panel.background=ggplot2::element_blank(),
          panel.grid.major=ggplot2::element_blank(),
          panel.grid.minor=ggplot2::element_blank(),
          axis.line=ggplot2::element_line(colour="black",size=1),
          panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
          text=ggplot2::element_text(face="bold",color="black",size=15),
          axis.ticks=ggplot2::element_line(colour="black",size=1),
          strip.background=ggplot2::element_rect(fill="#0033cc"),
          strip.text = ggplot2::element_text(face="bold",color="white")
        ) +
        ggplot2::geom_density(ggplot2::aes_string(fill="Category"),stat="density",alpha=0.2,lwd=0) +
        ggplot2::geom_line(ggplot2::aes_string(fill="Category"),stat="density",size=3) +
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1,strip.position="top") + ggplot2::coord_flip() +
        ggplot2::scale_fill_discrete(guide=FALSE)
      if(loglevel=="Fold"){ p1 <- p1 + ggplot2::xlab("Log-fold-change") } else p1 <- p1 + ggplot2::xlab("Log-read-count")
    }
    return(p1)
  }
  
  plot.scores <- function(scored,anno,formula="~Condition+Replicate"){
    
    if(!"Category"%in%colnames(scored)) scored$Category <- "Other"
    my_data <- reshape2::melt(scored[scored$Category!="Nontargeting",],id.vars=c("Gene","Category"),value.name="Score",variable.name="Sample")
    my_data <- merge(my_data,anno)
    if(input$ea.style=="Density plot") {
      p1 <- ggplot2::ggplot(my_data[my_data$Category!="Other",],ggplot2::aes_string("Score",color="Category")) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
          legend.position="top",
          legend.title=ggplot2::element_blank(),
          panel.background=ggplot2::element_blank(),
          panel.grid.major=ggplot2::element_blank(),
          panel.grid.minor=ggplot2::element_blank(),
          axis.line=ggplot2::element_line(colour="black",size=1),
          panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
          text=ggplot2::element_text(face="bold",color="black",size=15),
          axis.ticks=ggplot2::element_line(colour="black",size=1),
          strip.background=ggplot2::element_rect(fill="#0033cc"),
          strip.text = ggplot2::element_text(face="bold",color="white")
        ) +
        ggplot2::geom_density(ggplot2::aes_string(fill="Category"),stat="density",alpha=0.2,lwd=0) +
        ggplot2::geom_line(ggplot2::aes_string(fill="Category"),stat="density",size=3) +
        ggplot2::geom_line(data=my_data[my_data$Category=="Other",],stat="density",linetype="dashed",size=2) +
        ggplot2::ylab("Density") +
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1,strip.position="top") + ggplot2::coord_flip() +
        ggplot2::scale_fill_discrete(guide=FALSE)
    } else {
      p1 <- ggplot2::ggplot(my_data,ggplot2::aes_string("Category","Score",col="Category")) +
        ggplot2::xlab("") +
        ggplot2::geom_boxplot(lwd=2,fatten=1/4) +
        ggplot2::stat_summary(fun.y=mean,geom="point",shape=18,size=2) +
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
          legend.position="top",
          text=ggplot2::element_text(face="bold",color="black",size=15),
          panel.background=ggplot2::element_blank(),
          panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
          legend.title=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          strip.background=ggplot2::element_rect(fill="#0033cc"),
          strip.text = ggplot2::element_text(face="bold",color="white"),
          axis.text.x=ggplot2::element_blank())
    }
    print(p1)
  }
  
  plot.reproducibility <- function(dataset, title) {
    if(!"Category"%in%colnames(dataset)) dataset$Category <- "All"
    validate(need(ncol(dataset)>2,"This condition has no replicates."))
    # if(ncol(dataset)<3) return(NULL)
    colors <- table(dataset$Category)
    colors <- sort(colors,decreasing=T)
    colors1 <- 1:length(colors)
    names(colors1) <- names(colors)
    others <- dataset$Category==names(colors1)[1]
    graphics::pairs(dataset[,-which(colnames(dataset)%in%c("sgRNA","Gene","Category")),drop=FALSE],main="",
                    lower.panel = function(x,y) {
                      graphics::par(usr=c(0,1,0,1))
                      graphics::text(0.5,0.5,round(stats::cor(x,y),digits=3),cex=2)
                    },
                    upper.panel = function(x,y) {
                      graphics::abline(a=0,b=1,col="gray",lwd=5)
                      graphics::points(x[others],y[others],pch=19,col=colors1[dataset$Category[others]])
                      graphics::points(x[!others],y[!others],pch=19,col=colors1[dataset$Category[!others]])
                    },
                    oma=c(3,3,9,3)
    )
    par(xpd = TRUE)
    graphics::legend("topright",col=1:length(colors),legend=names(colors),pch=16,bty='n',ncol=length(colors),pt.cex=3,cex=1.5)
  }
  
  
  
  ### INPUT ### -----------------------------------------------------------
  volumes <- shinyFiles::getVolumes()
  paths <- reactiveValues(fastq=character(0),out=character(0))
  
  # Output file:
  
  # A) in the "FASTQ Alignment" tab:
  ## working directory
  observeEvent(input$dir_out,{
    shinyFiles::shinyDirChoose(input,'dir_out',roots=volumes(),session=session,restrictions=system.file(package='base'))
    paths$out <- shinyFiles::parseDirPath(volumes,input$dir_out)
  })
  output$dir_out_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<b>Select working directory:</b>")))
    } else return(HTML(paste0("<b>Working directory:</b>")))
  })
  output$dir_out_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    return(paths$out)
  })
  
  ## FASTQ folder
  observeEvent(input$dir_fastq,{
    shinyFiles::shinyDirChoose(input,'dir_fastq',roots=volumes(),session=session,restrictions=system.file(package='base'))
    paths$fastq <- shinyFiles::parseDirPath(volumes,input$dir_fastq)
  })
  output$dir_fastq_text1 <- renderUI({
    if(identical(paths$fastq,character(0))) {
      return(HTML(paste0("<b>Select FASTQ folder:</b><br/>")))
    } else return(HTML(paste0("<b>FASTQ folder:</b><br/>")))
  })
  output$dir_fastq_text2 <- renderText({
    req(paths$fastq!=0)
    return(paths$fastq)
  })
  
  ##
  output$q1 <- downloadHandler(
    filename = "sample_counts.txt",
    content = function(file) {
      load("data/example/dang_cck81.rda")
      write.table(dang_cck81,file,quote=F,row.names=F,col.names=T,sep='\t')
    }
  )
  
  output$q2 <- downloadHandler(
    filename = "sample_library.txt",
    content = function(file) {
      load("data/example/dang_cck81_library.rda")
      write.table(dang_cck81_library,file,quote=F,row.names=F,col.names=T,sep='\t')
    }
  )
  
  output$q3 <- downloadHandler(
    filename = "sample_library.txt",
    content = function(file) {
      load("data/example/dang_cck81_library.rda")
      write.table(dang_cck81_library,file,quote=F,row.names=F,col.names=T,sep='\t')
    }
  )
  
  # B) in the "Quality control" tab:
  observeEvent(input$dir_qc,{
    shinyFiles::shinyDirChoose(input,'dir_qc',roots=volumes(),session=session,restrictions=system.file(package='base'))
    paths$out <- shinyFiles::parseDirPath(volumes,input$dir_qc)
  })
  
  
  output$counts_path=reactive({
    req(!identical(paths$out,character(0)))
    if(file.exists(file.path(paths$out,"counts.txt"))){
      return(file.path(paths$out,"counts.txt"))
    }else{return(ready$counts_path)}
  })
  output$library_path=reactive({
    req(!identical(paths$out,character(0)))
    if(file.exists(file.path(paths$out,"library.txt"))){
      return(file.path(paths$out,"library.txt"))
    }else{return(ready$library_path)}
  })

  output$counts_path_txt <- renderText({
    req(!identical(paths$out,character(0)))
    return(HTML(paste0("<b>Counts file: </b>",file.path(paths$out,"counts.txt"))))
  })
  
  output$library_path_txt <- renderText({
    req(!identical(paths$out,character(0)))
    return(HTML(paste0("<b>sgRNA library file: </b>",file.path(paths$out,"library.txt"))))
  })
  
  output$dir_qc_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<b>Select a working directory:</b><br/>")))
    } else return(HTML(paste0("<b>Working directory:</b><br/>")))
  })
  output$dir_qc_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    return(paths$out)
  })
  output$counts_qc <- renderUI({
    if(ready$fastq==0) {
      return(NULL)
    } else return(HTML(paste0("<b>Table of read counts:</b><br/>",ready$fastq)))
  })
  output$library_qc <- renderUI({
    if(ready$library==0) {
      return(NULL)
    } else return(HTML(paste0("<b>sgRNA library file:</b><br/>",ready$library)))
  })
  
  # C) in the "Essentiality analysis" tab:
  output$dir_ea_text1 <- renderUI({
    if(identical(paths$out,character(0)) || (!input$save_anno1)) {
      return(HTML(paste0("<font color=#FF0000>Please <b><u>Check</u></b> and <b><u>Anntate</u></b> samples first.</font><br/>")))
    } else return(HTML(paste0("<b>Working directory selected:</b><br/>")))
  })
  output$dir_ea_text2 <- renderText({
    req(!identical(paths$out,character(0)) && (input$save_anno1))
    return(paths$out)
  })
  
  # D) in the "Differential essentiality analysis" tab:
    observeEvent(input$dir_da,{
    shinyFiles::shinyDirChoose(input,'dir_da',roots=volumes(),session=session,restrictions=system.file(package='base'))
    paths$out <- shinyFiles::parseDirPath(volumes,input$dir_da)
  })
  
  output$dir_da_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<font><b>Select working folder:</b></font><br/>")))
    } else return(HTML(paste0("<b>working directory selected:</b><br/>")))
  })
  output$dir_da_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    
    dir_list <- list.dirs(paths$out,recursive = F)
    cond_ls = basename(dir_list[grep("\\.vs\\.", dir_list)])
    
    if(length(cond_ls) < 2){
      output$dir_da_err <- renderUI({return(HTML("<font color='#FF0000'><b>Error: Less than 2 analysis were found in this folder !!!</b></font>"))})
    }else{
      updateSelectInput(session, "condition1",
                      label = "Select condition comparison:",
                      choices = cond_ls,
      )
      updateSelectInput(session, "condition2",
                        label = "Select condition comparison:",
                        choices = cond_ls
      )
      output$dir_da_err <- renderUI({return(HTML("<font color='#0000FF'><b>Comparisons are detected, run DA.</b></font>"))})
    }
    return(paths$out)
  })
  
  
  # E) in the "Network analysis" tab:
  output$dir_net_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<font color=#FF0000>Please <b><u>Check</u></b> and <b><u>Anntate</u></b> samples first.</font><br/>")))
    } else return(HTML(paste0("<b>Working directory selected:</b><br/>")))
  })
  output$dir_net_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    diff_out <- basename(list.files(path=paths$out,pattern="_differential",full.names=TRUE))
    comparisons <- substring(diff_out,1,regexpr("_differential",diff_out)-1)
    updateSelectInput(session, "comparisons",
                      label = "Select condition comparison:",
                      choices = comparisons,
                      selected = comparisons[1]
    )
    return(paths$out)
  })

  ## 
  ### FASTQ ### -----------------------------------------------------------
  
  # Read counts from FASTQ files:
  counts <- eventReactive(input$button_run, {
    withProgress(message='Input',value=0,detail="Initializing...",{
      
      validate(need(paths$out!=0,"Please select a working directory."))
      validate(need(paths$fastq!=0,"Please select a directory containing FASTQ files.")) #temporary

      fastq <- read.fastq(id=paths$out,
                         fastq=paths$fastq,
                         library=input$file.library$datapath,
                         spacer_start=input$spacer.range[1]-1,
                         spacer_length=input$spacer.range[2]-(input$spacer.range[1]-1),
                         reverse=input$spacer.rev,
                         complement=input$spacer.comp,
                         shinyF=incProgress)  #temporary
      incProgress(0.3,"Input","Finalizing...")
      
      counts <- fastq$Counts
      annotation <- fastq$Annotation
      library <- fastq$Library
      counts_mat <- as.matrix(counts[,sapply(counts,is.numeric)])
      rownames(counts_mat) <- counts$sgRNA
    })
    showNotification("Please fill out the sample information in the 'File annotation' tab.",action=NULL,duration=20,closeButton=TRUE,type="message")
    
    return(list(Counts=counts,Library=library,Annotation=annotation,Counts_mat=counts_mat))
  })
  
  output$mapped2 <- renderPlot({
    req(counts()$Annotation)
    plot.mapped2(counts()$Annotation)
  })
  output$mapped1 <- renderPlot({
    req(counts()$Annotation)
    plot.mapped1(counts()$Annotation)
  })
  
  output$anno_table <- rhandsontable::renderRHandsontable({
    rtable <- counts()$Annotation
    rownames(rtable) <- NULL
    # rtable$Condition <- c(rep("Plasmid",4),rep("DMSO",12),rep("PRMT5",12)) #temporary
    # rtable$Replicate <- c(rep("A",4),rep(rep(c("A","B","C"),each=4),2)) #temporary
    rtable$Use <- "X"
    rhandsontable::rhandsontable(rtable) # converts the R dataframe to rhandsontable object
  })
  
  output$text_anno1 <- renderUI({
    if(!input$save_anno) {
      return(HTML(paste0("Please fill out the table below. Remove the 'X' on the FASTQ files you do not wish to use. Then click on the red button in order to merge the FASTQ files belonging to the same sample. <br/>")))
    } else if(!is.null(counts1()$Counts)) return(HTML(paste0("Merged table of counts saved at <br/>",ready$fastq)))
  })
  
  output$data_input <- DT::renderDataTable({ #DT is necessary to make the row names appear
    if(input$fastq.level=="FASTQ files") {
      req(counts())
      return(counts()$Counts)
    } else if(input$fastq.level=="Conditions") {
      req(counts1())
      return(counts1()$Counts)
    } else return(NULL)
  }, options = list(pageLength = 30), rownames=FALSE)
  
  
  ### COUNTS ### -----------------------------------------------------------
  
  counts1 <- eventReactive(c(input$save_anno,input$button_check),{
    # Either the FASTQ module must have been run, or the RUN button pressed:
    if(input$save_anno==0 & input$button_check==0) return(NULL)
    
    withProgress(message='Read counts',value=0,detail="Initializing...",{
      validate(need(paths$out!=0,"Please enter an output directory first."))

      if(identical(paths$fastq,character(0))) {
        merged <- read.counts(id=paths$out,
                                     library=input$file.library2$datapath,
                                     counts=input$file.counts$datapath,
                                     shinyF=incProgress)
      } else {
        validate(need(!is.null(input$anno_table),"Please run the FASTQ alignment module first."))
        rtable <- rhandsontable::hot_to_r(input$anno_table)
        rtable1 <- rtable[!is.na(rtable[,1]),]
        rtable1$Condition <- stringr::str_replace_all(rtable1$Condition,"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
        validate(need(sum(rtable1[,-ncol(rtable1)]=="")==0,"Please fill out the FASTQ file annotation first."))
        if(sum(rtable1[,-ncol(rtable1)]=="")>0)
          showNotification("Please fill out the FASTQ file annotation first.",action=NULL,duration=10,closeButton=TRUE,type="message")
        write.csv(rtable1[rtable1$Use=="X",-ncol(rtable1)],file=paste0(paths$out,"/fastq_annotation.csv"),row.names=FALSE,quote=FALSE)
        merged <- read.counts(id=paths$out,
                                     library=counts()$Library,
                                     counts=counts()$Counts[,colnames(counts()$Counts)%in%c("sgRNA","Gene","Category",rtable1$File[rtable1$Use=="X"]),drop=F],
                                     annotation=rtable1[rtable1$Use=="X",-ncol(rtable1)],
                                     shinyF=incProgress)
        ready$fastq <- paste0(paths$out,"/counts.txt")
        ready$library <- paste0(paths$out,"/library.txt")
      }
      
      counts <- merged$Counts
      annotation <- merged$Annotation
      library <- merged$Library
      counts_mat <- as.matrix(counts[,sapply(counts,is.numeric)])
      rownames(counts_mat) <- counts$sgRNA
    })
    return(list(Counts=counts,Library=library,Annotation=annotation,Counts_mat=counts_mat))
  })
  
  output$classification <- renderPlot({
    req(counts1())
    plot.classification(counts1()$Library,counts1()$Counts_mat)
  })
  output$sizes <- renderPlot({
    req(counts1())
    plot.sizes(counts1()$Library,counts1()$Counts_mat)
  })
  output$gini <- renderPlot({
    req(counts1()$Annotation)
    plot.gini(counts1()$Annotation)
  })
  output$depths <- renderPlot({
    req(counts1())
    plot.depths(counts1()$Counts)
  })
  output$zeros <- renderPlot({
    req(counts1())
    plot.zeros(counts1()$Counts)
  })
  
  ### Check and Annotation ### -----------------------------------------------
  
Annos = eventReactive(input$save_anno1, {
    withProgress(message='Quality Control',value=0,detail="Initializing...",{
      validate(need(paths$out!=0,"Please enter an output directory."))
      print("Saving annotation...")
      rtable <- rhandsontable::hot_to_r(input$anno1_table)
      rtable1 <- rtable[!is.na(rtable[,1]),]  
      write.csv(rtable1,file=paste0(paths$out,"/annotation.csv"),row.names=FALSE,quote=FALSE)
    
      annos <- sample.anno(id=paths$out,
                          counts=counts1()$Counts,
                          annotation=rtable1,
                          shinyF=incProgress)
      
      conditions <- unique(annos$Annotation$Condition)
      
      updateSelectInput(session, "samples.case",
                        choices = conditions
      )
      output$use_reps_case <- renderUI({
        req(input$samples.case)
        if(!input$use_reps){
          choices <- annos$Annotation$Sample[annos$Annotation$Condition %in% input$samples.case]
          selectizeInput("use_reps1",label="Choose case samples for analysis:",choices=choices,multiple=TRUE,selected=1,width="100%")
        }
      })
      
      
      updateSelectInput(session, "samples.ctrl",
                        choices = conditions
      )
      output$use_reps_ctrl <- renderUI({
        req(input$samples.ctrl)
        if(!input$use_reps_ctrl){
          choices <- annos$Annotation$Sample[annos$Annotation$Condition %in% input$samples.ctrl]
          selectizeInput("use_reps1_ctrl",label="Choose control samples for analysis:",choices=choices,multiple=TRUE,selected=1,width="100%")
        }
      })
    
    })
    
  
    print("Saving annotation...Done!")
    return(annos)
  })
  observe(Annos())
################################################################################
  output$text_anno2 <- renderUI({
    if(!input$save_anno1) {
      return(HTML(paste0("Please fill out the table and then click on the 'Save Annotation' when finished.<br/>")))
    } else return(HTML(paste0("The annotation data were updated and saved to ",paths$out,"/annotation.csv")))
  })
  
  output$anno1_table <- rhandsontable::renderRHandsontable({
    req(counts1()$Annotation)
    rtable <- counts1()$Annotation
    rownames(rtable) <- NULL
    # rtable$Condition <- c("CBP30","Day0A","DMSO","CBP30","Day0B","DMSO","CBP30","Day0C","DMSO") # temporary
    # rtable$Replicate <- rep(c("A","B","C"),each=3) # temporary
    # rtable$Control <- rep(c("Day0A","Day0B","Day0C"),each=3) # temporary
    
    # rtable$Condition <- c("Plasmid",rep("DMSO",3),rep("PRMT5",3))
    # rtable$Replicate <- c("R1",rep(c("R1","R2","R3"),2))
    # rtable$Control <- c("",rep("Plasmid",6))
    
    rhandsontable::rhandsontable(rtable) # converts the R dataframe to rhandsontable object
  })
  
  observeEvent(input$log.lfc,{
    if(input$log.lfc=="Log-fold-change"){
      conditions <- unique(qc()$QC$Annotation$Condition)
    } else conditions <- unique(qc()$QC$Annotation$Condition)
    updateSelectInput(session, "lfc.condition",
                      label = "Condition:",
                      choices = c("All conditions",conditions),
                      selected = "All conditions")
  })
  
  output$folds <- renderPlot({
    req(qc(),input$lfc.condition,input$folds.level,input$log.lfc)
    annotation2 <- qc()$QC$Annotation
    
    if(input$lfc.condition=="All conditions") {
      annotation2 <- annotation2[,c("Sample","Condition","Replicate")]
    } else annotation2 <- annotation2[annotation2$Condition==input$lfc.condition,c("Sample","Condition","Replicate")]
    
    # Dataset selection:
    loglevel = "Fold"
    if(input$folds.level=="sgRNA"){
      if(input$log.lfc=="Log-fold-change") {
        dataset <- qc()$QC$sgRNA_lfc_reps
      } else {
        dataset <- qc()$QC$sgRNA_log_reps
        loglevel = "Log"
      }
      if(!"Category"%in%colnames(dataset)) dataset$Category <- "All"
      dataset <- reshape2::melt(dataset,id.vars=c("sgRNA","Gene","Category"),value.name="Value",variable.name="Sample")
    } else {
      if(input$log.lfc=="Log-fold-change") {
        dataset <- qc()$QC$gene_lfc_reps
      } else {
        dataset <- qc()$QC$gene_log_reps
        loglevel = "Log"
      }
      if(!"Category"%in%colnames(dataset)) dataset$Category <- "All"
      dataset <- reshape2::melt(dataset,id.vars=c("Gene","Category"),value.name="Value",variable.name="Sample")
    }
    dataset <- merge(dataset,annotation2)
    plot.folds(dataset, loglevel=loglevel)
  }, height = function() { session$clientData$output_folds_width * 0.5 })
  
  observeEvent(input$corr_pca,{
    if(input$corr_pca=="Correlation") {
      updateRadioButtons(session,"corr_pca_method",label="Method:",choices=c("pearson","spearman"),selected="pearson")
    } else updateRadioButtons(session,"corr_pca_method",label="Components:",choices=c("2nd vs 1st","3rd vs 1st"),selected="2nd vs 1st")
  })
  
  output$correlation <- renderPlot({
    req(input$corr_pca=="Correlation")
    req(input$corr_pca_method%in%c("pearson","spearman"))
    validate(need(ncol(qc()$sgRNA_reps_mat)>1,"At least two samples are required in order to compute correlations."))
    if(input$corr.level=="sgRNA") {
      dataset <- qc()$sgRNA_reps_mat
    } else dataset <- qc()$gene_reps_mat
    plot.correlation(dataset,method=input$corr_pca_method)
  }, height = function() { session$clientData$output_correlation_width })
  
  output$pca1 <- renderPlot({
    req(input$corr_pca=="Principal component analysis")
    validate(need(ncol(qc()$sgRNA_reps_mat)>2,"More samples are required in order to compute PCA"))
    if(input$corr_pca_method=="2nd vs 1st")  {
      Y <- "Z2"
      Ylab <- "Second principal component"
    } else if(input$corr_pca_method=="3rd vs 1st")  {
      Y <- "Z3"
      Ylab <- "Third principal component"
    } else return(NULL)
    if(input$corr.level=="sgRNA") {
      dataset <- qc()$PCA1
    } else dataset <- qc()$PCA2
    ggplot2::ggplot(dataset, ggplot2::aes_string("Z1",Y,label="Sample")) +
      ggrepel::geom_text_repel() +
      ggplot2::geom_point() +
      ggplot2::ggtitle(input$corr_pca_method) +
      ggplot2::xlab("First principal component") +
      ggplot2::ylab(Ylab) +
      ggplot2::theme(
        title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
        legend.text = ggplot2::element_text(size=15,color="black"),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15),
        panel.background=ggplot2::element_blank(),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
        legend.title=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank(),
        strip.background=ggplot2::element_rect(fill="#0033cc"),
        strip.text = ggplot2::element_text(face="bold",color="white"),
        axis.text.x=ggplot2::element_blank())
  },height = function() { session$clientData$output_pca1_width })
  
  output$pve <- renderPlot({
    req(input$corr_pca=="Principal component analysis")
    validate(need(ncol(qc()$sgRNA_reps_mat)>2,"More samples are required in order to compute PCA"))
    if(input$corr.level=="sgRNA") {
      pve_sample <- qc()$PVE1
    } else pve_sample <- qc()$PVE2
    pve <- data.frame(PC=seq(1:length(pve_sample)),pve_sample=pve_sample)
    cpve <- data.frame(PC=seq(1:length(pve_sample)),Cpve_sample=cumsum(pve_sample))
    p1 <- ggplot2::ggplot(pve, ggplot2::aes(PC,pve_sample)) +
      ggplot2::geom_point(col="darkblue",shape="O",size=5) +
      ggplot2::geom_line(col="blue") +
      ggplot2::xlab("Principal Component") +
      ggplot2::ylab("Proportion of Variance Explained") +
      ggplot2::theme(
        title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
        legend.text = ggplot2::element_text(size=15,color="black"),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15),
        panel.background=ggplot2::element_blank(),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
        legend.title=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank(),
        strip.background=ggplot2::element_rect(fill="#0033cc"),
        strip.text = ggplot2::element_text(face="bold",color="white"),
        axis.text.x=ggplot2::element_blank())
    p2 <- ggplot2::ggplot(cpve, ggplot2::aes(PC,Cpve_sample)) +
      ggplot2::geom_line(col="blue") +
      ggplot2::geom_point(col="darkblue",shape="O",size=5) +
      ggplot2::xlab("Principal Component") +
      ggplot2::ylab("Cumulative PVE") +
      ggplot2::theme(
        title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
        legend.text = ggplot2::element_text(size=15,color="black"),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15),
        panel.background=ggplot2::element_blank(),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
        legend.title=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank(),
        strip.background=ggplot2::element_rect(fill="#0033cc"),
        strip.text = ggplot2::element_text(face="bold",color="white"),
        axis.text.x=ggplot2::element_blank())
    gridExtra::grid.arrange(p1,p2,nrow=1)
  })
  
  observeEvent(input$reproducibility.stage,{
    if(input$reproducibility.stage=="Log-fold-change"){
      conditions <- unique(qc()$QC$Annotation$Condition)
    } else conditions <- unique(qc()$QC$Annotation$Condition)
    updateSelectInput(session, "qc.condition",
                      label = "Visualize replicates of condition:",
                      choices = conditions,
                      selected = conditions[1])
  })
  
  output$reproducibility <- renderPlot({
    req(input$qc.condition)
    validate(need(ncol(qc()$sgRNA_reps_mat)>2,"More samples are required in order to show reproducibility."))
    if(input$reproducibility.stage=="Log-fold-change") {
      annotation2 <- qc()$QC$Annotation[qc()$QC$Annotation$Condition==input$qc.condition,]
      grp <- split(annotation2$Sample,annotation2[c("Condition")],drop=TRUE)
    } else {
      annotation2 <- qc()$QC$Annotation[qc()$QC$Annotation$Condition==input$qc.condition,]
      grp <- split(annotation2$Sample,annotation2[c("Condition")],drop=TRUE)
    }
    if(input$reproducibility.level=="sgRNA"){
      if(input$reproducibility.stage=="Log-fold-change") {
        cols <- which(colnames(qc()$QC$sgRNA_lfc_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$sgRNA_lfc_reps[,cols,drop=F]
      } else {
        cols <- which(colnames(qc()$QC$sgRNA_log_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$sgRNA_log_reps[,cols,drop=F]
      }
    } else {
      if(input$reproducibility.stage=="Log-fold-change") {
        cols <- which(colnames(qc()$QC$gene_lfc_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$gene_lfc_reps[,cols,drop=F]
      } else {
        cols <- which(colnames(qc()$QC$gene_log_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$gene_log_reps[,cols,drop=F]
      }
    }
    plot.reproducibility(dataset,names(grp))
  },height = function() { session$clientData$output_reproducibility_width })
  
  output$data_qc <- DT::renderDataTable({ #DT is necessary to make the row names appear
    req(qc())
    if(input$qc.level=="sgRNA") {
      if(input$qc.stage=="Log-fold-change") {
        if(input$qc.averaged==FALSE) return(qc()$QC$sgRNA_lfc_reps) else return(qc()$QC$sgRNA_lfc)
      } else {
        if(input$qc.averaged==FALSE) return(qc()$QC$sgRNA_log_reps) else return(qc()$QC$sgRNA_log)
      }
    } else {
      if(input$qc.stage=="Log-fold-change") {
        if(input$qc.averaged==FALSE) return(qc()$QC$gene_lfc_reps) else return(qc()$QC$gene_lfc)
      } else {
        if(input$qc.averaged==FALSE) return(qc()$QC$gene_log_reps) else return(qc()$QC$gene_log)
      }
    }
  }, options = list(pageLength = 30), rownames=FALSE)
  
  
  ### LFC and 2-Tail RRA ### -----------------------------------------------------------
  ranges <- reactiveValues(x=NULL,y=NULL,xd=NULL,yd=NULL,xr=NULL,yr=NULL) 
  qc <- eventReactive(input$button_ea, {
    print("Running analysis....")
    withProgress(message='Running LFC...',value=0,detail="Initializing...",{
      validate(need(paths$out!=0,"Please enter an output directory."))
      if(input$pseudo_default==TRUE) {
        pseudocount <- 0.1
      } else pseudocount <- input$pseudo.count
      
      counts = counts1()$Counts
      annotation = Annos()$Annotation
      # Check replicates:
      case_reps = input$use_reps1
      ctrl_reps = input$use_reps1_ctrl
      
      if(is.null(case_reps)) {
        case_reps =annotation$Sample[annotation$Condition %in% input$samples.case]
      } 
      
      if(is.null(ctrl_reps)) {
        ctrl_reps =annotation$Sample[annotation$Condition %in% input$samples.ctrl]
      }
      
      counts = counts[,c("sgRNA","Gene","Category",case_reps,ctrl_reps)]
      annotation = annotation[annotation$Sample %in% c(case_reps,ctrl_reps),]
      output_dir = file.path(paths$out,paste0(input$samples.case,".vs.",input$samples.ctrl))
      dir.exists(output_dir) || dir.create(output_dir, recursive = T)
      
      
      if(input$program_option_ea=='MAGeCKRRA'){
        incProgress(0.1,"Running MAGeCK test...","Running...")
        if(Sys.which("mageck")=="")
          stop("command 'mageck' is not found. You may check your MAGeCK installation.")
        
        ## MAGeCK RRA==========================================================
        mageck_params = list(norm.method = input$MAGeCKRRA.norm)
        
        qc <- mageck.test(
          id=output_dir,
          counts=file.path(path$out,"counts.txt"),
          annotation=file.path(path$out,"annotation.csv"),
          case_reps = case_reps,
          ctrl_reps = ctrl_reps,
          mageck_params = mageck_params,
          report=input$qc.report,
          shinyF=incProgress
        )
        incProgress(0.4,"Running MAGeCK test...","Finalizing...")
        return(list(QC=qc,sgRNA_reps_mat=sgRNA_reps_mat,gene_reps_mat=gene_reps_mat,PCA1=pca1,PVE1=pve1_sample,PCA2=pca2,PVE2=pve2_sample,rra=rra))
      }
      
      ##================================================================
      qc <- quality.control(id=output_dir,
                            counts=counts,
                            annotation=annotation,
                            case_reps = case_reps,
                            ctrl_reps = ctrl_reps,
                            pseudocount=pseudocount,
                            report=input$qc.report,
                            shinyF=incProgress)
      
      incProgress(0.7,"Running LFC...","Finalizing...")
      
      sgRNA_reps_mat <- as.matrix(qc$sgRNA_lfc_reps[,sapply(qc$sgRNA_lfc_reps,is.numeric),drop=FALSE])
      rownames(sgRNA_reps_mat) <- qc$sgRNA_lfc_reps$sgRNA
      if(ncol(sgRNA_reps_mat)>2){
        pca1_sample <- stats::prcomp(t(sgRNA_reps_mat),scale=F)
        pca1 <- data.frame(Z1=pca1_sample$x[,1],Z2=pca1_sample$x[,2],Z3=pca1_sample$x[,3],Sample=colnames(sgRNA_reps_mat))
        pve1_sample <- 100*pca1_sample$sdev^2/sum(pca1_sample$sdev^2)
      } else pca1 <- pve1_sample <- NULL
      
      gene_reps_mat <- as.matrix(qc$gene_lfc_reps[,sapply(qc$gene_lfc_reps,is.numeric),drop=FALSE])
      rownames(gene_reps_mat) <- qc$gene_lfc_reps$Gene
      if(ncol(sgRNA_reps_mat)>2){
        pca2_sample <- stats::prcomp(t(sgRNA_reps_mat),scale=F)
        pca2 <- data.frame(Z1=pca2_sample$x[,1],Z2=pca2_sample$x[,2],Z3=pca2_sample$x[,3],Sample=colnames(sgRNA_reps_mat))
        pve2_sample <- 100*pca2_sample$sdev^2/sum(pca2_sample$sdev^2)
      } else pca2 <- pve2_sample <- NULL
      
    print("Running MoPAC LFC....Done!")
   
    incProgress(0.5,"RRA","Running RRA ...")
  ## RRA --------------------
    print("Running MoPAC RRA....")
      validate(need(paths$out!=0,"Please enter an output directory."))
      
      rra <- RRA.2tail(id=output_dir,
                              replicates=case_reps,
                              pvalue=input$pvalue,
                              annotation=qc$Annotation,
                              shinyF=incProgress)
      
      incProgress(0.7,"Running LFC...","Finalizing...")
      
      conditions <- names(rra$p.value)
      
      updateSelectInput(session, "lfc.condition",
                        label = "Condition:",
                        choices =conditions,
                        selected = conditions[1]
      )
      updateSelectInput(session, "qc.condition",
                        label = "Visualize replicates of condition:",
                        choices = conditions,
                        selected = conditions[1]
      )
      
      updateSelectInput(session, "rra.condition",
                        label = "Visualize essentiality p values of condition:",
                        choices = conditions,
                        selected = conditions[1])
      
      updateSelectizeInput(session=session,inputId='rragene',choices=rra$p.value[[1]]$Gene,selected=NULL,server=TRUE)
      
      updateSelectInput(session, "ea.condition",
                        label = "Condition:",
                        choices = conditions,
                        selected = conditions[1]
      )
      updateSelectInput(session, "condition1",
                        label = "More essential in condition:",
                        choices = conditions,
                        selected = conditions[1]
      )
      updateSelectInput(session, "condition2",
                        label = "Less essential in condition:",
                        choices = conditions,
                        selected = conditions[1]
      )
      
    })
    print("Running MoPAC RRA....Done!")
    return(list(QC=qc,sgRNA_reps_mat=sgRNA_reps_mat,gene_reps_mat=gene_reps_mat,PCA1=pca1,PVE1=pve1_sample,PCA2=pca2,PVE2=pve2_sample,rra=rra))
  })
  observe(qc())
  
  plot.rra <- function() {
    req(qc())
    pvalue <- qc()$rra$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    if(!"Category"%in%colnames(pvalue)) pvalue$Category <- "All"
    largest <- table(pvalue$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)] <- "X"
    if(input$rracolor){
      p1 <- ggplot2::ggplot(pvalue,ggplot2::aes_string("X","p",col="Category")) +
        ggplot2::geom_point(data=pvalue[pvalue$Category==largest,],size=3,pch=1) +
        ggplot2::geom_point(data=pvalue[pvalue$Category!=largest,],size=3,pch=16)
    } else {
      p1 <- ggplot2::ggplot(pvalue,ggplot2::aes_string("X","p")) +
        ggplot2::geom_point(data=pvalue,size=3,pch=1,col="darkgray")
    }
    p1 <- p1 + ggplot2::geom_point(data=pvalue[pvalue$Gene%in%input$rragene,],size=4,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=pvalue[pvalue$Gene==flags$click1,],size=5,pch=21,col="black",fill="goldenrod") +
      ggplot2::ggtitle(paste0(input$rra.condition," (select area)")) +
      ggplot2::xlab("Log-Fold-Change") + ggplot2::ylab("-log10(P-value)") +
      # ggplot2::guides(col=ggplot2::guide_legend(nrow=2,byrow=FALSE)) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p1)
  }
  
  observeEvent(c(qc(),input$rra.condition,input$rragene,input$clickrrazoom$x,input$rracolor), {
    output$rra <- renderPlot({
      req(qc()$rra,input$rra.condition)
      isolate(plot.rra())
    },height = function() {session$clientData$output_rra_width })
  })
  
  # After hovering on the plot:
  output$plotrrainfo <- renderUI({
    req(qc()$rra,input$rra.condition)
    pvalue <- qc()$rra$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    hover <- input$rra_hover
    dataset <- pvalue
    point <- nearPoints(dataset,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)],yvar="p")
    if(nrow(point)==0) return(NULL)
    info <- dataset[dataset$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  # After clicking on the plot:
  flags <- reactiveValues(click1 = 0, click2 = 0) # to control the appearance/dissapearance of widgets
  
  observeEvent(input$clickrrazoom$x, {
    pvalue <- qc()$rra$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    if(!"Category"%in%colnames(pvalue)) pvalue$Category <- "All"
    point <- nearPoints(pvalue,input$clickrrazoom,threshold=10,maxpoints=1,xvar=colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)],yvar="p")
    flags$click1 <- point$Gene
    output$inforra <- renderUI({
      req(qc()$rra,flags$click1)
      gene <- flags$click1 #the sample on which the user clicked
      info <- paste0("<font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      info <- paste(info,paste0("Category: ",point$Category),sep="<br/>")
      for(condition in names(qc()$rra$p.value)){
        info <- paste(info,paste0("<br/><font color=#0000FF><b>Condition: ",condition,"</b></font>"),sep="<br/>")
        info <- paste(info,paste0("Fold-change = ",qc()$rra$p.value[[condition]][qc()$rra$p.value[[condition]]$Gene==flags$click1,3]),sep="<br/>")
        info <- paste(info,paste0("P-value = ",qc()$rra$p.value[[condition]][qc()$rra$p.value[[condition]]$Gene==flags$click1,"p"]),sep="<br/>")
        info <- paste(info,paste0("FDR = ",qc()$rra$p.value[[condition]][qc()$rra$p.value[[condition]]$Gene==flags$click1,"FDR"]),sep="<br/>")
        info <- paste(info,paste0("Depletion rank = ",qc()$rra$p.value[[condition]][qc()$rra$p.value[[condition]]$Gene==flags$click1,"depleted"]),sep="<br/>")
        info <- paste(info,paste0("Enrichment rank = ",qc()$rra$p.value[[condition]][qc()$rra$p.value[[condition]]$Gene==flags$click1,"enriched"]),sep="<br/>")
      }
      return(HTML(info))
    })
  })
  
  # Zoom in on Gene distributions:
  output$rrazoom <- renderPlot({
    req(qc()$rra,input$rra.condition)
    pvalue <- qc()$rra$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    if(!"Category"%in%colnames(pvalue)) pvalue$Category <- "All"
    largest <- table(pvalue$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)] <- "X"
    if(!is.null(input$plotrra_brush)) {
      brush <- input$plotrra_brush
      pvalue1 <- pvalue[pvalue$X>brush$xmin & pvalue$X<brush$xmax & pvalue$p>brush$ymin & pvalue$p<brush$ymax,]
      if(input$rracolor){
        p <- ggplot2::ggplot(pvalue1,ggplot2::aes_string("X","p",col="Category")) +
          ggplot2::geom_point(data=pvalue1[pvalue1$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=pvalue1[pvalue1$Category!=largest,],size=3,pch=16)
      } else {
        p <- ggplot2::ggplot(pvalue1,ggplot2::aes_string("X","p")) +
          ggplot2::geom_point(data=pvalue1,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=pvalue1[pvalue1$Gene%in%input$rragene,],size=4,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=pvalue1[pvalue1$Gene==flags$click1,],size=5,pch=21,col="black",fill="goldenrod")
    } else {
      if(input$rracolor){
        p <- ggplot2::ggplot(pvalue,ggplot2::aes_string("X","p",col="Category")) +
          ggplot2::geom_point(data=pvalue[pvalue$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=pvalue[pvalue$Category!=largest,],size=3,pch=16)
      } else {
        p <- ggplot2::ggplot(pvalue,ggplot2::aes_string("X","p")) +
          ggplot2::geom_point(data=pvalue,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=pvalue[pvalue$Gene%in%input$rragene,],size=4,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=pvalue[pvalue$Gene==flags$click1,],size=5,pch=21,col="black",fill="goldenrod")
    }
    p <- p + ggplot2::ggtitle(paste0(input$rra.condition," (click on a point)")) +
      ggplot2::xlab("Log-Fold-Change") + ggplot2::ylab("-log10(P-value)") +
      # ggplot2::guides(col=ggplot2::guide_legend(nrow=2,byrow=FALSE)) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p)
  }, height = function() { session$clientData$output_rrazoom_width })
  
  # After hovering on the plot:
  output$plotrrazoominfo <- renderUI({
    req(qc()$rra,input$rra.condition)
    pvalue <- qc()$rra$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    hover <- input$rrazoom_hover
    dataset <- pvalue
    point <- nearPoints(dataset,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=colnames(dataset)[which(colnames(pvalue)==input$rra.condition)],yvar="p")
    if(nrow(point)==0) return(NULL)
    info <- dataset[dataset$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  output$data_ea <- DT::renderDataTable({ #DT is necessary to make the row names appear
    req(qc()$rra,input$ea.condition)
    return(qc()$rra$p.value[[input$ea.condition]])
  }, options = list(pageLength = 30), rownames=FALSE)
  
  
  ### Differential essentiality ### ---------------------------------------------
  
  diff <- eventReactive(input$button_da,{
    withProgress(message='Differential analysis',value=0,detail="Initializing...",{
      diff <- analyze.differential(id=paths$out,
                                      condition1=input$condition1,
                                      condition2=input$condition2,
                                      report=input$da.report,
                                      read.unenriched=TRUE,
                                      shinyF=incProgress)
      incProgress(0.2,'Differential analysis',"Finalizing...")
      sgRNA <- diff$sgRNA
      Gene <- diff$Gene
      if(!"Category"%in%colnames(Gene))
        Gene <- data.frame(Gene=Gene$Gene,Category="All",Gene[,-1],check.names=FALSE,stringsAsFactors=FALSE)
      if(!"P(replicates)"%in%colnames(Gene))
        Gene$`P(replicates)` <- 2*pnorm(-abs(Gene$Z))
      # Gene$logP <- -log10(Gene$`P(sgRNAs)`)
    })
    diff_out <- basename(list.files(path=paths$out,pattern="_differential",full.names=TRUE))
    comparisons <- substring(diff_out,1,regexpr("_differential",diff_out)-1)
    updateSelectInput(session, "comparisons",
                      label = "Condition comparisons:",
                      choices = comparisons,
                      selected = comparisons[1]
    )
    updateSelectizeInput(session=session,inputId='dagene',choices=Gene$Gene,selected=NULL,server=TRUE)
    return(list(sgRNA=sgRNA,Gene=Gene))
  })
  
  observe(diff())
  
  # Plotting distributions:
  output$plot1 <- renderPlot({
    req(diff())
    dataset <- diff()$Gene
    largest <- table(dataset$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    if(input$gscatter){
      condition1 <- colnames(dataset)[5]
      condition2 <- colnames(dataset)[6]
    } else {
      condition1 <- colnames(dataset)[3]
      condition2 <- colnames(dataset)[4]
    }
    if(input$diffcolor){
      p <- ggplot2::ggplot(dataset,ggplot2::aes_string(condition1,condition2,col="Category")) +
        ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
        ggplot2::geom_point(data=dataset[dataset$Category==largest,],size=3,pch=1) +
        ggplot2::geom_point(data=dataset[dataset$Category!=largest,],size=3,pch=1)
    } else {
      p <- ggplot2::ggplot(dataset,ggplot2::aes_string(condition1,condition2)) +
        ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
        ggplot2::geom_point(data=dataset,size=3,pch=1,col="darkgray")
    }
    p <- p + ggplot2::geom_point(data=dataset[dataset$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=dataset[dataset$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod") +
      ggplot2::ggtitle("Essentiality (select area)") +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p)
  }, height = function() { session$clientData$output_plot1_width })
  
  # After hovering on the plot:
  output$plot1info <- renderUI({
    req(diff())
    if(input$gscatter){
      condition1 <- colnames(diff()$Gene)[5]
      condition2 <- colnames(diff()$Gene)[6]
    } else {
      condition1 <- colnames(diff()$Gene)[3]
      condition2 <- colnames(diff()$Gene)[4]
    }
    hover <- input$plot1_hover
    point <- nearPoints(diff()$Gene,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=condition1,yvar=condition2)
    if(nrow(point)==0) return(NULL)
    info <- diff()$Gene[diff()$Gene$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    # paste0("<b>Category: </b>",info$Category),
    # paste0("<b>Z score: </b>",info$Z),
    # paste0("<b>P value: </b>",info$P),
    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  
  observeEvent(input$click2$x, {
    req(diff())
    if(input$gscatter){
      condition1 <- colnames(diff()$Gene)[5]
      condition2 <- colnames(diff()$Gene)[6]
    } else {
      condition1 <- colnames(diff()$Gene)[3]
      condition2 <- colnames(diff()$Gene)[4]
    }
    point <- nearPoints(diff()$Gene,input$click2,threshold=5,maxpoints=1,addDist=TRUE,xvar=condition1,yvar=condition2)
    flags$click2 <- point$Gene
    output$infoscore <- renderUI({
      req(diff()$Gene,flags$click2)
      gene <- flags$click2 #the sample on which the user clicked
      info <- paste0("<br/><font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      info <- paste(info,paste0("Category: ",point$Category),sep="<br/>")
      info <- paste(info,paste0("Z score: ",formatC(point$Z,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("P value: ",formatC(point$P,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("FDR: ",formatC(point$FDR,format="e",digits=2)),sep="<br/>")
      return(HTML(info))
    })
  })
  
  # Zoom in on Gene distributions:
  output$plot2 <- renderPlot({
    req(diff())
    if(input$gscatter){
      condition1 <- colnames(diff()$Gene)[5]
      condition2 <- colnames(diff()$Gene)[6]
    } else {
      condition1 <- colnames(diff()$Gene)[3]
      condition2 <- colnames(diff()$Gene)[4]
    }
    dataset <- diff()$Gene
    largest <- table(dataset$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    if(!is.null(input$plot1_brush)) {
      brush <- input$plot1_brush
      dataset1 <- dataset[dataset[,condition1]>brush$xmin & dataset[,condition1]<brush$xmax & dataset[,condition2]>brush$ymin & dataset[,condition2]<brush$ymax,]
      if(input$diffcolor){
        p <- ggplot2::ggplot(dataset,ggplot2::aes_string(condition1,condition2,col="Category")) +
          ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
          ggplot2::geom_point(data=dataset1[dataset1$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=dataset1[dataset1$Category!=largest,],size=3,pch=1)
      } else {
        p <- ggplot2::ggplot(dataset,ggplot2::aes_string(condition1,condition2)) +
          ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
          ggplot2::geom_point(data=dataset1,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=dataset1[dataset1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=dataset1[dataset1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    } else {
      if(input$diffcolor){
        p <- ggplot2::ggplot(dataset,ggplot2::aes_string(condition1,condition2,col="Category")) +
          ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
          ggplot2::geom_point(data=dataset[dataset$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=dataset[dataset$Category!=largest,],size=3,pch=1)
      } else {
        p <- ggplot2::ggplot(dataset,ggplot2::aes_string(condition1,condition2)) +
          ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
          ggplot2::geom_point(data=dataset,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=dataset[dataset$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=dataset[dataset$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    }
    p <- p + ggplot2::ggtitle("Essentiality (click on a point)") +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p)
  }, height = function() { session$clientData$output_plot2_width })
  
  # After hovering on the plot:
  output$plot2info <- renderUI({
    req(diff())
    if(input$gscatter){
      condition1 <- colnames(diff()$Gene)[5]
      condition2 <- colnames(diff()$Gene)[6]
    } else {
      condition1 <- colnames(diff()$Gene)[3]
      condition2 <- colnames(diff()$Gene)[4]
    }
    hover <- input$plot2_hover
    point <- nearPoints(diff()$Gene,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=condition1,yvar=condition2)
    if(nrow(point)==0) return(NULL)
    info <- diff()$Gene[diff()$Gene$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    # paste0("<b>Category: </b>",info$Category),
    # paste0("<b>Z score: </b>",info$Z),
    # paste0("<b>P value: </b>",info$P),
    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  ### Volcano plots ###
  
  output$plot3 <- renderPlot({
    req(diff())
    tmp <- diff()$Gene
    if(input$p.volcano=="sgRNA reproducibility") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    largest <- table(tmp$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    if(input$diffcolor){
      p <- ggplot2::ggplot(tmp,ggplot2::aes(Z,logP,col=Category)) +
        ggplot2::geom_point(data=tmp[tmp$Category==largest,],size=3,pch=1) +
        ggplot2::geom_point(data=tmp[tmp$Category!=largest,],size=3,pch=1)
    } else {
      p <- ggplot2::ggplot(tmp,ggplot2::aes(Z,logP)) +
        ggplot2::geom_point(data=tmp,size=3,pch=1,col="darkgray")
    }
    p <- p + ggplot2::geom_point(data=tmp[tmp$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=tmp[tmp$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod") +
      ggplot2::ggtitle("Differential essentiality (select area)") +
      ggplot2::xlab("Z score") + ggplot2::ylab("-log10(P-value)") +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p)
  },height = function() {session$clientData$output_plot3_width })
  
  # After hovering on the plot:
  output$plot3info <- renderUI({
    req(diff())
    hover <- input$plot3_hover
    tmp <- diff()$Gene
    if(input$p.volcano=="sgRNA reproducibility") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    point <- nearPoints(tmp,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar="Z",yvar="logP")
    if(nrow(point)==0) return(NULL)
    info <- tmp[tmp$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    # paste0("<b>Category: </b>",info$Category),
    # paste0("<b>Z score: </b>",info$Z),
    # paste0("<b>P value: </b>",info$P),
    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  # Zoom in on Gene distributions:
  output$plot4 <- renderPlot({
    req(diff())
    tmp <- diff()$Gene
    if(input$p.volcano=="sgRNA reproducibility") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    largest <- table(tmp$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    if(!is.null(input$plot3_brush)) {
      brush <- input$plot3_brush
      dataset1 <- tmp[tmp$Z>brush$xmin & tmp$Z<brush$xmax & tmp$logP>brush$ymin & tmp$logP<brush$ymax,]
      if(input$diffcolor){
        p <- ggplot2::ggplot(tmp,ggplot2::aes(Z,logP,col=Category)) +
          ggplot2::geom_point(data=dataset1[dataset1$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=dataset1[dataset1$Category!=largest,],size=3,pch=1)
      } else {
        p <- ggplot2::ggplot(tmp,ggplot2::aes(Z,logP)) +
          ggplot2::geom_point(data=dataset1,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=dataset1[dataset1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=dataset1[dataset1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    } else {
      if(input$diffcolor){
        p <- ggplot2::ggplot(tmp,ggplot2::aes(Z,logP,col=Category)) +
          ggplot2::geom_point(data=tmp[tmp$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=tmp[tmp$Category!=largest,],size=3,pch=1)
      } else {
        p <- ggplot2::ggplot(tmp,ggplot2::aes(Z,logP)) +
          ggplot2::geom_point(data=tmp,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=tmp[tmp$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=tmp[tmp$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    }
    p <- p + ggplot2::ggtitle("Differential essentiality (click on a point)") +
      ggplot2::xlab("Z score") + ggplot2::ylab("-log10(P-value)") +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p)
  },height = function() { session$clientData$output_plot3_width })
  
  # After hovering on the plot:
  output$plot4info <- renderUI({
    req(diff())
    hover <- input$plot4_hover
    tmp <- diff()$Gene
    if(input$p.volcano=="sgRNA reproducibility") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    point <- nearPoints(tmp,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar="Z",yvar="logP")
    if(nrow(point)==0) return(NULL)
    info <- tmp[tmp$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    # paste0("<b>Category: </b>",info$Category),
    # paste0("<b>Z score: </b>",info$Z),
    # paste0("<b>P value: </b>",info$P),
    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  observeEvent(input$click4$x, {
    req(diff())
    tmp <- diff()$Gene
    if(input$p.volcano=="sgRNA reproducibility") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    point <- nearPoints(tmp,input$click4,threshold=5,maxpoints=1,addDist=TRUE,xvar="Z",yvar="logP")
    point1 <- tmp[tmp$Gene==point$Gene,]
    flags$click2 <- point$Gene
    output$infoscore <- renderUI({
      req(tmp,flags$click2)
      gene <- flags$click2 #the sample on which the user clicked
      info <- paste0("<br/><font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      info <- paste(info,paste0("Category: ",point1$Category),sep="<br/>")
      info <- paste(info,paste0("Z score: ",formatC(point1$Z,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>",colnames(tmp)[3]," score: ",formatC(as.numeric(point1[3]),format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("P(-): ",formatC(point1$`P(-)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("FDR(-): ",formatC(point1$`FDR(-)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>",colnames(tmp)[4]," score: ",formatC(as.numeric(point1[4]),format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("P(+): ",formatC(point1$`P(+)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("FDR(+): ",formatC(point1$`FDR(+)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>Replicate P value: ",formatC(point1$`P(replicates)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>sgRNA P value: ",formatC(point1$`P(sgRNAs)`,format="e",digits=2)),sep="<br/>")
      return(HTML(info))
    })
  })
  
  ### LOG-FOLD-CHANGES
  
  output$plot5 <- renderPlot({
    req(diff())
    cond_group1 = unlist(stringr::str_split(input$condition1, "\\.vs\\."))[1]
    cond_group2 = unlist(stringr::str_split(input$condition2, "\\.vs\\."))[1]
    if(input$sscatter){
      condition1 <- cond_group1
      condition2 <- cond_group2
    } else {
      condition1 <- paste0(cond_group1,"_biased")
      condition2 <- paste0(cond_group2,"_biased")
    }
    lfc1 <- diff()$sgRNA
    if(!"Category"%in%colnames(lfc1)) lfc1$Category <- "All"
    largest <- table(lfc1$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    QQ <- as.data.frame(stats::qqplot(x=lfc1[,condition1],y=lfc1[,condition2],plot.it=FALSE))
    if(input$diffcolor){
      p <- ggplot2::ggplot(lfc1,ggplot2::aes_string(condition1,condition2,col="Category")) +
        ggplot2::geom_point(data=lfc1[lfc1$Category==largest,],size=3,pch=1) +
        ggplot2::geom_point(data=lfc1[lfc1$Category!=largest,],size=3,pch=1)
    } else {
      p <- ggplot2::ggplot(lfc1,ggplot2::aes_string(condition1,condition2)) +
        ggplot2::geom_point(data=lfc1,size=3,pch=1,col="darkgray")
    }
    p <- p + ggplot2::geom_point(data=QQ,ggplot2::aes(x,y,col="QQ"),size=1) +
      ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
      ggplot2::geom_point(data=lfc1[lfc1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=lfc1[lfc1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod") +
      ggplot2::ggtitle("sgRNA log-fold-changes (select area)") +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p)
  }, height = function() { session$clientData$output_plot5_width })
  
  # After hovering on the plot:
  output$plot5info <- renderUI({
    req(diff())
    cond_group1 = unlist(stringr::str_split(input$condition1, "\\.vs\\."))[1]
    cond_group2 = unlist(stringr::str_split(input$condition2, "\\.vs\\."))[1]
    if(input$sscatter){
      condition1 <- cond_group1
      condition2 <- cond_group2
    } else {
      condition1 <- paste0(cond_group1,"_biased")
      condition2 <- paste0(cond_group2,"_biased")
    }
    hover <- input$plot5_hover
    point <- nearPoints(diff()$sgRNA,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=condition1,yvar=condition2)
    if(nrow(point)==0) return(NULL)
    info <- diff()$sgRNA[diff()$sgRNA$sgRNA==point$sgRNA,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),
                    paste0("<b>sgRNA: </b>",info$sgRNA),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  output$plot6 <- renderPlot({
    req(diff())
    cond_group1 = unlist(stringr::str_split(input$condition1, "\\.vs\\."))[1]
    cond_group2 = unlist(stringr::str_split(input$condition2, "\\.vs\\."))[1]
    if(input$sscatter){
      condition1 <- cond_group1
      condition2 <- cond_group2
    } else {
      condition1 <- paste0(cond_group1,"_biased")
      condition2 <- paste0(cond_group2,"_biased")
    }
    lfc1 <- diff()$sgRNA
    if(!"Category"%in%colnames(lfc1)) lfc1$Category <- "All"
    largest <- table(lfc1$Category)
    largest <- sort(largest,decreasing=T)
    largest <- names(largest)[1]
    QQ <- as.data.frame(stats::qqplot(x=lfc1[,condition1],y=lfc1[,condition2],plot.it=FALSE))
    if(!is.null(input$plot5_brush)) {
      brush <- input$plot5_brush
      dataset1 <- lfc1[lfc1[,condition1]>brush$xmin & lfc1[,condition1]<brush$xmax & lfc1[,condition2]>brush$ymin & lfc1[,condition2]<brush$ymax,]
      QQ1 <- QQ[QQ$x>brush$xmin & QQ$x<brush$xmax & QQ$y>brush$ymin & QQ$y<brush$ymax,]
      if(input$diffcolor){
        p <- ggplot2::ggplot(dataset1,ggplot2::aes_string(condition1,condition2,col="Category")) +
          ggplot2::geom_point(data=dataset1[dataset1$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=dataset1[dataset1$Category!=largest,],size=3,pch=1)
      } else {
        p <- ggplot2::ggplot(lfc1,ggplot2::aes_string(condition1,condition2)) +
          ggplot2::geom_point(data=dataset1,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=QQ1,ggplot2::aes(x,y,col="QQ"),size=1) +
        ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
        ggplot2::geom_point(data=dataset1[dataset1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=dataset1[dataset1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    } else {
      if(input$diffcolor){
        p <- ggplot2::ggplot(lfc1,ggplot2::aes_string(condition1,condition2,col="Category")) +
          ggplot2::geom_point(data=lfc1[lfc1$Category==largest,],size=3,pch=1) +
          ggplot2::geom_point(data=lfc1[lfc1$Category!=largest,],size=3,pch=1)
      } else {
        p <- ggplot2::ggplot(lfc1,ggplot2::aes_string(condition1,condition2,col="Category")) +
          ggplot2::geom_point(data=lfc1,size=3,pch=1,col="darkgray")
      }
      p <- p + ggplot2::geom_point(data=QQ,ggplot2::aes(x,y,col="QQ"),size=1) +
        ggplot2::geom_abline(ggplot2::aes(slope=1,col="X=Y",intercept=0)) +
        ggplot2::geom_point(data=lfc1[lfc1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=lfc1[lfc1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    }
    p <- p + ggplot2::ggtitle("sgRNA log-fold-changes (click on a point)") +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(size=15,color="black"),
                     axis.text = ggplot2::element_text(size=15,color="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    return(p)
  }, height = function() { session$clientData$output_plot6_width })
  
  # After hovering on the plot:
  output$plot6info <- renderUI({
    req(diff())
    cond_group1 = unlist(stringr::str_split(input$condition1, "\\.vs\\."))[1]
    cond_group2 = unlist(stringr::str_split(input$condition2, "\\.vs\\."))[1]
    if(input$sscatter){
      condition1 <- cond_group1
      condition2 <- cond_group2
    } else {
      condition1 <- paste0(cond_group1,"_biased")
      condition2 <- paste0(cond_group2,"_biased")
    }
    hover <- input$plot6_hover
    point <- nearPoints(diff()$sgRNA,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=condition1,yvar=condition2)
    if(nrow(point)==0) return(NULL)
    info <- diff()$sgRNA[diff()$sgRNA$sgRNA==point$sgRNA,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),
                    paste0("<b>sgRNA: </b>",info$sgRNA),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220
  
  observeEvent(input$click6$x, {
    req(diff())
    cond_group1 = unlist(stringr::str_split(input$condition1, "\\.vs\\."))[1]
    cond_group2 = unlist(stringr::str_split(input$condition2, "\\.vs\\."))[1]
    if(input$sscatter){
      condition1 <- cond_group1
      condition2 <- cond_group2
    } else {
      condition1 <- paste0(cond_group1,"_biased")
      condition2 <- paste0(cond_group2,"_biased")
    }
    point <- nearPoints(diff()$sgRNA,input$click6,threshold=5,maxpoints=1,addDist=TRUE,xvar=condition1,yvar=condition2)
    flags$click2 <- point$Gene
    point1 <- diff()$Gene[diff()$Gene$Gene==point$Gene,]
    output$infoscore <- renderUI({
      req(diff(),flags$click2)
      gene <- flags$click2 #the sample on which the user clicked
      info <- paste0("<br/><font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      if("Category"%in%colnames(diff()$Gene)) info <- paste(info,paste0("Category: ",point1$Category),sep="<br/>")
      info <- paste(info,paste0("Z score: ",formatC(point1$Z,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("Replicate P value: ",formatC(point1$`P(replicates)`,format="e",digits=2)),sep="<br/>")
      return(HTML(info))
    })
  })
  
  ### TABLE OF SCORES ###
  
  output$table_da <- DT::renderDataTable({ #DT is necessary to make the row names appear
    req(diff())
    if(input$data_da=="Gene"){
      return(
        DT::datatable(diff()$Gene) %>% 
          DT::formatRound(purrr::map_lgl(.$x$data, is.numeric), digits = 4))
    } else return(DT::datatable(diff()$sgRNA) %>% 
                    DT::formatRound(purrr::map_lgl(.$x$data, is.numeric), digits = 4))
  }, options = list(pageLength = 30), rownames=FALSE) 
  
  
  
}

