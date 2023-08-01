

sample.anno <- function(id=NULL,counts=NULL,annotation=NULL,shinyF=NULL) {
  #######################################################################################################################################
  ########################### INPUT #####################################################################################################
  #######################################################################################################################################
  if(!is.null(id)) { if(!dir.exists(id)) stop("The working directory selected was not found") }
  # Read counts -------------------------------------------------------------------------
  print("Importing data..."); if(!is.null(shinyF)) shinyF(0.1,"Check and annotate","Importing data...")
  if(is.null(id) & is.null(counts)) { # insufficient information provided
    stop("Please specify a working directory or provide a table of read counts.")
  } else if(!is.null(id) & is.null(counts)) { # is a table of read counts available in the working directory?
    if(file.exists(paste0(id,"/counts.txt"))) {
      counts0 <- utils::read.table(paste0(id,"/counts.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Please provide a table of read counts or run the FASTQ alignment module first.")
  } else if(is.data.frame(counts)) { # was a data frame provided as a table of read counts?
    counts0 <- counts
    counts0[,sapply(counts0,is.factor)] <- sapply(counts0[,sapply(counts0,is.factor)],as.character)
  } else if(is.character(counts)) { # was a file path provided as a table of read counts?
    if(file.exists(counts)) {
      if(tools::file_ext(counts)=="xlsx") {
        counts0 <- openxlsx::read.xlsx(counts)
      } else if(tools::file_ext(counts)=="csv") {
        counts0 <- utils::read.csv(counts,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"")
      } else counts0 <- utils::read.table(counts,header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Counts file not found in the specified address.")
  } else stop("Please provide a valid table of read counts.")
  colnames(counts0)[grep("sgRNA|sgrna|SGRNA",colnames(counts0))] <- "sgRNA"
  colnames(counts0)[grep("Gene|gene|GENE",colnames(counts0))] <- "Gene"
  colnames(counts0)[grep("Category|category|CATEGORY",colnames(counts0))] <- "Category"
  if(!"sgRNA"%in%colnames(counts0) | !"Gene"%in%colnames(counts0))
    stop("Please run the input module first.(sgRNA or Gene column is missing!)")
  # Read annotation -------------------------------------------------------------------------
  if(is.null(id) & is.null(annotation)) { # insufficient information provided
    stop("Please specify a working directory or provide a data frame of annotation.")
  } else if(!is.null(id) & is.null(annotation)) { # is a table of annotation available in the working directory?
    if(file.exists(paste0(id,"/annotation.csv"))) {
      annotation0 <- utils::read.csv(paste0(id,"/annotation.csv"),header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
      annotation0[annotation0 == ""] <- NA   
      if(sum(is.na(annotation0))>0)
        annotation0 <- openxlsx::read.xlsx(paste0(id,"/annotation.xlsx"),check.names=FALSE,na.strings="")
      annotation0[annotation0 == ""] <- NA   
      if(sum(is.na(annotation0))>0)
        stop("Please fill out the sample file annotation before proceeding.")
    } else stop("Please provide a table of annotation or run the input module first.")
  } else if(is.data.frame(annotation)) { # was a data frame provided as a table of read counts?
    annotation0 <- annotation
    annotation0[annotation0 == ""] <- NA     
    if(sum(is.na(annotation0))>0)
      stop("Please fill out the sample file annotation before proceeding.")
    annotation0[,sapply(annotation0,is.factor)] <- sapply(annotation0[,sapply(annotation0,is.factor)],as.character)
  }  else if(is.character(annotation)) {
    if(file.exists(annotation)) {
      if(tools::file_ext(annotation)=="xlsx") {
        annotation0 <- openxlsx::read.xlsx(annotation,check.names=FALSE,na.strings="")
      } else if(tools::file_ext(annotation)=="csv") {
        annotation0 <- utils::read.csv(annotation,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
      } else annotation0 <- utils::read.table(annotation,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
    } else stop("Annotation file not found in the specified address.")
    
    annotation0[annotation0 == ""] <- NA   
    if(sum(is.na(annotation0))>0)
      stop("Please fill out the sample file annotation before proceeding.")
    } else stop("The annotation object provided could not be read.")
  if(sum(!c("Sample","Condition","Replicate")%in%colnames(annotation0))>0)
    stop("Please make sure the annotation contains the columns: Sample, Condition, Replicate")
  
  print("Done.")
  return(list(Counts=counts0, Annotation=annotation0))
}
