

quality.control <- function(id=NULL,counts=NULL,annotation=NULL,case_reps=NULL,ctrl_reps=NULL, pseudocount=0.1,report=FALSE,shinyF=NULL,mageck_params=NULL) {
  # id="~/test/test_mar26";counts=NULL;annotation=NULL;pseudocount=0.1;report=FALSE;shinyF=NULL
  # id=NULL;counts=NULL;annotation=NULL;pseudocount=0.1;report=FALSE;shinyF=NULL
  # counts=counts$Counts; annotation=annotation; report=FALSE
  #######################################################################################################################################
  ########################### INPUT #####################################################################################################
  #######################################################################################################################################
  if(!is.null(id)) { if(!dir.exists(id)) stop("The working directory selected was not found") }
  # Read counts -------------------------------------------------------------------------
  print("Importing data..."); if(!is.null(shinyF)) shinyF(0.1,"Quality Control","Importing data...")
  if(is.null(id) & is.null(counts)) { # insufficient information provided
    stop("Please specify a working directory or provide a table of read counts.")
  } else if(!is.null(id) & is.null(counts)) { # is a table of read counts available in the working directory?
    if(file.exists(paste0(id,"/counts.txt"))) {
      counts = paste0(id,"/counts.txt")
    } else stop("Please provide a table of read counts or run the FASTQ alignment module first.")
  } else if(is.character(counts)) { # was a file path provided as a table of read counts?
    if(!file.exists(counts)) {
       stop(paste0("Counts file not found at: ",counts))
    }
  } else stop("Please provide a valid table of read counts.")
  
  # Read annotation -------------------------------------------------------------------------
  if(is.null(id) & is.null(annotation)) { # insufficient information provided
    stop("Please specify a working directory or provide a data frame of annotation.")
  } else if(!is.null(id) & is.null(annotation)) { # is a table of annotation available in the working directory?
    if(file.exists(paste0(id,"/annotation.csv"))) {
      annotation = paste0(id,"/annotation.csv")
    } else stop("Please provide a table of sample annotations.")
  } else if(is.character(annotation)) {
    if(!file.exists(annotation)) {  stop(paste0("Counts file not found at: ",annotation))}
  } else stop("The annotation object provided could not be read.")
  
  #######################################################################################################################################
  ########################### PROCESSING ################################################################################################
  #######################################################################################################################################
  print("Processing data...")
  if(!is.null(shinyF))
    shinyF(0.2,"MAGeCK...","Processing data...")
  
  mageck_cmd = paste("mageck test -k",counts,"-t",case_reps,"-c",ctrl_reps,"-n",paste0(id,"/rra"),"--adjust-method",mageck_params$norm.method,"--normcounts-to-file")
  if(report==TRUE) mageck_cmd = paste(mageck_cmd,"--pdf-report")
  system(mageck_cmd)
  
  norm_counts <- utils::read.table(paste0(id,"/rra.normalized.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
  
  #############
  annotation0 <- utils::read.csv(annotation,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
  groups <- split(annotation0[,"Sample",drop=F],annotation0$Condition)
  ctrl_group = unique(annotation0$Condition[annotation0$Sample %in% ctrl_reps])
  annotation1 = annotation0[!(annotation0$Condition %in% ctrl_group),]
  
  sgRNA_log_reps = norm_counts
  
  control_log_reps = as.data.frame(unname(lapply(groups,function(y)sgRNA_log_reps[,intersect(y$Sample,ctrl_reps),drop=T])))
  control_log_means = rowMeans(control_log_reps)
  
  sgRNA_lfc_reps <- cbind(sgRNA_log_reps[,which(colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category")),drop=F],do.call(cbind,unname(lapply(groups,function(y)
    sgRNA_log_reps[,y$Sample,drop=T] - control_log_means))))
  sgRNA_lfc_reps = sgRNA_lfc_reps[,!(names(sgRNA_lfc_reps) %in% ctrl_reps)]
  
  # Gene level -------------------------------------------------------------------------------------------------
  # Remove genes with a single sgRNA:
  NsgRNA <- table(sgRNA_log_reps$Gene)
  gene_log_reps <- sgRNA_log_reps[sgRNA_log_reps$Gene%in%names(NsgRNA)[NsgRNA>1],]
  # Average sgRNA log2 reads per gene:
  gene_log_reps <- do.call(rbind,lapply(split(sgRNA_log_reps[,-which(colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category")),drop=FALSE],sgRNA_log_reps$Gene),colMeans))
  gene_log_reps <- data.frame(Gene=rownames(gene_log_reps),gene_log_reps,check.names=FALSE,stringsAsFactors=FALSE)
  if("Category"%in%colnames(sgRNA_log_reps))
    gene_log_reps <- merge(unique(sgRNA_log_reps[,which(colnames(sgRNA_log_reps)%in%c("Gene","Category")),drop=FALSE]),gene_log_reps)
  
  # Gene log-fold-change:
  control_log_reps_gene = as.data.frame(unname(lapply(groups,function(y) gene_log_reps[,intersect(y$Sample,ctrl_reps),drop=T])))
  control_log_means_gene = rowMeans(control_log_reps_gene)
  
  gene_lfc_reps <- cbind(gene_log_reps[,which(colnames(gene_log_reps)%in%c("Gene","Category")),drop=F],do.call(cbind,unname(lapply(groups,function(y)
    gene_log_reps[,y$Sample,drop=T] - control_log_means_gene))))
  gene_lfc_reps = gene_lfc_reps[,!(names(gene_lfc_reps) %in% ctrl_reps)]
  
  # Replicate average -------------------------------------------------------------------------------------------
  # Log reads:
  groups <- split(as.character(annotation0$Sample),annotation0[,c("Condition"),drop=F],drop=T)
  sgRNA_log <- do.call(cbind,lapply(groups,function(y)rowMeans(sgRNA_log_reps[,y,drop=F])))
  sgRNA_log <- cbind(sgRNA_log_reps[,which(colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category")),drop=F],sgRNA_log)
  gene_log <- do.call(cbind,lapply(groups,function(y)rowMeans(gene_log_reps[,y,drop=F])))
  gene_log <- cbind(gene_log_reps[,which(colnames(gene_log_reps)%in%c("Gene","Category")),drop=FALSE],gene_log)
  
  # Fold changes:
  groups <- split(as.character(annotation1$Sample),annotation1[,c("Condition"),drop=F],drop=T)
  sgRNA_lfc <- do.call(cbind,lapply(groups,function(y)rowMeans(sgRNA_lfc_reps[,y,drop=F])))
  sgRNA_lfc <- cbind(sgRNA_lfc_reps[,which(colnames(sgRNA_lfc_reps)%in%c("sgRNA","Gene","Category")),drop=F],sgRNA_lfc)
  gene_lfc <- do.call(cbind,lapply(groups,function(y)rowMeans(gene_lfc_reps[,y,drop=F])))
  gene_lfc <- cbind(gene_lfc_reps[,which(colnames(gene_lfc_reps)%in%c("Gene","Category")),drop=F],gene_lfc)
  #######################################################################################################################################
  ########################### OUTPUT ####################################################################################################
  #######################################################################################################################################
  rowSort <- order(sgRNA_log_reps$Gene,sgRNA_log_reps$sgRNA) #careful: some sgRNAs are repeated for multiple genes!
  infoCols <- colnames(sgRNA_log_reps)[colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category")]
  sgRNA_log_reps <- sgRNA_log_reps[rowSort,c(infoCols,annotation0$Sample[order(annotation0$Condition)])]
  sgRNA_lfc_reps <- sgRNA_lfc_reps[rowSort,c(infoCols,annotation1$Sample[order(annotation1$Condition)])]
  
  sgRNA_log <- sgRNA_log[rowSort,c(infoCols,gtools::mixedsort(colnames(sgRNA_log)[-which(colnames(sgRNA_log)%in%c("sgRNA","Gene","Category"))]))]
  sgRNA_lfc <- sgRNA_lfc[rowSort,c(infoCols,gtools::mixedsort(colnames(sgRNA_lfc)[-which(colnames(sgRNA_lfc)%in%c("sgRNA","Gene","Category"))]))]
  
  infoCols <- colnames(gene_lfc_reps)[colnames(gene_lfc_reps)%in%c("Gene","Category")]
  gene_log_reps <- gene_log_reps[,c(infoCols,gtools::mixedsort(colnames(gene_log_reps)[-which(colnames(gene_log_reps)%in%c("Gene","Category"))]))]
  gene_lfc_reps <- gene_lfc_reps[,c(infoCols,gtools::mixedsort(colnames(gene_lfc_reps)[-which(colnames(gene_lfc_reps)%in%c("Gene","Category"))]))]
  
  gene_log <- gene_log[,c(infoCols,gtools::mixedsort(colnames(gene_log)[-which(colnames(gene_log)%in%c("Gene","Category"))]))]
  gene_lfc <- gene_lfc[,c(infoCols,gtools::mixedsort(colnames(gene_lfc)[-which(colnames(gene_lfc)%in%c("Gene","Category"))]))]
  
  if(!is.null(id)) {
    print("Generating output...")
    if(!is.null(shinyF))
      shinyF(0.3,"Quality Control","Generating output...")
    utils::write.csv(annotation0,paste0(id,"/annotation.csv"),quote=FALSE,row.names=FALSE)
    utils::write.table(counts0,paste0(id,"/counts.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_log_reps,paste0(id,"/sgRNA_log_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_lfc_reps,paste0(id,"/sgRNA_lfc_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_log,paste0(id,"/sgRNA_log.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_lfc,paste0(id,"/sgRNA_lfc.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_log_reps,paste0(id,"/gene_log_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_lfc_reps,paste0(id,"/gene_lfc_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_log,paste0(id,"/gene_log.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_lfc,paste0(id,"/gene_lfc.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    # Generate PDF file of quality control --------------------------------------------------------------------
  }
  print("Done.")
  return(list(Counts=counts0, Annotation=annotation1,
              sgRNA_log_reps=sgRNA_log_reps, sgRNA_lfc_reps=sgRNA_lfc_reps, sgRNA_log=sgRNA_log, sgRNA_lfc=sgRNA_lfc,
              gene_log_reps=gene_log_reps, gene_log=gene_log, gene_lfc_reps=gene_lfc_reps, gene_lfc=gene_lfc))
}
