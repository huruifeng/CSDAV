#' @title Gene essentiality 2-tail robust rank aggregation
#'
#' @description Performs 2-tail a-RRA analysis of gene essentiality. Within MoPAC's framework, this is necessary if there are no control genes available for the normalization of two conditions when computing the differential essentiality scores, or if the control genes should be filtered when computing the essentiality scores.
#'
#' @param id an optional identifier for the analysis. All input will be read from and all output will be stored in a directory with this name. If \code{NULL}, two data frames must be supplied containing the sgRNA-level and gene-level replicate-averaged log-fold-changes. Defaults to \code{NULL}.
#' @param sgRNA_lfc_reps an optional data frame containing the replicate-averaged log-fold-changes for each sgRNA. If \code{NULL}, the corresponding file in the identifier directory will be read.
#' @param gene_lfc_reps an optional data frame containing the replicate-averaged log-fold-changes for each gene If \code{NULL}, the corresponding file in the identifier directory will be read.
#' @param annotation an optional data frame or file path containing the sample annotation information. If \code{NULL}, the annotation file in the id directory will be read. Defaults to \code{NULL}.
#' @param replicates the replicates selected to continue the analysis
#' @param pvalue numeric. P-value threshold for enrichment. Defaults to 0.05.
#' @param shinyF logical. An internal parameter used by the graphical user interface. If using command-line leave as \code{NULL}.
#'
#' @details
#' A modified version of the alpha-robust-rank-aggregation algorithm from MAGeCK-RRA is implemented whereby the skewness of the sgRNA distributions is simultaneously assessed for both negative and positive selection.
#' A gene with a significant p value may thus be either significantly depleted or significantly enriched in the screen, and it is therefore necessary to consider the sign of its log-fold-change as well in order to distinguish between negative and positive selection.
#'
#' If \code{id} is not \code{NULL}, the following files are produced in the directory with name equal to the identifier:
#' \itemize{
#'  \item{\code{RRA.txt} files, one for each condition.}
#'  \item{\code{unenriched.txt}, the list of unenriched genes.}
#' }
#'
#' @return A list containing the following components:
#' \item{p.value}{List of tables of p values and FDR for all genes of each condition, indicating the ranking of depletion/enrichment.}
#' \item{Unenriched}{List of genes which were not found to be significantly enriched in any of the conditions analyzed.}
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords RRA essentiality
#'
#' @references Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol. 15, 554 (2014).
#'
#' @examples
#'
#' ## Example I: Two-tail RRA using an identifier to produce output:
#'
#' ## 1. Prepare read counts information:
#' counts <- read.counts(id="DANG_CCK81", counts=MoPAC::dang_cck81, library=MoPAC::sgRNA_library)
#'
#' ## 2. Manually fill out either of the two annotation files produced by read.counts.
#' ## In this example we will fill it out using R:
#' annotation <- counts$Annotation
#' annotation$Condition <- c("Plasmid0","Plasmid1",rep("DANG",4),rep("CCK81",3))
#' annotation$Replicate <- c("A","A","A","B","C","D","A","B","C")
#' annotation$Control <- c("","",rep("Plasmid0",4),rep("Plasmid1",3))
#' utils::write.csv(annotation,"DANG_CCK81/annotation.csv",quote=FALSE,row.names=FALSE)
#'
#' ## 3. Get pre-processed fold changes and Finally get the two-tail RRA essentiality analysis:
#' significance <- RRA.2tail(id="DANG_CCK81", annotation=annotation)
#'
#' ## Example II: Two-tail RRA without using an identifier:
#'
#' ## 1. Prepare read counts information:
#' counts <- read.counts(counts=MoPAC::dang_cck81, library=MoPAC::sgRNA_library)
#'
#' ## 2. Prepare sample annotation as a data frame:
#' annotation <- counts$Annotation
#' annotation$Condition <- c("Plasmid0","Plasmid1",rep("DANG",4),rep("CCK81",3))
#' annotation$Replicate <- c("A","A","A","B","C","D","A","B","C")
#' annotation$Control <- c("","",rep("Plasmid0",4),rep("Plasmid1",3))
#'
#' ## 3. Get pre-processed fold changes:
#' qc <- quality.control(counts=counts$Counts, annotation=annotation)
#'
#' ## 4. Finally get the two-tail RRA essentiality analysis:
#' significance <- RRA.2tail(sgRNA_lfc_reps=qc$sgRNA_lfc_reps,gene_lfc_reps=qc$gene_lfc_reps,annotation=annotation)
#'
#' ##--------------------------------------------------------------------------------------
#'
#' @export
#'
library(Rcpp)
sourceCpp("src/rra2tail.cpp")

RRA.2tail <- function(id=NULL,sgRNA_lfc_reps=NULL,gene_lfc_reps=NULL,annotation=NULL,replicates=NULL,pvalue=0.05,shinyF=NULL) {
  # id=NULL;sgRNA_lfc_reps=NULL;gene_lfc_reps=NULL;annotation=NULL;replicates=NULL;pvalue=0.05;shinyF=NULL
  # sgRNA_lfc_reps=qc$sgRNA_lfc_reps;gene_lfc_reps=qc$gene_lfc_reps;annotation=annotation
  set.seed(123)
  #######################################################################################################################################
  ########################### INPUT #####################################################################################################
  #######################################################################################################################################
  print("Importing data...")
  # Import sgRNA-level log-fold-changes ----------------------------------------------------
  if(is.null(id) & is.null(sgRNA_lfc_reps)) { # insufficient information provided
    stop("Please specify a working directory or provide the sgRNA log-fold-changes per replicate.")
  } else if(!is.null(id) & is.null(sgRNA_lfc_reps)) { # is a table of read sgRNA_reps available in the working directory?
    if(file.exists(paste0(id,"/sgRNA_lfc_reps.txt"))) {
      sgRNA_lfc_reps <- utils::read.table(paste0(id,"/sgRNA_lfc_reps.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Please provide the sgRNA log-fold-changes generated by the quality control module.")
  } else if(is.data.frame(sgRNA_lfc_reps)) {
    sgRNA_lfc_reps[,sapply(sgRNA_lfc_reps,is.factor)] <- sapply(sgRNA_lfc_reps[,sapply(sgRNA_lfc_reps,is.factor)],as.character)
  } else if(is.character(sgRNA_lfc_reps)) { # was a file path provided?
    if(file.exists(sgRNA_lfc_reps)) {
      sgRNA_lfc_reps <- utils::read.table(sgRNA_lfc_reps,header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("sgRNA log-fold-change file not found in the specified address.")
  } else stop("Please provide a valid table of sgRNA log-fold-changes.")
  if(!"sgRNA"%in%colnames(sgRNA_lfc_reps) | !"Gene"%in%colnames(sgRNA_lfc_reps))
    stop("Please run the quality control module first.")
  # Import Gene-level log-fold-changes ----------------------------------------------------
  if(is.null(id) & is.null(gene_lfc_reps)) { # insufficient information provided
    stop("Please specify a working directory or provide the gene log-fold-changes.")
  } else if(!is.null(id) & is.null(gene_lfc_reps)) { # is a table of read gene_reps available in the working directory?
    if(file.exists(paste0(id,"/gene_lfc_reps.txt"))) {
      gene_lfc_reps <- utils::read.table(paste0(id,"/gene_lfc_reps.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Please provide the gene log-fold-changes generated by the quality control module.")
  } else if(is.data.frame(gene_lfc_reps)) {
    gene_lfc_reps[,sapply(gene_lfc_reps,is.factor)] <- sapply(gene_lfc_reps[,sapply(gene_lfc_reps,is.factor)],as.character)
  } else if(is.character(gene_lfc_reps)) { # was a file path provided?
    if(file.exists(gene_lfc_reps)) {
      gene_lfc_reps <- utils::read.table(gene_lfc_reps,header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Gene log-fold-change file not found in the specified address.")
  } else stop("Please provide a valid table of gene log-fold-changes.")
  if(!"Gene"%in%colnames(gene_lfc_reps))
    stop("Please run the quality control module first.")
  # Check replicates:
  if(!is.null(replicates)) {
    if(sum(replicates%in%colnames(sgRNA_lfc_reps))<length(replicates))
      stop("Selected replicates not found in the column names.")
  } else replicates <- colnames(sgRNA_lfc_reps)[-which(colnames(sgRNA_lfc_reps)%in%c("sgRNA","Gene","Category"))]
  # Replicate average -------------------------------------------------------------------------------------------
  annotation <- annotation[annotation$Sample%in%replicates,]
  gene_lfc_reps <- gene_lfc_reps[,which(colnames(gene_lfc_reps)%in%c("Gene","Category",replicates)),drop=F]
  groups <- split(as.character(annotation$Sample),annotation[,c("Condition"),drop=F],drop=T)
  sgRNA_lfc <- do.call(cbind,lapply(groups,function(y)rowMeans(sgRNA_lfc_reps[,y,drop=F])))
  sgRNA_lfc <- cbind(sgRNA_lfc_reps[,which(colnames(sgRNA_lfc_reps)%in%c("sgRNA","Gene","Category")),drop=F],sgRNA_lfc)
  gene_lfc <- do.call(cbind,lapply(groups,function(y)rowMeans(gene_lfc_reps[,y,drop=F])))
  gene_lfc <- cbind(gene_lfc_reps[,which(colnames(gene_lfc_reps)%in%c("Gene","Category")),drop=F],gene_lfc)
  conditions <- unique(annotation$Condition)
  # Two-tail RRA --------------------------------------------------------------------------------------------------------
  # Two-tail p values and FDR of essentiality:
  p <- list()
  for(C in conditions) {
    print(paste0("Computing essentiality P-value on ",C,"..."))
    if(!is.null(shinyF)) shinyF(1.0/length(conditions),message='2-tail RRA',detail=paste0("Computing essentiality P-value on ",C))
    p[[C]] <- merge(gene_lfc[,which(colnames(gene_lfc)%in%c("Gene","Category",C))],getPValues(sgRNA_lfc$Gene,sgRNA_lfc$sgRNA,sgRNA_lfc[,C],0.5))
    p[[C]] <- p[[C]][order(p[[C]]$p,p[[C]]$rho),]
    p[[C]]$FDR <- stats::p.adjust(p[[C]]$p,method="fdr")
    p[[C]]$depleted[p[[C]][,C]<0] <- 1:sum(p[[C]][,C]<0)
    p[[C]]$enriched[p[[C]][,C]>0] <- 1:sum(p[[C]][,C]>0)
    if(!is.null(id)) utils::write.table(p[[C]],paste0(id,"/",C,".RRA.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
  }
  # Control gene and non-differential gene filtering ---------------------------------------------------------------------
  # Find significant genes in each condition:
  # a) Small p value and negative fold change = "low" (essential).
  # b) Large p value = "mid" (nonessential).
  # c) Small p value and positive fold change = "high" (enriched).
  rownames(gene_lfc) <- gene_lfc$Gene
  low <- mid <- high <- gene_lfc[,conditions,drop=FALSE]*0
  for(C in colnames(low)) {
    low[,C] <- as.numeric(rownames(low) %in% p[[C]]$Gene[p[[C]][,C]<0 & p[[C]]$p<pvalue])
    mid[,C] <- as.numeric(rownames(mid) %in% p[[C]]$Gene[p[[C]]$p>pvalue])
    high[,C] <- as.numeric(rownames(high) %in% p[[C]]$Gene[p[[C]][,C]>0 & p[[C]]$p<pvalue])
  }
  # Find genes which are not enriched in any condition:
  high$Any <- apply(high[,,drop=F],1,function(x)sum(x)>0)
  unenriched <- rownames(high)[high$Any==FALSE]
  # Output -------------------------------------------------------------------------------------------------------------------
  if(!is.null(id)) {
    utils::write.table(data.frame(Unenriched=unenriched),paste0(id,"/unenriched.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
  }
  print("Done.")
  return(list(sgRNA_lfc=sgRNA_lfc, gene_lfc_reps=gene_lfc_reps, gene_lfc=gene_lfc, p.value=p, Unenriched=unenriched))
}
