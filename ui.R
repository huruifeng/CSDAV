#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(visNetwork)
library(shinyBS)
library(igraph)
library(shinyWidgets)



# Define UI for application that draws a histogram
ui = fluidPage(
  includeCSS("www/mystyle.css"),
  navbarPage(
    title = "CSDAV",
    id = "CSDAV",
    theme = shinythemes::shinytheme("cerulean"),
    
    # Tab: FASTQ alignment ==========================
    tabPanel("FASTQ alignment",
             style = "padding:0 10px;",
             fluidRow(
               column(width = 4,
                 div(
                   style = "display:inline-block",
                   shinyWidgets::dropdownButton(
                     inputId = "mydropdown1",
                     label = "Objectives",
                     status = "primary",
                     width = "400px",
                     circle = FALSE,
                     tags$ul(
                       tags$li("To map the reads in all FASTQ files contained in a specified directory to a custom sgRNA library."),
                       tags$li("To visualize the sequencing depth and percentage of reads mapped in each file."),
                       tags$li("To annotate the conditions and replicates of the experiment.")
                     )
                   )
                 ),
                 div(
                   style = "display:inline-block",
                   shinyWidgets::dropdownButton(
                     inputId = "mydropdown2",
                     label = "Instructions",
                     status = "primary",
                     width = "400px",
                     circle = FALSE,
                     tags$ol(
                       tags$li("Select a working directory in which all output is to be stored."),
                       tags$li("Select the directory containing the FASTQ files to be mapped."),
                       tags$li("Load an sgRNA library file containing the following two columns: 'sgRNA' and 'Gene'."),
                       "If the library includes gene categories, please add a third column called 'Category'.",
                       tags$li(
                         "Select the starting and ending position of the spacer in the FASTQ file sequences, and specify whether the spacer is to be reversed and/or complemented."
                       ),
                       tags$li("Click on RUN."),
                       tags$li("Once finished, go to the 'File annotation' tab and follow the instructions shown.")
                       # fill out the condition and replicate information in the table provided. Finally, click on the red button above the table to save it.")
                     )
                   )
                 ),
                 wellPanel(
                   style = "margin-top:10px;",
                   tags$style(".popover{max-width: 100%;}"),
                   fluidRow(column(7, uiOutput("dir_out_text1"), textOutput("dir_out_text2")),
                            column(3,
                              shinyFiles::shinyDirButton(
                                "dir_out",
                                label = "Browse...",
                                title = "Please select a working directory:",
                                class = "btn-secondary"
                              )
                            )),
                   HTML("<br>"),
                   fluidRow(column(7, uiOutput("dir_fastq_text1"),textOutput("dir_fastq_text2")),
                            column(3,
                              shinyFiles::shinyDirButton(
                                "dir_fastq",
                                label = "Browse...",
                                title = "Please select the directory containing FASTQ files:",
                                class = "btn-secondary"
                              )
                            )),
                   HTML("<br>"),
                   
                   fileInput(
                     "file.library",
                     label = p(
                       "Please load an sgRNA library file:",
                       downloadLink("q2", label = "Format references", class = "btn btn-primary btn-xs")
                     )
                   ),
                  
                   selectInput("program_option",
                                 p("Select the program for alignment:"), 
                                 choices = list("MAGeCK count" = "MAGeCK", "MoPAC" = "MoPAC"), selected = "MAGeCK"
                    ),
                 
                   conditionalPanel(
                     condition = "input.program_option == 'MAGeCK'",
                     sliderInput( "gRNA.length","sgRNA length range:",
                                  min = 10,
                                  max = 30,
                                  value = c(19, 21),
                                  step = 1
                     )
                   ),
              
                   conditionalPanel(
                     condition = "input.program_option == 'MoPAC'",
                     sliderInput( "spacer.range","Spacer location:",
                                  min = 1,
                                  max = 50,
                                  value = c(1, 19),
                                  step = 1
                     ),
                     fluidRow(
                       column(6, checkboxInput("spacer.rev", label = "Reverse.", value = T)),
                       column( 6, checkboxInput("spacer.comp", label = "Complement.", value = T))
                     ),
                   ),
                   
                  
                   
                  
                   fluidRow(
                     align = "center",
                     actionButton("button_run", label = "RUN", class = "btn-primary")
                   )
                 )
                 
               ), ## END Col4
               column(8,tabsetPanel(
                 tabPanel("Charts",
                          plotOutput("mapped2"),
                          plotOutput("mapped1")
                 ),
                 tabPanel("File annotation",
                          HTML("<br>"),uiOutput("text_anno1"),HTML("<br>"),
                          fluidRow(align="center",actionButton("save_anno","Click here when finished.",class="btn-danger")),
                          br(),
                          rhandsontable::rHandsontableOutput("anno_table")),
                 tabPanel("Table",
                          wellPanel(fluidRow(
                            column(6,radioButtons("fastq.level",label="Level:",choices=list("FASTQ files","Conditions"),selected="FASTQ files",inline=T)),
                            DT::dataTableOutput(outputId="data_input"))
                          ))
               )),  ## END Col8
               tags$head(tags$style(".shiny-output-error{color: blue;}")),
             ), ## fluidRow
             
             shinyBS::bsPopover("q1",content=HTML(paste("sgRNA &nbsp;&nbsp;&nbsp;&nbsp; Gene &nbsp;&nbsp;&nbsp;&nbsp; Category &nbsp;&nbsp;&nbsp;&nbsp; sample1 &nbsp;&nbsp;&nbsp;&nbsp; sample2 &nbsp;&nbsp;&nbsp;&nbsp;",
                                                        "sgrna1 &nbsp;&nbsp;&nbsp;&nbsp; gene1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; category1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; .............. &nbsp;&nbsp;&nbsp;&nbsp; ..............",
                                                        "sgrna2 &nbsp;&nbsp;&nbsp;&nbsp; gene2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; category2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; .............. &nbsp;&nbsp;&nbsp;&nbsp; ..............",
                                                        "sgrna3 &nbsp;&nbsp;&nbsp;&nbsp; gene3 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; category3 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; .............. &nbsp;&nbsp;&nbsp;&nbsp; ..............",
                                                        sep="<br/>")),title="Format of counts table.",placement="right",options=list(container="body")),
             shinyBS::bsPopover("q2",content=HTML(paste("sgRNA  &nbsp;&nbsp;&nbsp;      Gene &nbsp;&nbsp;&nbsp;&nbsp;      Category",
                                                        "sgrna1 &nbsp;&nbsp;&nbsp;      gene1 &nbsp;&nbsp;&nbsp;      category1",
                                                        "sgrna2 &nbsp;&nbsp;&nbsp;      gene2 &nbsp;&nbsp;&nbsp;      category2",
                                                        "sgrna3 &nbsp;&nbsp;&nbsp;      gene3 &nbsp;&nbsp;&nbsp;      category3",
                                                        sep="<br/>")),title="sgRNA Library Format",placement="right",options=list(container="body")),
             shinyBS::bsPopover("q3",content=HTML(paste("sgRNA  &nbsp;&nbsp;&nbsp;      Gene &nbsp;&nbsp;&nbsp;&nbsp;      Category",
                                                        "sgrna1 &nbsp;&nbsp;&nbsp;      gene1 &nbsp;&nbsp;&nbsp;      category1",
                                                        "sgrna2 &nbsp;&nbsp;&nbsp;      gene2 &nbsp;&nbsp;&nbsp;      category2",
                                                        "sgrna3 &nbsp;&nbsp;&nbsp;      gene3 &nbsp;&nbsp;&nbsp;      category3",
                                                        sep="<br/>")),title="sgRNA Library Format",placement="right",options=list(container="body")),
             ),
    ## END tab: FASTQ alignment
    
    # Tab: Quality control =======================
    tabPanel("Check and Annotate",
             style = "padding:0 10px;",
             fluidRow(
                column(4,
                    div(style="display:inline-block",
                        shinyWidgets::dropdownButton(
                          inputId = "mydropdown3",
                          label = "Objectives",
                          status = "primary",
                          width="400px",
                          circle = FALSE,
                          tags$ul(
                            tags$li("To visualize the quality of the experiment."),
                            tags$li("To annotate the control samples for each condition."),
                          )
                        )),
                    div(style="display:inline-block",
                        shinyWidgets::dropdownButton(
                          inputId = "mydropdown4",
                          label = "Instructions",
                          status = "primary",
                          width="400px",
                          circle = FALSE,
                          tags$ol(
                            tags$li("If you skipped the FASTQ Alignment module:"),
                            tags$ul(
                              tags$li("Select a working directory in which all output is to be stored."),
                              tags$li("Upload a table of read counts, which must include a column named 'sgRNA'."),
                              tags$li("If the table of read counts is missing the 'Gene' column, please also upload an sgRNA library file containing the following two columns: 'sgRNA' and 'Gene'."),
                              tags$li("If the library includes gene categories, please add a column called 'Category' to either the table of read counts or the library file."),
                              tags$li("Click on Check"),
                              tags$li("Once finished, go to the 'Condition annotation' tab and follow the instructions shown.")# fill out the condition and replicate information in the table provided. Finally, click on the red button above the table to save it.")
                            ),
                            tags$li("Otherwise, proceed directly to the 'Condition annotation' tab and follow the instructions shown.")# specify which is the control sample (i.e. 'Day 0') for each condition in the table provided. Finally, click on the red button above the table to save it.")
                          )
                        )),
                    wellPanel(
                      style = "margin-top:10px;",
                      tags$style(".popover{max-width: 100%;}"),
                      fluidRow(column(7,uiOutput("dir_qc_text1"), textOutput("dir_qc_text2")),
                               column(3,shinyFiles::shinyDirButton("dir_qc",label="Browse...",title="Select working directory",class="btn-secondary"))),
                      HTML("<br>"),
                      conditionalPanel(
                        condition="output.fastq == 0 && output.counts_path==0",
                        fileInput("file.counts",label=p("Read counts file:",downloadLink("q1",label = "Format references",class="btn btn-primary btn-xs"))),
                        fileInput("file.library2",label=p("sgRNA library file:",downloadLink("q3",label = "Format references" ,class="btn btn-primary btn-xs")))
                      ),
                      
                      conditionalPanel(
                        condition="output.counts_path!=0",
                        uiOutput("counts_path_txt"),HTML("<br>"),
                        uiOutput("library_path_txt"),HTML("<br>")
                      ),
                      
                      conditionalPanel(
                        condition="output.fastq != 0",
                        uiOutput("counts_qc"),HTML("<br>"),
                        uiOutput("library_qc"),HTML("<br>")
                      ),
                      conditionalPanel(
                        condition="output.fastq == 0",
                        fluidRow(align="center",actionButton("button_check",label="Check",class="btn-primary"))
                      )
                    )),
             column(8,tabsetPanel(
               tabPanel("Charts",
                        plotOutput("classification"),
                        plotOutput("sizes"),
                        plotOutput("depths"),
                        plotOutput("zeros"),
                        plotOutput("gini")
               ),
               tabPanel("Condition annotation",
                        HTML("<br>"),uiOutput("text_anno2"),HTML("<br>"),
                        rhandsontable::rHandsontableOutput("anno1_table"),
                        br(),
                        fluidRow(align="center",actionButton("save_anno1","Save Annotation",class="btn-danger")),
                        br()
              )
             ))
          )## fluidRow
      ),
    ## END Tab QC
    
    # Tab: Analysis: single condition =======================
    tabPanel("Analysis: single condition",
             column(4,
                    div(style="display:inline-block",
                        shinyWidgets::dropdownButton(
                          inputId = "mydropdown5",
                          label = "Objectives",
                          status = "primary",
                          width="400px",
                          circle = FALSE,
                          tags$ul(
                            tags$li("To compute the sgRNA and gene log-fold-changes."),
                            tags$li("To perform a 2-tail a-RRA analysis of gene essentiality."),
                            tags$li("To compute a condition-specific measure of gene essentiality score.")
                          )
                        )),
                    div(style="display:inline-block",
                        shinyWidgets::dropdownButton(
                          inputId = "mydropdown6",
                          label = "Instructions",
                          status = "primary",
                          width="400px",
                          circle = FALSE,
                          tags$ol(
                            tags$li("If a subset of the samples are to be used, uncheck the corresponding box."),
                            tags$li("If you wish to modify the RRA threshold to identify enriched genes, uncheck the corresponding box. The normalization is first carried out on the unenriched genes and then interpolated to all genes."),
                            tags$li("Click on RUN.")
                          )
                        )),
                    wellPanel(
                      style = "margin-top:10px;",
                      tags$style(".popover{max-width: 100%;}"),
                
                      uiOutput("dir_ea_text1"),
                      textOutput("dir_ea_text2"),HTML("<br>"),
                      
                      selectInput("program_option_ea",
                                  p("1. Select the program for analysis:"), 
                                  choices = list("MAGeCK RAA" = "MAGeCKRRA","MAGeCK MLE" = "MAGeCKMLE", "MoPAC LFC&RRA" = "MoPACRRA"), selected = "MAGeCKRRA"
                      ),
                      
                      conditionalPanel(condition="input.program_option_ea=='MoPACRRA'",
                                wellPanel( "MoPAC Parameters:",
                                       checkboxInput("pseudo_default",label="Use default pseudo-count = 0.1.",value=TRUE),
                                       conditionalPanel(condition="!input.pseudo_default",
                                                        sliderInput("pseudo.count",label="Choose pseudo-count:",value=0.1,min=0.01,max=0.3,step=0.01)),
                                       
                                       checkboxInput("rra_default",label="Use default RRA P-value threshold = 0.05",value=TRUE),
                                       conditionalPanel(condition="!input.rra_default",
                                                        # checkboxInput("read.filtered",label="Filter control genes with 2-tail RRA.",value=F),
                                                        sliderInput("pvalue",label="P value threshold:",value=0.05,min=0.0,max=0.1,step=0.01)),
                                       HTML("<br>"),
                                       checkboxInput("qc.report",label="Generate PDF report (requires pandoc installation).",value=FALSE)
                                ),
                      ),
                     conditionalPanel(condition="input.program_option_ea=='MAGeCKRRA'",
                              wellPanel("MAGeCK RRA Parameters:",
                                      selectInput("MAGeCKRRA.norm",
                                                  "Normalization method:", 
                                                  choices = list("None" = "None","Median"= "Median","Log" = "Log"), selected = "Median"
                                      ),
                              ),
                     ),
                     conditionalPanel(condition="input.program_option_ea=='MAGeCKMLE'",
                              wellPanel("MAGeCK MLE:",
                                      selectInput("MAGeCKMLE.param",
                                                  "Params:", 
                                                  choices = list("X" = "X","Y"= "Y","Z" = "Z"), selected = "X"
                                      ),
                              ),
                     ),
                     
                     selectInput("samples.case","2. Select case group:",choices = NULL, multiple=TRUE, selected = 1),
                     checkboxInput("use_reps",label="Use all case replicates.",value=T),
                     uiOutput("use_reps_case"),
                     
                     
                     HTML("<br>"),
                     selectInput("samples.ctrl","3. Select control group:",choices = NULL, multiple=TRUE, selected = 1),
                     checkboxInput("use_reps_ctrl",label="Use all control replicates.",value=T),
                     uiOutput("use_reps_ctrl"),
                     
                     HTML("<br>"),
                     fluidRow(align="center",actionButton("button_ea",label="Run Analysis",class="btn-primary")),
                     HTML("<br>"),
                    ),
                    uiOutput("inforra")
             ),
             column(8,tabsetPanel(
               tabPanel("Volcano",
                        wellPanel(
                          selectInput("rra.condition",choices="",label="Visualize essentiality p values of condition:"),
                          HTML("<hr>"),
                          selectizeInput("rragene",label="Highlight genes:",choices=NULL,multiple=TRUE,selected=1,width="100%"),
                          checkboxInput("rracolor",label="Highlight categories.",value=F),
                          HTML("<br>")
                        ),
                        fluidRow(
                          column(6,div(style="position:relative",plotOutput("rra",width="100%",height="auto",hover=hoverOpts("rra_hover",delay=10,delayType="debounce"),click="clickrra",brush=brushOpts(id="plotrra_brush",resetOnNew=FALSE)),uiOutput("plotrrainfo"))),
                          column(6,div(style="position:relative",plotOutput("rrazoom",width="100%",height="auto",hover=hoverOpts("rrazoom_hover",delay=10,delayType="debounce"),click="clickrrazoom"),uiOutput("plotrrazoominfo")))
                        )),
               tabPanel("Distributions",
                        wellPanel(
                          selectInput("lfc.condition",choices="",label="Condition:"),
                          fluidRow(
                            column(4,radioButtons("lfc.style",label="Visualization:",choices=c("Box plot","Density plot"),selected="Box plot",inline=T)),
                            column(4,radioButtons("folds.level",label="Level:",choices=c("sgRNA","Gene"),selected="Gene",inline=T)),
                            column(4,radioButtons("log.lfc",label="Measurement:",choices=c("Log-read-count","Log-fold-change"),selected="Log-fold-change",inline=T))
                          )),
                        plotOutput("folds",width="100%",height="auto")#,
                        # plotOutput("ssmd")
               ),
               tabPanel("Correlation/PCA",
                        wellPanel(fluidRow(
                          column(4,radioButtons("corr_pca",label="Analysis:",choices=list("Correlation","Principal component analysis"),selected="Correlation",inline=F)),
                          column(4,radioButtons("corr_pca_method",label="Method:",choices=c("pearson","spearman"),selected="pearson",inline=F)),
                          column(4,radioButtons("corr.level",label="Level:",choices=c("sgRNA","Gene"),selected="Gene",inline=F))
                        )),
                        conditionalPanel(condition="input.corr_pca=='Correlation'",
                                         plotOutput("correlation",width="100%",height="auto")),
                        conditionalPanel(condition="input.corr_pca=='Principal component analysis'",
                                         plotOutput("pca1",width="100%",height="auto"),
                                         plotOutput("pve"))
               ),
               tabPanel("Reproducibility",
                        wellPanel(selectInput("qc.condition",choices="",label="Condition:"),
                                  fluidRow(
                                    column(6,radioButtons("reproducibility.level",label="Level:",choices=c("sgRNA","Gene"),selected="Gene",inline=T)),
                                    column(6,radioButtons("reproducibility.stage",label="Measurement:",choices=c("Log-read-count","Log-fold-change"),selected="Log-fold-change",inline=T))
                                  )),
                        plotOutput("reproducibility",width="100%",height="auto")
               ),
               tabPanel("LFC Table",
                        wellPanel(fluidRow(
                          column(4,radioButtons("qc.level",label="Level:",choices=list("sgRNA","Gene"),selected="sgRNA",inline=T)),
                          column(4,radioButtons("qc.stage",label="Measurement:",choices=c("Log-read-count","Log-fold-change"),selected="Log-fold-change",inline=T)),
                          column(4,checkboxInput("qc.averaged",label="Replicate-averaged.",value=T))
                        )),
                        DT::dataTableOutput(outputId="data_qc")
               ),
               tabPanel("Table",
                        wellPanel(fluidRow(
                          selectInput("ea.condition",choices="",label="Condition:")
                          # column(6,checkboxInput("ea.averaged",label="Replicate-averaged.",value=T))
                        )),
                        DT::dataTableOutput(outputId="data_ea")
               ))
             )
    ),
    ## END Tab Analysis: single condition
    
    # Tab: Analysis: differential =======================
    tabPanel("Analysis: differential",
             column(4,
                    div(style="display:inline-block",
                        shinyWidgets::dropdownButton(
                          inputId = "mydropdown7",
                          label = "Objectives",
                          status = "primary",
                          width="400px",
                          circle = FALSE,
                          tags$ul(
                            tags$li("To normalize the distribution of one pair of conditions."),
                            tags$li("To generate differential gene essentiality scores for the pair of conditions selected.")
                          )
                        )),
                    div(style="display:inline-block",
                        shinyWidgets::dropdownButton(
                          inputId = "mydropdown8",
                          label = "Instructions",
                          status = "primary",
                          width="400px",
                          circle = FALSE,
                          tags$ol(
                            tags$li("Select the two conditions to analyze."),
                            tags$li("Choose whether you wish to generate an html report and have 'plotly' installed."),
                            tags$li("Click on RUN.")
                          )
                        )),
                    wellPanel(
                      style = "margin-top:10px;",
                      tags$style(".popover{max-width: 100%;}"),
                    
                      fluidRow(column(7, uiOutput("dir_da_text1"),textOutput("dir_da_text2")),
                               column(3,shinyFiles::shinyDirButton("dir_da",label="Browse...",title="Select working directory",class="btn-secondary"))),
                      uiOutput("dir_da_err"),
                       HTML("<br>"),
                      
                      selectInput("condition1",choices="",label="Essential analysis 1:"),
                      selectInput("condition2",choices="",label="Essential analysis 2:"),
                      br(),
                      checkboxInput("da.report",label="Generate HTML report (requires pandoc and plotly installation).",value=FALSE),
                      fluidRow(align="center",actionButton("button_da",label="RUN",class="btn-primary")),
                      hr(),
                      
                      
                      selectizeInput("dagene",label="Highlight genes:",choices=NULL,multiple=TRUE,selected=1,width="100%"),
                      checkboxInput("diffcolor",label="Highlight categories.",value=F),
                      # radioButtons("normalization",label="Normalization based on:",choices=list("Nondifferentials","Controls"),selected="Nondifferentials",inline=T),
                     
                    ),
                    uiOutput("infoscore")
             ),
             column(8,tabsetPanel(
               tabPanel("Volcano",
                        wellPanel(
                          radioButtons("p.volcano",label="P-value source:",choices=c("sgRNA reproducibility","Biological reproducibility"),selected="sgRNA reproducibility",inline=T)
                          ),
                        fluidRow(
                          column(6,div(style="position:relative",plotOutput("plot3",width="100%",height="auto",hover=hoverOpts("plot3_hover",delay=10,delayType="debounce"),click="click3",brush=brushOpts(id="plot3_brush",resetOnNew=FALSE)),uiOutput("plot3info"))),
                          column(6,div(style="position:relative",plotOutput("plot4",width="100%",height="auto",hover=hoverOpts("plot4_hover",delay=10,delayType="debounce"),click="click4"),uiOutput("plot4info")))
                        )),
               tabPanel("Gene scatter",
                        wellPanel(checkboxInput("gscatter",label="Use normalized values.",value=T)),
                        fluidRow(
                          column(6,div(style="position:relative",plotOutput("plot1",width="100%",height="auto",hover=hoverOpts("plot1_hover",delay=10,delayType="debounce"),click="click1",brush=brushOpts(id="plot1_brush",resetOnNew=FALSE)),uiOutput("plot1info"))),
                          column(6,div(style="position:relative",plotOutput("plot2",width="100%",height="auto",hover=hoverOpts("plot2_hover",delay=10,delayType="debounce"),click="click2"),uiOutput("plot2info")))),
                        plotOutput("weights"),
                        plotOutput("scaling"),
                        plotOutput("shifting")
               ),
               tabPanel("sgRNA scatter",
                        wellPanel(checkboxInput("sscatter",label="Use normalized values.",value=T)),
                        fluidRow(
                          column(6,div(style="position:relative",plotOutput("plot5",width="100%",height="auto",hover=hoverOpts("plot5_hover",delay=10,delayType="debounce"),click="click5",brush=brushOpts(id="plot5_brush",resetOnNew=FALSE)),uiOutput("plot5info"))),
                          column(6,div(style="position:relative",plotOutput("plot6",width="100%",height="auto",hover=hoverOpts("plot6_hover",delay=10,delayType="debounce"),click="click6"),uiOutput("plot6info")))
                        )),
               tabPanel("Table",
                        wellPanel(radioButtons("data_da",label="Level:",choices=c("sgRNA","Gene"),selected="Gene",inline=T)),
                        DT::dataTableOutput(outputId="table_da")
               )
             ))
    ),
    
    
    
    
    
  ) ## END navbarPage
)


