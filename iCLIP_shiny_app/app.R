#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(RColorBrewer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(htmltools)

theme_set(theme_classic()) 
theme_update(plot.background = element_blank(), 
             panel.background = element_blank(), 
             legend.background = element_blank(), 
             legend.key = element_blank(),
             axis.line = element_blank())

load("CLIP_tool_data.RData")

js <- '
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
    Shiny.onInputChange("keyPressed", Math.random());
  }
});
'

# Define UI for application that draws a histogram
ui <- fluidPage(
    
  
    tags$script(js), #records enter being pressed
    # Application title
    titlePanel("CLIP data visualisation"),
    
    br(),
    br(),

    wellPanel(fluidRow(
      column(width=4,
        textInput("gene",label = "Please provide mouse gene name or Ensembl ID:"),
        actionButton("update_gene", "Submit"),
        uiOutput("transcript")
        ),
      column(width=4,
        radioButtons("CLIP_type", "Select which CLIP datasets to display:", choiceNames = unique(CLIP_datasets$Type), choiceValues = unique(CLIP_datasets$Type), selected = unique(CLIP_datasets$Type)[1]),
        uiOutput("CLIP_datasets"),
    checkboxInput("CLIPscale", "Display all CLIP datasets on the same scale", value = T),
    checkboxInput("mergeReplicates", "Merge replicate datasets", value = T),
    checkboxInput("clusters", "Display clusters below each track", value = F)
      ),
      column(width=4,
        checkboxInput("FilterFDR", "Apply FDR/score thresholds", value = F),
    uiOutput("FDRthresh"),
    uiOutput("Scorethresh"),
    uiOutput("nReplicates"),
    uiOutput("ShowAll"),
    splitLayout(cellWidths = c("70%", "30%"),
                p(strong("Enter motif(s) to display below plot")),
            actionButton("titleBtId", "", icon=icon('question-circle'),class = "btn-xs", title = "Info")
            ),
        textAreaInput("motif", label = NULL, rows = 3),
        actionButton("update_motif", "Update motif(s)"),
        checkboxInput("rev_scale", "Reverse orientation of plots", value = F)
      ),
    )),
    

    
    h2("Complete gene"),
    div(  htmlDependency("font-awesome", 
                         "5.13.0", "www/shared/fontawesome", package = "shiny", 
                         stylesheet = c("css/all.min.css", "css/v4-shims.min.css")),
          fluidRow(
      column(width=8,
        textOutput("selected_gene"),
    plotOutput("genePlot", height="auto")
      ),
    column(width=2,
           wellPanel(
             radioButtons("CLIPcolour", "Colour for CLIP crosslinks", choiceNames =list(HTML('<div style="display:flex"><div class="fa fa-plus-square"
                                         style="color:#cd2626;margin-top:3px;"></div><div class="fa fa-minus-square"
                                         style="color:#3a5fcd;margin-top:3px;"></div></div>'),
                                          HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#cd2626;margin-top:3px;"></div></div>'),
                                          HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#228b22;margin-top:3px;"></div></div>'),
                                          HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#3a5fcd;margin-top:3px;"></div></div>'),
                                          HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#000080;margin-top:3px;"></div></div>'),
                                          HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:black;margin-top:3px;"></div></div>'), "Other"), 
                                choiceValues = c("default", "firebrick3", "forestgreen", "royalblue3", "navy", "black", "other")),
             uiOutput("CLIPcol_text"),
             radioButtons("GENEcolour", "Colour for transcripts", choiceNames =list(HTML('<div style="display:flex"><div class="fa fa-plus-square"
                                         style="color:#cd2626;margin-top:3px;"></div><div class="fa fa-minus-square"
                                         style="color:#3a5fcd;margin-top:3px;"></div></div>'),
                                                                                           HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#cd2626;margin-top:3px;"></div></div>'),
                                                                                           HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#228b22;margin-top:3px;"></div></div>'),
                                                                                           HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#3a5fcd;margin-top:3px;"></div></div>'),
                                                                                           HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#000080;margin-top:3px;"></div></div>'),
                                                                                           HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:black;margin-top:3px;"></div></div>'), "Other"), 
                          choiceValues = c("default", "firebrick3", "forestgreen", "royalblue3", "navy", "black", "other")),
             uiOutput("GENEcol_text")
           )
           ),
      column(width=2,
      wellPanel(
        radioButtons("CLUSTERcolour", "Colour for clusters", choiceNames =list(HTML('<div style="display:flex"><div class="fa fa-plus-square"
                                         style="color:#cd2626;margin-top:3px;"></div><div class="fa fa-minus-square"
                                         style="color:#3a5fcd;margin-top:3px;"></div></div>'),
                                                                                      HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#cd2626;margin-top:3px;"></div></div>'),
                                                                                      HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#228b22;margin-top:3px;"></div></div>'),
                                                                                      HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#3a5fcd;margin-top:3px;"></div></div>'),
                                                                                      HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:#000080;margin-top:3px;"></div></div>'),
                                                                                      HTML('<div style="display:flex"><div class="fa fa-square"
                                         style="color:black;margin-top:3px;"></div></div>'), "Other"), 
                     choiceValues = c("default", "firebrick3", "forestgreen", "royalblue3", "navy", "black", "other")),
        uiOutput("CLUSTERcol_text"),
        sliderInput("pWidth", "Specify plot width (mm)", min =100, max = 500, value = 230),
          uiOutput("height_slider"),
          downloadButton("PDFcomplete", "Export PDF"),
          p(),
          downloadButton("PNGcomplete", "Export PNG")
        )
      )
    )),
    
    
    br(),
    br(),
    
    wellPanel(fluidRow(
      column(width=6,
        radioButtons("feature_zoom", "Select feature for zoomed view (new selection will override coordinates). Note that feature view will comprise the selected feature of all selected isoforms.", choiceNames = c("None", "5'UTR", "3'UTR", "CDS"), choiceValues = c("None", "five_prime_utr", "three_prime_utr", "CDS"), selected = "None"),
        br(),
        p("Note that if clusters are displayed, only those entirely contained within the zoomed in view will be shown")
      ),
      column(width=6,
        p(strong("Specify coordinates for zoomed view (update will override feature selection)")),
        numericInput("coordMin_zoom", "Minimum:", value = 0),
        numericInput("coordMax_zoom", "Maximum:", value = 0),
        actionButton("update_coord", "Update coordinates")    
      ),
    )),
    
    
    #actionButton("browser", "browser"),
    br(),
    h2("Zoomed in view"),
    fluidRow(
      column(width=10,
        plotOutput("genePlotZoom", height="auto")
      ),
      column(width=2,
          wellPanel(
            downloadButton("PDFzoom", "Export PDF"),
            p(),
            downloadButton("PNGzoom", "Export PNG")
            )
      ),
    ),
    br(),
    #tableOutput("tableTest")
    
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  rv <- reactiveValues(
    zoom_max = 0,
    zoom_min = 0,
    motif_lengths = numeric(),
    motifs = character(),
    motif_cols = brewer.pal(8, "Set2")[c(2,3,1,4:8)],
    chromosome = "",
    gene = "",
    strand = "+",
    use_coord = F, 
    input_gene = "",
    gene_transcripts = slice(gene_names_ids,0),
    gene_summary_subset = slice(gene_summary, 0),
    invalid_gene = F,
    plotting_colour = "firebrick3",
    gene_colour = "firebrick3",
    cluster_colour = "firebrick3",
    FDRthreshold = 1,
    ScoreThreshold = 0,
    display_all = F, 
    replicate_threshold = 2,
    CLIPpath = "CLIP_data/"
  )
  
  observeEvent(input$CLIP_type, {
    if(input$CLIP_type == "ZFP36-family CLIP") {
      rv$CLIPpath <- "CLIP_data/"
    } else {
      rv$CLIPpath <- "m6A_data/"
    }
  })
  
  observeEvent(input$titleBtId, {
    showModal(modalDialog(
      title = "Note",
      "Standard degenerate base codes can be used. Flexible numbers of bases can also be included within curly brackets, eg T{1-2} will allow 1 or 2 Ts at that position; N{0-3} will allow an insert of up to 3 bases of any type. By default the longest possible match will be returned - to return the shortest add ? eg N{0-3}?. Overlapping matches will be found.",
      easyClose = TRUE
    ))
  })
  observeEvent(input$update_gene | input$keyPressed, {
    if(input$gene != rv$input_gene) {
      if(grepl("^ENSMUSG", input$gene, ignore.case = T)) {
        gene_names_ids %>%
          filter(toupper(gene_id) == toupper(input$gene)) %>%
          unite(transcript_id_name, transcript_id, transcript_name, transcript_biotype, sep = "; ", remove = F) -> rv$gene_transcripts
      } else {
        gene_names_ids %>%
          filter(toupper(gene_name) == toupper(input$gene)) %>%
          unite(transcript_id_name, transcript_id, transcript_name, transcript_biotype, sep = "; ", remove = F) -> rv$gene_transcripts
      }
      gene_summary %>%
        filter(gene_name == rv$gene_transcripts$gene_name[1]) -> rv$gene_summary_subset
      if(nrow(rv$gene_summary_subset) > 0) {
        rv$invalid_gene <- F
        if(rv$gene_summary_subset$gene_name != rv$gene) { 
          rv$use_coord <- F
          rv$gene <- rv$gene_summary_subset$gene_name 
        }
        if(rv$gene_summary_subset$chromosome != rv$chromosome) { rv$chromosome <- rv$gene_summary_subset$chromosome }
        if(rv$gene_summary_subset$strand != rv$strand) { rv$strand <- rv$gene_summary_subset$strand }
        if(input$CLIPcolour == "default") {
          rv$plotting_colour <- if_else(rv$strand == "+", "firebrick3", "royalblue3")
        }
        if(input$GENEcolour == "default") {
          rv$gene_colour <- if_else(rv$strand == "+", "firebrick3", "royalblue3")
        }
        if(input$CLUSTERcolour == "default") {
          rv$cluster_colour <- if_else(rv$strand == "+", "firebrick3", "royalblue3")
        }
      } else { rv$invalid_gene <- T }
      rv$input_gene <- input$gene
    }
    
  })
  
    output$selected_gene <- renderText({
      if(rv$invalid_gene == T) {
        if(nchar(isolate(input$gene)) > 0) {
                paste(isolate(input$gene), "not found. Please select a valid Ensembl gene name or ID (GRCm38v101). Only genes with at least one isoform with CCDS and/or transcript support level 1 will be recognised")
            } 
        }
        
    })
    
 
    output$transcript <- renderUI({
        rv$gene
        if(nchar(isolate(input$gene)) == 0) { } else {
            checkboxGroupInput("transcript_selection", "Select transcript isoform(s) to plot:", choiceNames = rv$gene_transcripts$transcript_id_name, choiceValues = rv$gene_transcripts$transcript_id, selected = rv$gene_transcripts$transcript_id)
        }
    })
    
    output$CLIP_datasets <- renderUI({
      checkboxGroupInput("CLIPdatasets", NULL, choiceNames = unique(CLIP_datasets$dataset[CLIP_datasets$Type == input$CLIP_type]), choiceValues = unique(CLIP_datasets$dataset[CLIP_datasets$Type == input$CLIP_type]), selected = if(input$CLIP_type == "ZFP36-family CLIP") { unique(CLIP_datasets$dataset[CLIP_datasets$Type == input$CLIP_type])[1:3] } else { unique(CLIP_datasets$dataset[CLIP_datasets$Type == input$CLIP_type])[1]})
    })
    
    output$height_slider <- renderUI({
      sliderInput("pHeight", "Specify plot height (mm)", min =max(50, round((plot_height()/2.83)-200)), max = round((plot_height()/2.83)+200), value = round((plot_height()/2.83)))
    })
    
    output$CLIPcol_text <- renderUI({
      if(input$CLIPcolour == "other") {
        textInput("CLIPother_col", "Enter in hexadecimal format:")
      }
    })
    output$GENEcol_text <- renderUI({
      if(input$GENEcolour == "other") {
        textInput("GENEother_col", "Enter in hexadecimal format:")
      }
    })
    output$CLUSTERcol_text <- renderUI({
      if(input$GENEcolour == "other") {
        textInput("CLUSTERother_col", "Enter in hexadecimal format:")
      }
    })
    
    transcript_coord_gene <- reactive({
        gene_coord %>%
            filter(transcript_id %in% rv$gene_transcripts$transcript_id)
    })
    
    transcript_coord <- reactive({
        transcript_coord_gene() %>%
            filter(transcript_id %in% input$transcript_selection) %>%
            mutate(transcript_offset = as.numeric(as.factor(transcript_name))+0.3)
        
        
    })
    
    intron_coord_gene <- reactive({
        lapply(rv$gene_transcripts$transcript_id, function(x) {
            transcript_coord_gene() %>%
                filter(transcript_id == x) %>%
                dplyr::rename(intron_start = end) -> intron_temp
          if(nrow(intron_temp) > 1) {
            intron_temp %>% 
                mutate(intron_end = c(start[2:length(start)],NA)) %>%
                filter(!is.na(intron_end)) %>%
                filter(intron_start != intron_end)
          } else {
            intron_temp %>% 
              add_column(intron_end = NA_real_) %>%
              filter(!is.na(intron_end)) %>%
              filter(intron_start != intron_end)
          }
        }) %>%
            bind_rows()
    })
    
    intron_coord <- reactive({
        intron_coord_gene() %>%
            filter(transcript_id %in% input$transcript_selection) %>%
            mutate(transcript_offset = as.numeric(as.factor(transcript_name))+0.3)
    })
    
    observeEvent(input$CLIPcolour, {
      if(input$CLIPcolour == "default") {
        rv$plotting_colour <- if_else(rv$strand == "+", "firebrick3", "royalblue3")
      } else if(input$CLIPcolour != "other") { rv$plotting_colour <- input$CLIPcolour }
    })
    
    observeEvent(input$CLIPother_col, {
      if(grepl("^#[A-Fa-f0-9]{6}$", input$CLIPother_col)) {
        rv$plotting_colour <- input$CLIPother_col
      }
    })
    
    observeEvent(input$GENEcolour, {
      if(input$GENEcolour == "default") {
        rv$gene_colour <- if_else(rv$strand == "+", "firebrick3", "royalblue3")
      } else if(input$GENEcolour != "other") { rv$gene_colour <- input$GENEcolour }
    })
    
    observeEvent(input$GENEother_col, {
      if(grepl("^#[A-Fa-f0-9]{6}$", input$GENEother_col)) {
        rv$gene_colour <- input$GENEother_col
      }
    })
    
    observeEvent(input$CLUSTERcolour, {
      if(input$CLUSTERcolour == "default") {
        rv$cluster_colour <- if_else(rv$strand == "+", "firebrick3", "royalblue3")
      } else if(input$CLUSTERcolour != "other") { rv$cluster_colour <- input$CLUSTERcolour }
    })
    
    observeEvent(input$CLUSTERother_col, {
      if(grepl("^#[A-Fa-f0-9]{6}$", input$CLUSTERother_col)) {
        rv$cluster_colour <- input$CLUSTERother_col
      }
    })
    
    
    gene_coord_min <- reactive({
        min(rv$gene_summary_subset$start, rv$gene_summary_subset$end)
    })
    
    gene_coord_max <- reactive({
        max(rv$gene_summary_subset$start, rv$gene_summary_subset$end)
    })
    
    
    observeEvent(input$feature_zoom, {
        rv$use_coord <- F
        if(input$feature_zoom %in% transcript_coord()$feature) {
            transcript_coord() %>%
                filter(feature == input$feature_zoom) %>%
                select(start, end) %>% 
                min()+0.5 -> rv$zoom_min
            transcript_coord() %>%
                filter(feature == input$feature_zoom) %>%
                select(start, end) %>% 
                max()-0.5 -> rv$zoom_max
        } else {
            rv$zoom_max <- 0
            rv$zoom_min <- 0
        }
    })
    
    
    
    observeEvent(input$transcript_selection, {
        if(rv$use_coord == F) {
            if(input$feature_zoom %in% transcript_coord()$feature) {
                transcript_coord() %>%
                filter(feature == input$feature_zoom) %>%
                select(start, end) %>% 
                min()+0.5 -> rv$zoom_min
            transcript_coord() %>%
                filter(feature == input$feature_zoom) %>%
                select(start, end) %>% 
                max()-0.5 -> rv$zoom_max
            } else {
            rv$zoom_max <- 0
            rv$zoom_min <- 0
            }
        }
    })
    
    
    observeEvent(input$update_coord, {
        rv$use_coord <- T
        rv$zoom_min <- min(input$coordMin_zoom, input$coordMax_zoom)
        rv$zoom_max <- max(input$coordMin_zoom, input$coordMax_zoom)
    })
    
    observeEvent(input$update_motif, {
        rv$motifs <- unique(toupper(strsplit(input$motif, "[\n ,;]+")[[1]]))
        rv$motifs <- rv$motifs[rv$motifs != ""] # probably no longer necessary now strsplit will treat multiple consecutive separating characters as single sep
        rv$motif_lengths <- sapply(rv$motifs, nchar)
        if(length(rv$motifs) > length(rv$motif_cols)) {
            rv$motif_cols <- rep(rv$motif_cols, ceiling(length(rv$motifs)/length(rv$motif_cols)))
        }
        })
    
    zoom_length <- reactive(rv$zoom_max - rv$zoom_min)
    coord_min <- reactive(min(select(transcript_coord(), start, end))+0.5) # previously adjusted for plotting
    coord_max <- reactive(max(select(transcript_coord(), start, end))-0.5)
    coord_length <- reactive(coord_max() - coord_min())
    
    complete_seq <- reactive({ 
        Seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, names = paste0("chr", rv$chromosome), start = gene_coord_min(), end = gene_coord_max(), strand = rv$strand, as.character = T)
        
        if(rv$strand == "+") {
            tibble(Base = strsplit(Seq, split= "")[[1]]) %>%
                mutate(Base = factor(Base, levels = c("T", "C","A","G"))) %>%
                add_column(position = gene_coord_min():gene_coord_max()) -> seq_df
        } else {
            tibble(Base = strsplit(Seq,split= "")[[1]]) %>%
                mutate(Base = factor(Base, levels = c("T", "C","A","G"))) %>%
                add_column(position = gene_coord_max():gene_coord_min()) -> seq_df
        }
        if(length(rv$motifs) > 0) {
            for(motif in rv$motifs) {
                motif %>%
                  gsub("U", "T", .) %>%
                  gsub("R", "[AG]", .) %>%
                  gsub("Y","[CT]", .) %>%
                  gsub("M","[AC]", .) %>%
                  gsub("K","[GT]", .) %>%
                  gsub("S","[CG]", .) %>%
                  gsub("W","[AT]", .) %>%
                  gsub("H","[ACT]", .) %>%
                  gsub("B","[CGT]", .) %>%
                  gsub("V","[ACG]", .) %>%
                  gsub("D","[AGT]", .) %>%
                  gsub("N","[ACGT]", .) %>%
                  gsub("-",",", .) -> expanded_motif # last (-) should only occur if flexible numbers of bases used
              if(!grepl("\\{", motif)) {
                motif_indices <- gregexpr(paste0("(?=",expanded_motif,")"), Seq, perl=T)[[1]]
                seq_df %>%
                    mutate(!!motif := row_number() %in% motif_indices) %>%
                    mutate(!!paste0(motif, "_length") := if_else(get(motif) == T, rv$motif_lengths[which(rv$motifs == motif)], NA_integer_)) -> seq_df
              } else {
                motif_indices_all <- tibble(index = gregexpr(paste0("(?=",expanded_motif,")"), Seq, perl=T)[[1]])
                motif_indices_nonOverlap <- gregexpr(expanded_motif, Seq)[[1]] # necessary because notation above returns all lengths as 0
                motif_indices_all %>%
                  mutate(length = attributes(motif_indices_nonOverlap)$match.length[match(index, motif_indices_nonOverlap)]) -> motif_indices_all
                for(i in motif_indices_all$index[is.na(motif_indices_all$length)]) {
                  motif_indices_all$length[motif_indices_all$index == i] <- attributes(regexpr(expanded_motif, substr(Seq, i, nchar(Seq))))$match.length
                }
                seq_df %>%
                  mutate(!!motif := row_number() %in% motif_indices_all$index) %>%
                  mutate(!!paste0(motif, "_length") := motif_indices_all$length[match(row_number(), motif_indices_all$index)]) -> seq_df
              }
            }
        }
        return(seq_df)
    })
    
    output$FDRthresh <- renderUI({
      if(input$FilterFDR == T) {
        sliderInput("FDRthreshold", "FDR <=", min = 0, max = 1, value = 1, step= 0.005)
      }
    })
    
    observeEvent(input$FDRthreshold, {
      if(input$FDRthreshold != rv$FDRthreshold) {
        rv$FDRthreshold  <- input$FDRthreshold
      }
    })
    
    output$Scorethresh <- renderUI({
      if(input$FilterFDR == T) {
        sliderInput("Scorethreshold", "Score >=", min = 0, max = 30, value = 0, step= 1)
      }
    })
    
    observeEvent(input$Scorethreshold, {
      if(input$Scorethreshold != rv$ScoreThreshold) {
        rv$ScoreThreshold  <- input$Scorethreshold
      }
    })
    
    output$ShowAll <- renderUI({
      if(input$FilterFDR == T) {
        checkboxInput("ShowFiltered", "Display excluded sites", value = F)
      }
    })
    
    observeEvent(input$ShowFiltered, {
      if(input$ShowFiltered != rv$display_all) {
        rv$display_all  <- input$ShowFiltered
      }
    })
    
    output$nReplicates <- renderUI({
      if(input$FilterFDR == T) {
        sliderInput("replicates", "Number of replicates (if merged) that must meet threshold:", value = 2, min =1, max = 5)
      }
    })
    
    observeEvent(input$replicates, {
      if(input$replicates != rv$replicate_threshold) {
        rv$replicate_threshold  <- input$replicates
      }
    })
    
    CLIPdata_chrom <- reactive({ 
        if(nchar(rv$chromosome) > 0) {
            readRDS(paste0(rv$CLIPpath, "merged_CLIP_chr", rv$chromosome, ".Rds")) 
        }
    })
    
    Clusterdata_chrom <- reactive({ 
      if(nchar(rv$chromosome) > 0) {
        readRDS(paste0("CLIP_data/merged_clusters_chr", rv$chromosome, ".Rds")) 
      }
    })
    
    all_samples <- reactive({
      if(input$mergeReplicates == T) {
        input$CLIPdatasets
      } else {
        CLIP_datasets$Sample[CLIP_datasets$dataset %in% input$CLIPdatasets]
      }
    })
    
    CLIPdata_subset <- reactive({ 
        CLIPdata_chrom() %>%
            filter(position >= coord_min() & position <= coord_max()) %>%
            filter(strand == rv$strand) %>%
            left_join(select(CLIP_datasets, Sample, dataset, replicates), by = "Sample") %>%
            filter(dataset %in% input$CLIPdatasets) %>%
            mutate(score_filter = crosslinks) -> CLIPdata_temp
        if(input$mergeReplicates == T){
          CLIPdata_temp %>%
            arrange(FDR) %>%
            group_by(chromosome, position, strand, dataset) %>%
            mutate(replicate_threshold = if_else(replicates == T, as.numeric(rv$replicate_threshold), 1)) %>%
            summarise(crosslinks = sum(crosslinks), FDR = FDR[replicate_threshold[1]], score_filter = score_filter[order(score_filter, decreasing = T)][replicate_threshold[1]], .groups = "drop") %>%
            arrange(position) %>%
            rename(Sample = dataset) -> CLIPdata_temp
        }
        CLIPdata_temp %>%
          mutate(hit_start = position-0.5, hit_end = position+0.5) -> CLIPdata_temp
        if(input$FilterFDR == T) {
          CLIPdata_temp %>% 
            mutate(display = FDR <= rv$FDRthreshold & score_filter >= rv$ScoreThreshold) %>%
            mutate(display = if_else(is.na(display), F, display)) -> CLIPdata_temp 
          if(rv$display_all == F) {
            CLIPdata_temp %>%
              filter(display == T) -> CLIPdata_temp
          }
        } else {
          CLIPdata_temp %>%
            add_column(display = T) -> CLIPdata_temp
        }
        for(d in all_samples()[!all_samples() %in% CLIPdata_temp$Sample]) {
            CLIPdata_temp %>%
                add_row(crosslinks = 0, Sample = d) -> CLIPdata_temp
        }
        CLIPdata_temp %>%
            mutate(CLIP_offset = as.numeric(as.character(factor(Sample, levels = all_samples(), labels = seq(0, by = 3, length.out = length(all_samples())))))) -> CLIPdata_temp
        if(input$CLIPscale == T) {
            CLIPdata_temp %>%
                mutate(CLIP_y = 2.5*crosslinks/max(5,ceiling(max(crosslinks))))
        } else {
            CLIPdata_temp %>%
                group_by(Sample) %>%
                mutate(CLIP_y = 2.5*crosslinks/max(5,ceiling(max(crosslinks)))) %>%
                ungroup()
        }
    })
    
    CLIPdata_zoom <- reactive({
        CLIPdata_chrom() %>%
            filter(position >= rv$zoom_min & position <= rv$zoom_max) %>%
            filter(strand == rv$strand) %>%
            left_join(select(CLIP_datasets, Sample, dataset, replicates), by = "Sample") %>%
            filter(dataset %in% input$CLIPdatasets) %>%
            mutate(score_filter = crosslinks) -> CLIPdata_temp
      if(input$mergeReplicates == T){
        CLIPdata_temp %>%
          arrange(FDR) %>%
          group_by(chromosome, position, strand, dataset) %>%
          mutate(replicate_threshold = if_else(replicates == T, as.numeric(rv$replicate_threshold), 1)) %>%
          summarise(crosslinks = sum(crosslinks), FDR = FDR[replicate_threshold[1]], score_filter = score_filter[order(score_filter, decreasing = T)][replicate_threshold[1]], .groups = "drop") %>%
          arrange(position) %>%
          rename(Sample = dataset) -> CLIPdata_temp
      }
      CLIPdata_temp %>%
        mutate(hit_start = position-0.5, hit_end = position+0.5) -> CLIPdata_temp
      if(input$FilterFDR == T) {
        CLIPdata_temp %>% 
          mutate(display = FDR <= rv$FDRthreshold & score_filter >= rv$ScoreThreshold) %>%
          mutate(display = if_else(is.na(display), F, display)) -> CLIPdata_temp 
        if(rv$display_all == F) {
          CLIPdata_temp %>%
            filter(display == T) -> CLIPdata_temp
        }
      } else {
        CLIPdata_temp %>%
          add_column(display = T) -> CLIPdata_temp
      }
      for(d in all_samples()[!all_samples() %in% CLIPdata_temp$Sample]) {
        CLIPdata_temp %>%
          add_row(crosslinks = 0, Sample = d) -> CLIPdata_temp
      }
      CLIPdata_temp %>%
        mutate(CLIP_offset = as.numeric(as.character(factor(Sample, levels = all_samples(), labels = seq(0, by = 3, length.out = length(all_samples())))))) -> CLIPdata_temp
      if(input$CLIPscale == T) {
            CLIPdata_temp %>%
                mutate(CLIP_y = 2.5*crosslinks/max(5,ceiling(max(crosslinks))))
        } else {
            CLIPdata_temp %>%
                group_by(Sample) %>%
                mutate(CLIP_y = 2.5*crosslinks/max(5,ceiling(max(crosslinks)))) %>%
                ungroup()
        }
    })
    
    Clusterdata_subset <- reactive({ 
      if(input$mergeReplicates == T) {
        Clusterdata_chrom() %>%
          filter(start >= coord_min() & end <= coord_max()) %>%
          filter(strand == rv$strand) %>%
          filter(dataset %in% all_samples()) %>%
          mutate(CLIP_offset = CLIPdata_subset()$CLIP_offset[match(dataset, CLIPdata_subset()$Sample)])
      } else {
        Clusterdata_chrom() %>%
          mutate(dataset = sample) %>%
          filter(start >= coord_min() & end <= coord_max()) %>%
          filter(strand == rv$strand) %>%
          filter(dataset %in% all_samples()) %>%
          mutate(CLIP_offset = CLIPdata_subset()$CLIP_offset[match(dataset, CLIPdata_subset()$Sample)])
      }
    })
    Clusterdata_zoom <- reactive({ 
      if(input$mergeReplicates == T) {
        Clusterdata_chrom() %>%
          filter(start >= rv$zoom_min & end <= rv$zoom_max) %>%
          filter(strand == rv$strand) %>%
          filter(dataset %in% all_samples()) %>%
          mutate(CLIP_offset = CLIPdata_subset()$CLIP_offset[match(dataset, CLIPdata_subset()$Sample)])
      } else {
        Clusterdata_chrom() %>%
          mutate(dataset = sample) %>%
          filter(start >= rv$zoom_min & end <= rv$zoom_max) %>%
          filter(strand == rv$strand) %>%
          filter(dataset %in% all_samples()) %>%
          mutate(CLIP_offset = CLIPdata_subset()$CLIP_offset[match(dataset, CLIPdata_subset()$Sample)])
      }
    })
    
    
    CLIPdata_labels <- reactive({
        CLIPdata_subset() %>% 
            group_by(Sample, CLIP_offset) %>% 
            summarise(max_all = ceiling(max(crosslinks)), .groups = "drop") -> labels_temp
        CLIPdata_zoom() %>%
            group_by(Sample, CLIP_offset) %>% 
            summarise(max_zoom = ceiling(max(crosslinks)), .groups = "drop") %>%
            full_join(labels_temp) %>%
            mutate(label_y = CLIP_offset + 1.2) %>%
            mutate(tick_y = CLIP_offset + 2.5) -> labels_temp
        if(input$CLIPscale == T) {
            labels_temp %>%
                add_column(max_zoom_scaled = max(5,labels_temp$max_zoom)) %>%
                add_column(max_all_scaled = max(5,labels_temp$max_all)) -> labels_temp 
        } else {
            labels_temp %>%
                mutate(max_zoom_scaled = if_else(max_zoom < 5, 5, max_zoom)) %>%
                mutate(max_all_scaled = if_else(max_all < 5, 5, max_all)) -> labels_temp 
        }
        if(input$mergeReplicates == T) {
          labels_temp %>%
            mutate(label_plot = sub(", ", "\n", Sample))
        } else {
          labels_temp %>%
            mutate(label_plot = sub(", ", "\n", CLIP_datasets$label[match(Sample, CLIP_datasets$Sample)]))
        }
    })
    
    
    motif_positions <- reactive({
        if(length(rv$motifs) > 0){
            mp <- list()
            for(i in 1:length(rv$motifs)) {
                motif <- rv$motifs[i]
                complete_seq() %>%
                    filter(get(motif) == T) %>%
                    add_column(motif_offset = i) -> mp[[motif]]
                if(rv$strand == "+") {
                    mp[[motif]] %>%
                        mutate(motif_start = position-0.5) %>%
                        mutate(motif_end = position-0.5+get(paste0(motif, "_length"))) -> mp[[motif]]
                } else {
                    mp[[motif]] %>%
                        mutate(motif_end = position+0.5) %>%
                        mutate(motif_start = position+0.5-get(paste0(motif, "_length"))) -> mp[[motif]]
                }
            }
            return(mp)
            
        }
    })
    
    
    plot_height <- function(){
      if(input$mergeReplicates == T) {
        if_else(nrow(transcript_coord())>0,75+(20*length(input$transcript_selection))+(60*length(input$CLIPdatasets)), 200) -> ph
      } else {
        if_else(nrow(transcript_coord())>0,75+(20*length(input$transcript_selection))+(60*nrow(CLIP_datasets[CLIP_datasets$dataset %in% input$CLIPdatasets,])), 200) -> ph
      }
        #if(zoom_length() > 0 & zoom_length() <= 100 & length(rv$motifs) ==0) { ph + 20 -> ph }
        if(length(rv$motifs) > 0) { ph + (20*length(rv$motifs)) -> ph }
        return(ph)
    }
    
    pHeight_ui <- function(){
      req(input$pHeight)
      round(input$pHeight*2.83)
    }
    pWidth_ui <- function(){
      round(input$pWidth*2.83)
    }
    
    plot_ymin <- reactive({
        -(length(input$transcript_selection) + length(rv$motifs)) -> min_y
        if(length(rv$motifs) == 0 & zoom_length() > 0 & zoom_length() <= 100) { min_y - 1 -> min_y }
        return(min_y)
    })
    
    plot_data <- function(min_x, max_x, length_x, max_y_scaled, CLIP_data_toPlot, Cluster_data_toPlot) {
      transcript_coord() %>%
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = -transcript_offset, ymax = -transcript_offset+0.8, fill = feature_type, colour = feature_type)) +
        scale_fill_manual(guide = F, breaks = c("CDS", "UTR"), limits = c("CDS", "UTR", T,F),values = c(rv$gene_colour, NA, rv$plotting_colour, "grey70"), drop = F) +
        scale_colour_manual(guide  = F, breaks = c("CDS", "UTR"), limits = c("CDS", "UTR", T,F),values = c(rv$gene_colour, rv$gene_colour, rv$plotting_colour, "grey70"), drop = F) +
        geom_segment(data = intron_coord(), aes(x = intron_start, xend = intron_end, y = -transcript_offset+0.2, yend = -transcript_offset+0.4), colour = rv$gene_colour) +
        geom_segment(data = intron_coord(), aes(x = intron_start, xend = intron_end, y = -transcript_offset+0.6, yend = -transcript_offset+0.4), colour = rv$gene_colour) +
        labs(x = paste("Chromosome", rv$chromosome), y = NULL, title = paste0(transcript_coord()$gene_name[1], "; ", transcript_coord()$gene_id[1])) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(hjust = 0.65), plot.title = element_text(size = 15, face = "bold"), plot.margin = unit(c(1, 4, 1, 1), 'lines')) +
        geom_rect(data = CLIP_data_toPlot, aes(xmin = hit_start, xmax = hit_end, ymin = CLIP_offset, ymax = CLIP_offset+CLIP_y, fill = display, colour = display)) +
        annotate("rect", ymin = plot_ymin()-1, ymax = 0, xmin = max_x +0.005*length_x, xmax = max_x+0.5*length_x, fill = "white", colour = NA) +
        annotate("rect", ymin = plot_ymin()-1, ymax = 0, xmin = min_x-0.5*length_x, xmax = min_x-0.005*length_x, fill = "white", colour = NA) +
        annotate("segment",y = CLIPdata_labels()$CLIP_offset,yend = CLIPdata_labels()$CLIP_offset, x = min_x-0.002*length_x, xend = max_x+0.002*length_x) -> gp
      if(input$clusters == T) {
        gp + geom_rect(data = Cluster_data_toPlot, aes(xmin = start, xmax = end, ymin= CLIP_offset-0.35, ymax = CLIP_offset-0.15), fill = rv$cluster_colour) -> gp
      }
      #ggplot_build(gp) -> gbp
      #gbp$layout$panel_params[[1]]$x.major_source -> all_breaks
      #if(is.null(all_breaks)) {
      #  gbp$layout$panel_params[[1]]$x$breaks -> all_breaks
      #}
      #all_breaks[all_breaks > min_x & all_breaks < max_x] -> use_breaks
      pretty(min_x:max_x, n=3) -> use_breaks
      use_breaks[use_breaks >= min_x & use_breaks <= max_x] -> use_breaks
      if(input$rev_scale == F) { 
        if(length(rv$motifs) > 0) {
          for(motif in rv$motifs) {
            motif_positions()[[motif]] %>%
              filter(motif_start <= max_x & motif_end >= min_x) %>%
              mutate(motif_start = if_else(motif_start < min_x, as.numeric(min_x), motif_start)) %>%
              mutate(motif_end = if_else(motif_end > max_x, as.numeric(max_x), motif_end)) -> motif_data
            offset <- which(rv$motifs == motif)
            gp +
              geom_rect(data = motif_data, aes(xmin = motif_start, xmax = motif_end, ymin = -max(transcript_coord()$transcript_offset)-motif_offset, ymax = -max(transcript_coord()$transcript_offset)-(motif_offset-0.8)), colour = rv$motif_cols[offset], fill = rv$motif_cols[offset], size = 0.2) +
              annotate("text", x = min_x-0.01*length_x, y = -max(transcript_coord()$transcript_offset)-(offset-0.4), label = motif, hjust = 1, size = 3.5, colour = rv$motif_cols[offset]) -> gp
          }
        }
        if(rv$strand == "+") {
          gp +
            annotate("segment", x = min_x+0.4*length_x, xend = min_x+0.55*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.15, yend = max(CLIPdata_labels()$CLIP_offset) +3.15, arrow = arrow(length = unit(0.25,"cm"), type = "closed")) +
            annotate("text", x = min_x+0.39*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.2, label = "Transcription", hjust = 1) -> gp
        } else {
          gp +
            annotate("segment", x = min_x+0.4*length_x, xend = min_x+0.25*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.15, yend = max(CLIPdata_labels()$CLIP_offset) +3.15, arrow = arrow(length = unit(0.25,"cm"), type = "closed")) +
            annotate("text", x = min_x+0.41*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.2, label = "Transcription", hjust = 0) -> gp
        }
        gp +
          annotate("text", x = min_x-0.01*length_x, y = unique(-transcript_coord()$transcript_offset+0.4), label = unique(transcript_coord()$transcript_id), hjust = 1, size = 3.5) +
          annotate("text", x = max_x+0.03*length_x, y = mean(CLIPdata_labels()$label_y), label = "Number\nof CLIP\ncrosslinks", size = 3.5, hjust = 0) +
          annotate("segment", x = max_x+0.015*length_x, xend = max_x+0.015*length_x, y = min(CLIPdata_labels()$CLIP_offset)+0.5, yend = max(CLIPdata_labels()$CLIP_offset)+2, arrow = arrow(length = unit(0.25,"cm"), type = "closed")) +
          annotate("text", x = c(min_x+0.75*length_x, min_x+0.9*length_x), y = max(CLIPdata_labels()$CLIP_offset) +3.2, label = c("CDS", "UTR"), size = 3.5, hjust = 0) +
          annotate("rect", xmin = c(min_x+0.7*length_x, min_x+0.85*length_x), xmax = c(min_x+0.74*length_x, min_x+0.89*length_x), ymin = max(CLIPdata_labels()$CLIP_offset) +2.8, ymax = max(CLIPdata_labels()$CLIP_offset) +3.5, colour = rv$gene_colour, fill = c(rv$gene_colour, NA)) +
          geom_text(data = CLIPdata_labels(), aes(x = min_x-0.012*length_x, y = label_y, label = label_plot), hjust = 1, size = 3.5) +
          geom_text(data = CLIPdata_labels(), aes(x = min_x-0.012*length_x, y = tick_y, label = max_y_scaled), hjust = 1, size = 3) +
          geom_segment(data = CLIPdata_labels(), aes(x = min_x-0.008*length_x, xend = min_x-0.002*length_x, y = tick_y, yend = tick_y)) +
          coord_cartesian(xlim = c(min_x-0.35*length_x, max_x+0.5), ylim = c(plot_ymin()-0.6,max(CLIPdata_labels()$CLIP_offset)+3.2 ), clip = "off") +
          scale_x_continuous(expand = c(0,0), breaks = use_breaks) +
          scale_y_continuous(expand = c(0,0)) -> gp
      } else {
        if(length(rv$motifs) > 0) {
          for(motif in rv$motifs) {
            motif_positions()[[motif]] %>%
              filter(motif_start <= max_x & motif_end >= min_x) %>%
              mutate(motif_start = if_else(motif_start < min_x, as.numeric(min_x), motif_start)) %>%
              mutate(motif_end = if_else(motif_end > max_x, as.numeric(max_x), motif_end)) -> motif_data
            offset <- which(rv$motifs == motif)
            gp +
              geom_rect(data = motif_data, aes(xmin = motif_start, xmax = motif_end, ymin = -max(transcript_coord()$transcript_offset)-motif_offset, ymax = -max(transcript_coord()$transcript_offset)-(motif_offset-0.8)), colour = rv$motif_cols[offset], fill = rv$motif_cols[offset], size = 0.2) +
              annotate("text", x = max_x+0.01*length_x, y = -max(transcript_coord()$transcript_offset)-(offset-0.4), label = motif, hjust = 1, size = 3.5, colour = rv$motif_cols[offset]) -> gp
          }
        }
        if(rv$strand == "+") {
          gp +
            annotate("segment", x = min_x+0.6*length_x, xend = min_x+0.75*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.15, yend = max(CLIPdata_labels()$CLIP_offset) +3.15, arrow = arrow(length = unit(0.25,"cm"), type = "closed")) +
            annotate("text", x = min_x+0.59*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.2, label = "Transcription", hjust = 0) -> gp
        } else {
          gp +
            annotate("segment", x = min_x+0.6*length_x, xend = min_x+0.45*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.15, yend = max(CLIPdata_labels()$CLIP_offset) +3.15, arrow = arrow(length = unit(0.25,"cm"), type = "closed")) +
            annotate("text", x = min_x+0.61*length_x, y = max(CLIPdata_labels()$CLIP_offset) +3.2, label = "Transcription", hjust = 1) -> gp
        }
        gp +
          annotate("text", x = max_x+0.01*length_x, y = unique(-transcript_coord()$transcript_offset+0.4), label = unique(transcript_coord()$transcript_id), hjust = 1, size = 3.5) +
          annotate("text", x = min_x-0.03*length_x, y = mean(CLIPdata_labels()$label_y), label = "Number\nof CLIP\ncrosslinks", size = 3.5, hjust = 0) +
          annotate("segment", x = min_x-0.015*length_x, xend = min_x-0.015*length_x, y = min(CLIPdata_labels()$CLIP_offset)+0.5, yend = max(CLIPdata_labels()$CLIP_offset)+2, arrow = arrow(length = unit(0.25,"cm"), type = "closed")) +
          annotate("text", x = c(min_x+0.25*length_x, min_x+0.1*length_x), y = max(CLIPdata_labels()$CLIP_offset) +3.2, label = c("CDS", "UTR"), size = 3.5, hjust = 0) +
          annotate("rect", xmin = c(min_x+0.3*length_x, min_x+0.15*length_x), xmax = c(min_x+0.26*length_x, min_x+0.11*length_x), ymin = max(CLIPdata_labels()$CLIP_offset) +2.8, ymax = max(CLIPdata_labels()$CLIP_offset) +3.5, colour = rv$gene_colour, fill = c(rv$gene_colour, NA)) +
          geom_text(data = CLIPdata_labels(), aes(x = max_x+0.012*length_x, y = label_y, label = label_plot), hjust = 1, size = 3.5) +
          geom_text(data = CLIPdata_labels(), aes(x = max_x+0.012*length_x, y = tick_y, label = max_y_scaled), hjust = 1, size = 3) +
          geom_segment(data = CLIPdata_labels(), aes(x = max_x+0.008*length_x, xend = max_x+0.002*length_x, y = tick_y, yend = tick_y)) +
          coord_cartesian(xlim = c(max_x+0.35*length_x, min_x-0.5), ylim = c(plot_ymin()-0.6,max(CLIPdata_labels()$CLIP_offset)+3.2 ), clip = "off") +
          scale_x_reverse(expand = c(0,0), breaks = use_breaks) +
          scale_y_continuous(expand = c(0,0)) -> gp
      }
      if(length_x <= 100) {
        complete_seq() %>%
          filter(position %in% c(min_x:max_x)) -> zoom_seq
        gp +
          geom_text(data = zoom_seq, aes(x = position, y = -max(transcript_coord()$transcript_offset)-0.6, label = Base, colour = Base), size = min(4.5,250/length_x)) +
          scale_colour_manual(guide = F, breaks = c("CDS", "UTR"), limits = c("CDS", "UTR", T,F,"T", "C","A","G"), values = c(rv$gene_colour, rv$gene_colour, rv$plotting_colour, "grey70", brewer.pal(4,"Set1")), drop = F) -> gp
      }
      
      return(gp)
    }
    
    complete_plot <- reactive({
      if(nrow(transcript_coord()) > 0) {
        plot_data(min_x=coord_min(), max_x=coord_max(), length_x=coord_length(), max_y_scaled=CLIPdata_labels()$max_all_scaled, CLIP_data_toPlot=CLIPdata_subset(), Cluster_data_toPlot=Clusterdata_subset())
      } 
    })
    
    zoom_plot <- reactive({
      if(nrow(transcript_coord()) > 0 & rv$zoom_max != 0) {
        plot_data(min_x=rv$zoom_min, max_x=rv$zoom_max, length_x=zoom_length(), max_y_scaled=CLIPdata_labels()$max_zoom_scaled, CLIP_data_toPlot=CLIPdata_zoom(), Cluster_data_toPlot=Clusterdata_zoom())
      }
    })
    
    output$genePlot <- renderPlot({
        complete_plot()
    }, width= pWidth_ui, height= pHeight_ui)
    
    output$genePlotZoom <- renderPlot({
        zoom_plot()
    }, width=pWidth_ui, height= pHeight_ui)
    
    output$PDFcomplete <- downloadHandler(
      filename = function() { paste0(input$gene, "_complete_gene_plot.pdf") },
      content = function(file) {
        ggsave(file, complete_plot(), width = input$pWidth, height=input$pHeight, units = "mm")
      }
    )
    output$PNGcomplete <- downloadHandler(
      filename = function() { paste0(input$gene, "_complete_gene_plot.png") },
      content = function(file) {
        ggsave(file, complete_plot(), width = input$pWidth, height=input$pHeight, units = "mm")
      }
    )
    
    output$PDFzoom <- downloadHandler(
      filename = function() { paste0(input$gene, "_zoom_plot.pdf") },
      content = function(file) {
        ggsave(file, zoom_plot(), width = input$pWidth, height=input$pHeight, units = "mm")
      }
    )
    output$PNGzoom <- downloadHandler(
      filename = function() { paste0(input$gene, "_zoom_plot.png") },
      content = function(file) {
        ggsave(file, zoom_plot(), width = input$pWidth, height=input$pHeight, units = "mm")
      }
    )
    
    #output$tableTest <- renderTable(if(nrow(transcript_coord())>0) {CLIPdata_subset() %>% head(20)})
    
    #observeEvent(input$browser, { browser() })
}

# Run the application 
shinyApp(ui = ui, server = server)
