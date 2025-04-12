# Loading libraries needed for DNA GC content 
library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(Biostrings)

# Increasing the file size limit for uploads
options(shiny.maxRequestSize = 100 * 1024^2)

# This function calculates GC content from a DNA sequence
# It splits the sequence into individual characters counts the number of G and C bases,and returns the GC percentage rounded to two decimal places
calc_gc_content <- function(sequence) {
  seq_chars <- strsplit(as.character(sequence), "")[[1]]
  gc_count <- sum(seq_chars %in% c("G", "C"))
  total_bases <- length(seq_chars)
  round((gc_count / total_bases) * 100, 2)
}

# This function calculates the percentages of A, T, G, and C
# It creates a frequency table for each base and returns the values as percentages of the total base count
calc_base_percentages <- function(sequence) {
  seq_chars <- strsplit(as.character(sequence), "")[[1]]
  total_bases <- length(seq_chars)
  counts <- table(factor(seq_chars, levels = c("A", "T", "G", "C")))
  percentages <- round((counts / total_bases) * 100, 2)
  return(as.numeric(percentages))
}

# Defining the User Interface (UI)
# The sidebar lets users upload a FASTA file or paste a sequence
# The main panel displays plots, tables, and sequence info
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),
  titlePanel("DNA GC Content Analyzer"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload FASTA File", accept = c(".fasta", ".fa", ".fna", ".txt")),
      
      tags$hr(),
      strong("Or paste a DNA sequence below"),
      textAreaInput("dna_seq", NULL, placeholder = "ATGCGT..."),
      
      actionButton("analyze", "Analyze"),
      tags$hr(),
      
      selectInput("seq_selector", "Select a Sequence to View", choices = NULL),
      
      verbatimTextOutput("sequence_view"),
      
      downloadButton("download_fasta", "Download Selected Sequence - GC Content"),
      br(), br(),
      
      downloadButton("download_csv", "Download Summary CSV")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Nucleotide Plot", 
                 plotOutput("nucleotide_plot"),
                 br(),
                 downloadButton("download_plot", "Download Plot")
        ),
        tabPanel("Summary Table", DTOutput("summary_table"))
      )
    )
  )
)

# This is the server function
# It defines how data is processed and displayed
server <- function(input, output, session) {
  
  # This reactive block triggers when the Analyze button is clicked
  # It processes either a FASTA file or a pasted DNA string
  dna_data <- eventReactive(input$analyze, {
    if (!is.null(input$file)) {
      # If a file is uploaded, read it using Biostrings
      sequences <- readDNAStringSet(input$file$datapath, format = "fasta")
    } else {
      # If a sequence is pasted, clean it by removing non-ATGC characters
      raw <- toupper(gsub("[^ATGC]", "", input$dna_seq))
      if (nchar(raw) == 0) return(NULL)
      sequences <- DNAStringSet(setNames(DNAString(raw), "User_Input"))
    }
    
    # For each sequence, calculate GC content, base percentages, and length
    results <- lapply(names(sequences), function(name) {
      seq <- sequences[[name]]
      percentages <- calc_base_percentages(seq)
      list(
        name = name,
        sequence = as.character(seq),
        length = nchar(as.character(seq)),
        gc = calc_gc_content(seq),
        A = percentages[1],
        T = percentages[2],
        G = percentages[3],
        C = percentages[4]
      )
    })
    
    return(results)
  })
  
  # After processing, update the dropdown menu
  # This lets us to  select a sequence to view in detail
  observeEvent(dna_data(), {
    updateSelectInput(inputId = "seq_selector", choices = sapply(dna_data(), function(x) x$name))
  })
  
  # Here showing the summary stats of the selected sequence
  output$sequence_view <- renderText({
    req(dna_data())
    req(input$seq_selector)
    
    selected <- Filter(function(x) x$name == input$seq_selector, dna_data())
    if (length(selected) == 0) return("No matching sequence found.")
    selected <- selected[[1]]
    
    paste0(
      "Sequence: ", selected$name, "\n",
      "Length: ", selected$length, " bp\n",
      "GC Content: ", selected$gc, " %\n",
      "A: ", selected$A, " % | ",
      "T: ", selected$T, " % | ",
      "G: ", selected$G, " % | ",
      "C: ", selected$C, " %"
    )
  })
  
  # Render a summary table showing key stats for all sequences
  output$summary_table <- renderDT({
    data <- dna_data()
    req(data)
    
    df <- data.frame(
      Name = sapply(data, function(x) x$name),
      Length = sapply(data, function(x) x$length),
      GC_Content = sapply(data, function(x) x$gc),
      A_Percentage = sapply(data, function(x) x$A),
      T_Percentage = sapply(data, function(x) x$T),
      G_Percentage = sapply(data, function(x) x$G),
      C_Percentage = sapply(data, function(x) x$C),
      stringsAsFactors = FALSE
    )
    
    datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })
  
  # Render a bar plot of base percentages
  # Each sequence gets a separate panel
  output$nucleotide_plot <- renderPlot({
    data <- dna_data()
    req(data)
    
    plot_data <- do.call(rbind, lapply(data, function(x) {
      data.frame(
        Sequence = x$name,
        Nucleotide = c("A", "T", "G", "C"),
        Percentage = c(x$A, x$T, x$G, x$C)
      )
    }))
    
    ggplot(plot_data, aes(x = Nucleotide, y = Percentage, fill = Nucleotide)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ Sequence, scales = "free_y") +
      theme_minimal() +
      ggtitle("Nucleotide Composition (%) per Sequence")
  })
  
  # Allow users to download the summary table as a CSV file
  output$download_csv <- downloadHandler(
    filename = function() {
      paste("gc_summary_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      data <- dna_data()
      req(data)
      
      df <- data.frame(
        Name = sapply(data, function(x) x$name),
        Length = sapply(data, function(x) x$length),
        GC_Content = sapply(data, function(x) x$gc),
        A_Percentage = sapply(data, function(x) x$A),
        T_Percentage = sapply(data, function(x) x$T),
        G_Percentage = sapply(data, function(x) x$G),
        C_Percentage = sapply(data, function(x) x$C),
        stringsAsFactors = FALSE
      )
      
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # Download a selected sequence in FASTA format
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0(input$seq_selector, "_sequence.fasta")
    },
    content = function(file) {
      selected <- dna_data()[[input$seq_selector]]
      seq <- selected$sequence
      fasta_text <- paste0(">", input$seq_selector, "\n", paste(strwrap(seq, 80), collapse = "\n"))
      writeLines(fasta_text, file)
    }
  )
  
  # Save the plot as an image file in PNG format
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("nucleotide_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      data <- dna_data()
      req(data)
      
      plot_data <- do.call(rbind, lapply(data, function(x) {
        data.frame(
          Sequence = x$name,
          Nucleotide = c("A", "T", "G", "C"),
          Percentage = c(x$A, x$T, x$G, x$C)
        )
      }))
      
      p <- ggplot(plot_data, aes(x = Nucleotide, y = Percentage, fill = Nucleotide)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ Sequence, scales = "free_y") +
        theme_minimal() +
        ggtitle("Nucleotide Composition (%) per Sequence")
      
      ggsave(file, p, width = 10, height = 6)
    }
  )
}

# Starting the app
shinyApp(ui = ui, server = server)
