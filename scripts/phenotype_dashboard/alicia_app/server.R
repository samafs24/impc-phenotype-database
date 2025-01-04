# Install necessary packages if not already installed
# install.packages(c("shiny", "plotly", "cluster", "reshape2", "DBI", "RMySQL", "factoextra", "umap", "dendextend", "ggdendro"))

# Load libraries
library(shiny)
library(plotly)
library(cluster)
library(reshape2)
library(DBI)
library(RMySQL)
library(factoextra)
library(umap)
library(dendextend)
library(ggdendro)

# Define Server Logic
server <- function(input, output, session) {
  
  # Establish a connection to the MySQL database
  con <- dbConnect(
    RMySQL::MySQL(),
    dbname = "IMPCDb",
    host = "localhost",
    port = 3306,
    user = "root",
    password = "KCL2024!"
  )
  
  # Ensure the database connection is closed when the app stops
  onStop(function() {
    dbDisconnect(con)
  })
  
  # Dynamically populate dropdowns for Figure 1
  observe({
    # Populate knockout mouse options
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes ORDER BY gene_symbol ASC;")
    if (nrow(gene_choices) > 0) {
      updateSelectInput(session, "genotype_mouse", choices = gene_choices$gene_symbol)
    }
    
    # Populate mouse strain options
    mouse_strain <- dbGetQuery(con, "SELECT DISTINCT mouse_strain FROM Analyses;")
    if (nrow(mouse_strain) > 0) {
      updateSelectInput(session, "genotype_mouse_strain", choices = c("All", mouse_strain$mouse_strain))
    }
    
    # Populate life stage options
    life_stage <- dbGetQuery(con, "SELECT DISTINCT mouse_life_stage FROM Analyses;")
    if (nrow(life_stage) > 0) {
      updateSelectInput(session, "genotype_life_stage", choices = c("All", life_stage$mouse_life_stage))
    }
  })
  
  # Dynamically populate dropdowns for Figure 2
  observe({
    # Populate phenotype group options
    procedures <- dbGetQuery(con, "SELECT DISTINCT procedure_name FROM ProceduresTable ORDER BY procedure_name ASC;")
    procedures$procedure_name <- str_to_title(procedures$procedure_name)
    updateSelectInput(session, "procedure", choices = procedures$procedure_name)
  })
  
  observe({
    # Populate phenotype based on the selected group
    req(input$procedure)
    phenotypes <- dbGetQuery(con, sprintf("
    SELECT P.parameter_name FROM Parameters P
    JOIN ProceduresTable PT ON P.procedure_id = PT.procedure_id
    WHERE PT.procedure_name = '%s'
    ORDER BY P.parameter_name ASC;", input$procedure))
    phenotypes$parameter_name <- str_to_title(phenotypes$parameter_name)
    updateSelectInput(session, "phenotype", choices = phenotypes$parameter_name)
  })
  
  
  # Dynamically populate dropdowns for Figure 3
  observe({
    # Populate gene options for gene-specific clustering
    all_genes <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes ORDER BY gene_symbol ASC;")
    updateSelectizeInput(session, "user_genes", choices = all_genes$gene_symbol, server = TRUE)
  })
  
  ## Visualizations
  
  # Visualization 1: show phenotypes for selected genotype 
  # Render the UI container based on the selected plot type
  output$genotype_plot_container <- renderUI({
    if (input$genotype_plot_type == "All Phenotypes") {
      div(
        style = "overflow-x: auto; overflow-y: hidden; height: 700px;",
        plotlyOutput("mouse_genotype_plot", width = "5000px", height = "100%")
      )
    } else {
      plotlyOutput("mouse_genotype_plot", width = "100%", height = "700px")
    }
  })
  
  # Render the bar plot 
  output$mouse_genotype_plot <- renderPlotly({
    req(input$genotype_mouse)
    
    base_query <- sprintf("
    SELECT P.parameter_name, P.parameter_id, AVG(CASE WHEN A.p_value = 0 THEN 0.000001 ELSE A.p_value END) AS avg_p_value, COUNT(A.p_value) AS data_count
    FROM Analyses A
    JOIN Parameters P ON A.parameter_id = P.parameter_id
    WHERE A.gene_accession_id IN (
    SELECT gene_accession_id FROM Genes WHERE gene_symbol = '%s')
    AND A.p_value IS NOT NULL
                          AND P.parameter_name != 'NA'", input$genotype_mouse)
    
    # Add conditions for mouse strain and life stage if applicable
    if (input$genotype_mouse_strain != "All") {
      base_query <- paste0(base_query, " AND A.mouse_strain = '", input$genotype_mouse_strain, "'")
    }
    if (input$genotype_life_stage != "All") {
      base_query <- paste0(base_query, " AND A.mouse_life_stage = '", input$genotype_life_stage, "'")
    }
    
    # Finalize query 
    final_query <- paste0(base_query, "GROUP BY P.parameter_id, P.parameter_name
                          ORDER BY avg_p_value ASC;")
    
    cat(final_query)
    data <- dbGetQuery(con, final_query)
    
    output$no_data_message_genotype <- renderUI({
      if (is.null(data) || nrow(data) == 0) {
        tagList(
          h3("No data available for the selected parameters. Please adjust your filters.",
             style = "color: red; text-align: center;")
        )
      } else {
        NULL
      }
    })
    
    data <- data %>%
      filter(!is.na(parameter_name) & tolower(parameter_name) != "na") %>% # Remove rows with missing or invalid parameter names
      mutate(
        parameter_name = str_to_title(parameter_name), # Convert parameter names to title case
        order_var = paste0(parameter_name, "_", parameter_id), # Create a variable that differentiates parameters with the same name for ordering in the graph (based on ID)
        Threshold = ifelse(avg_p_value < input$genotype_threshold, "Significant", "Not Significant") # Determine significance
      )
    
      if (input$genotype_plot_type == "Top 25 Phenotypes") { # Subset the data for top 25 genes
        data <- data[1:min(25, nrow(data)), ]
      }
      
    # 
    p <- ggplot(data, aes(
      x = reorder(order_var, avg_p_value), # Order the data based on the order_var
      y = -log2(avg_p_value), # Transform the p-value for intuitive visualization
      fill = Threshold, 
      text = paste0("Parameter Name: ", parameter_name, "<br>P-value: ", signif(avg_p_value, digits = 3), "<br>Threshold: ", Threshold))) + # Add labels when hover 
      geom_bar(stat = "identity", width = 0.6, show.legend = TRUE) +
      scale_x_discrete(labels = data$parameter_name) +  # Display only parameter_name
      scale_fill_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) + 
      labs(
        title = paste(input$genotype_plot_type, "for", input$genotype_mouse),
        x = "Phenotype",
        y = "Significance (-log2(p-value))"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, face = "bold"),  
        axis.title.y = element_text(size = 12, face = "bold"),  
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
      ) +
      geom_hline(yintercept = -log2(input$genotype_threshold), linetype = "dashed", color = "black") # Add threshold line 
    ggplotly(p, tooltip = "text")
    })
  
  # Visualization 2: show genotypes for selected phenotype 
  # Render the UI container based on the selected plot type
  output$phenotype_plot_container <- renderUI({
    if (input$phenotype_plot_type == "All Genotypes") {
      div(
        style = "overflow-x: auto; overflow-y: hidden; height: 700px;",
        plotlyOutput("mouse_phenotype_plot", width = "5000px", height = "100%")
      )
    } else {
      plotlyOutput("mouse_phenotype_plot", width = "100%", height = "700px")
    }
  })
  
  output$mouse_phenotype_plot <- renderPlotly({
    req(input$phenotype, input$procedure)
    
    # Query to fetch p-values for the selected phenotype and procedure
    query <- sprintf("SELECT AVG(CASE WHEN A.p_value = 0 THEN 0.000001 ELSE A.p_value END) AS avg_p_value, COUNT(A.p_value) AS data_count, G.gene_symbol FROM Analyses A
                     JOIN Genes G ON A.gene_accession_id = G.gene_accession_id
                     JOIN Parameters P ON A.parameter_id = P.parameter_id
                     JOIN ProceduresTable PT ON P.procedure_id = PT.procedure_id
                     WHERE P.parameter_name = '%s' AND PT.procedure_name = '%s'
                     GROUP BY G.gene_symbol ORDER BY avg_p_value ASC;", 
                     input$phenotype, input$procedure)
    
    cat(query)
    
    # Fetch data from the database
    data <- dbGetQuery(con, query)
    
    # If no data is available, display a message
    output$no_data_message_phenotype <- renderUI({
      if (is.null(data) || nrow(data) == 0) {
        tagList(
          h3("No data available for the selected parameters. Please adjust your filters.",
             style = "color: red; text-align: center;")
        )
      } else {
        NULL
      }
    })
    
    # Group and summarize data for bar plot
    data <- data %>%
      filter(!is.na(gene_symbol) & tolower(gene_symbol) != "na") %>% # Remove rows with missing or invalid parameter names
      mutate(
        Threshold = ifelse(avg_p_value < input$phenotype_threshold, "Significant", "Not Significant") # Determine significance
      )
    
    if (input$phenotype_plot_type == "Top 25 Genotypes") { # Subset the data for top 25 genes
      data <- data[1:min(25, nrow(data)), ]
    }
    
    p <- ggplot(data, aes(x = reorder(gene_symbol, avg_p_value), y = -log10(avg_p_value), fill = Threshold, text = paste("p-value:", avg_p_value, "<br>Gene:", gene_symbol))) +
      geom_bar(stat = "identity", width = 0.6, show.legend = TRUE) +
      scale_fill_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) +
      labs(
        title = paste(input$phenotype_plot_type, "for", input$phenotype),
        subtitle = paste("Showing genes with p-value <= ", input$phenotype_threshold),
        x = "Genotype", 
        y = "Significance (-log10(p-value))"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(size = 10)
      ) +
      geom_hline(yintercept = -log10(input$phenotype_threshold), linetype = "dashed", color = "black")
    
    ggplotly(p, tooltip = "text")
  })
  
  }
