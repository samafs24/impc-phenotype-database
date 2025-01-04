# Install necessary packages if not already installed
# install.packages(c("shiny", "plotly", "tidyverse", "cluster", "reshape2", "DBI", "RMySQL", "factoextra", "umap", "dendextend", "ggdendro"))

# Load libraries
library(shiny)
library(plotly)
library(tidyverse)
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
    # Populate knockout mouse options (only propose genes for which data is available in the dropout)
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes G 
                               JOIN Analyses A ON G.gene_accession_id = A.gene_accession_id 
                               WHERE A.p_value IS NOT NULL 
                               ORDER BY gene_symbol ASC;")
    updateSelectInput(session, "genotype_mouse", choices = gene_choices$gene_symbol)
    })
  
  # Dynamically update mouse strain and life stage options based on selected gene
  observe({
    req(input$genotype_mouse)
    
    # Populate mouse strain options (only propose strains where data is available for the selected gene)
    mouse_strains <- dbGetQuery(con, sprintf("
        SELECT DISTINCT mouse_strain FROM Analyses A 
        JOIN Genes G ON G.gene_accession_id = A.gene_accession_id 
        WHERE G.gene_symbol = '%s' AND A.p_value IS NOT NULL;", input$genotype_mouse))
    updateSelectInput(session, "genotype_mouse_strain", choices = c("All", mouse_strains$mouse_strain))
    
    # Populate life stage options (only propose life stages for which data is available for the selected gene)
    life_stages <- dbGetQuery(con, sprintf("
        SELECT DISTINCT mouse_life_stage FROM Analyses A 
        JOIN Genes G ON G.gene_accession_id = A.gene_accession_id 
        WHERE G.gene_symbol = '%s' AND A.p_value IS NOT NULL;", input$genotype_mouse))
    updateSelectInput(session, "genotype_life_stage", choices = c("All", life_stages$mouse_life_stage))
    })
  
  # Dynamically populate dropdowns for Figure 2
  # Populate procedure options (only procedures with associated data in Analyses)
  observe({procedures <- dbGetQuery(con, "
        SELECT DISTINCT PT.procedure_name FROM ProceduresTable PT
        JOIN Parameters P ON PT.procedure_id = P.procedure_id
        JOIN Analyses A ON P.parameter_id = A.parameter_id
        WHERE A.p_value IS NOT NULL ORDER BY PT.procedure_name ASC;")
    procedures$procedure_name <- str_to_title(procedures$procedure_name)
    updateSelectInput(session, "procedure", choices = procedures$procedure_name)
  })
  
  # Populate parameter options based on the selected procedure (only parameters with associated data in Analyses)
  observe({
    req(input$procedure)
    parameters <- dbGetQuery(con, sprintf("
        SELECT DISTINCT P.parameter_name FROM Parameters P
        JOIN ProceduresTable PT ON P.procedure_id = PT.procedure_id
        JOIN Analyses A ON P.parameter_id = A.parameter_id
        WHERE PT.procedure_name = '%s' AND A.p_value IS NOT NULL
        ORDER BY P.parameter_name ASC;", input$procedure))
    parameters$parameter_name <- str_to_title(parameters$parameter_name)
    updateSelectInput(session, "phenotype", choices = parameters$parameter_name)
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
        style = "overflow-x: auto; overflow-y: hidden; height: 720px;",
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
    SELECT P.parameter_name, P.parameter_id, AVG(CASE WHEN A.p_value = 0 THEN 0.000001 ELSE A.p_value END) AS avg_p_value
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
      geom_bar(stat = "identity", width = 0.8, show.legend = TRUE) +
      scale_x_discrete(labels = data$parameter_name) +  # Display only parameter_name
      scale_fill_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) + 
      labs(
        title = paste(input$genotype_plot_type, "for", input$genotype_mouse),
        x = "Phenotype",
        y = "Significance (-log2(p-value))"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 0, vjust = 1, size = 9, face = "bold"),
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
        style = "overflow-x: auto; overflow-y: hidden; height: 550px;",
        plotlyOutput("mouse_phenotype_plot", width = "5000px", height = "100%")
      )
    } else {
      plotlyOutput("mouse_phenotype_plot", width = "100%", height = "700px")
    }
  })
  
  output$mouse_phenotype_plot <- renderPlotly({
    req(input$phenotype, input$procedure)
    
    # Query to fetch p-values for the selected phenotype and procedure
    query <- sprintf("SELECT AVG(CASE WHEN A.p_value = 0 THEN 0.000001 ELSE A.p_value END) AS avg_p_value, G.gene_symbol FROM Analyses A
                     JOIN Genes G ON A.gene_accession_id = G.gene_accession_id
                     JOIN Parameters P ON A.parameter_id = P.parameter_id
                     JOIN ProceduresTable PT ON P.procedure_id = PT.procedure_id
                     WHERE P.parameter_name = '%s' AND PT.procedure_name = '%s'
                     GROUP BY G.gene_symbol ORDER BY avg_p_value ASC;", 
                     input$phenotype, input$procedure)
    
    cat(query)
    
    # Fetch data from the database
    data <- dbGetQuery(con, query)
    
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
      geom_bar(stat = "identity", width = 0.8, show.legend = TRUE) +
      scale_fill_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) +
      labs(
        title = paste(input$phenotype_plot_type, "for", input$phenotype),
        subtitle = paste("Showing genes with p-value <= ", input$phenotype_threshold),
        x = "Genotype", 
        y = "Significance (-log10(p-value))"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(size = 10)
      ) +
      geom_hline(yintercept = -log10(input$phenotype_threshold), linetype = "dashed", color = "black")
    
    ggplotly(p, tooltip = "text")
  })
  
  #Visualisation 3: Clustering of Genes based on Similarity of Parameter Score 
  # Render the UI container based on the selected plot type
  output$clustering_plot_container <- renderUI({
    if (input$cluster_method == "Hierarchical") {
      div(
        style = "overflow-x: hidden; overflow-y: auto; height: 700px;",  # Enable vertical scrolling
        plotlyOutput("clustering_tab", width = "100%", height = "3000px")
      )
    } else {
      plotlyOutput("clustering_tab", width = "100%", height = "700px")
    }
  })
  
  # Execute the query
  data <- dbGetQuery(con, "SELECT gene_accession_id, parameter_id, ROUND(AVG(p_value), 6) AS avg_rounded_pvalue
                       FROM Analyses WHERE p_value IS NOT NULL GROUP BY gene_accession_id, parameter_id
                       ORDER BY avg_rounded_pvalue ASC;")
  
  gene_symbols <- dbGetQuery(con, "Select gene_symbol, gene_accession_id FROM Genes;")
  
  # Transform data: rows = genes, columns = parameters, cell values = avg p-values
  pca_matrix <- reactive({
    data %>%
      pivot_wider(
        names_from = parameter_id,
        values_from = avg_rounded_pvalue
      ) %>%
      column_to_rownames("gene_accession_id") 
  })
  
  # Render the cluster plots
  output$clustering_tab <- renderPlotly({
    req(pca_matrix(), gene_symbols)  # Ensure the PCA matrix is available
    data_wide <- pca_matrix()  # This is the wide gene-by-parameter matrix
    
    # Filter genes where no p-value is lower than 0.05
    if (input$gene_subset == "Genes with Significant Phenotypes (p<0.05)") {
      keep_rows <- apply(data_wide, 1, function(x) all(x >= 0.05, na.rm = TRUE))
      data_wide <- data_wide[keep_rows, , drop = FALSE]
    }

    # Standardization: scale the data so all features contribute equally to analysis
    mat_scaled <- scale(data_wide)
    
    # Hierarchical clustering
    if (input$cluster_method == "Hierarchical") {
      # Perform hierarchical clustering on scaled data
      hc <- hclust(dist(mat_scaled), method = "ward.D")
      #    Convert hclust -> dendrogram -> apply "hang" to shorten tips
      #    "hang" sets how far tips hang below the rest of the dendrogram.
      #    0.1 or 0.2 often works well. Smaller => shorter vertical lines to labels.
      # Convert clustering results to a dendrogram 
      dend <- as.dendrogram(hc) %>%
        hang.dendrogram(hang = 0.2)  # tweak 0.2 as desired
      
      # Convert the dendrogram into a ggplot-compatible structure
      dend_data <- dendro_data(dend, type = "rectangle")
      
      str(dend_data$labels)
      
      
      # Add gene symbols to dendrogram labels
      dend_data$labels <- dend_data$labels %>%
        left_join(gene_symbols, by = c("label" = "gene_accession_id")) %>%  # Merge gene symbols
        mutate(label = gene_symbol)  # Replace gene IDs with gene symbols
      
      set.seed(123)  # For reproducibility

      # Generate random colors for the branches in the dendrogram
      branch_colors <- sample(colors(), nrow(dend_data$segments), replace = TRUE)  # Random branch colors
  

      # Create dendrogram plot using ggplot
      p <- ggplot() +
        # Add dendrogram branches as line segments
        geom_segment(
          data = dend_data$segments,
          aes(x = x, y = y, xend = xend, yend = yend),
          color = branch_colors,  # Assign random branch colors
          size = 0.8              # Set line thickness for branches
        ) +
        
        # Add gene labels to the dendrogram
        geom_text(
          data = dend_data$labels %>%
            mutate(y = y - max(dend_data$segments$y) * 0.2),  # Shift labels further down
          aes(x = x, y = y, label = label),
          size = 3,
          hjust = 0.5,  # Center the text horizontally
          color = "black"
        ) +
        
        # Flip coordinates to display horizontally
        coord_flip() +
        
        # Reverse y-axis to make clusters grow upward
        scale_y_reverse(expand = c(0.2, 0)) +
        theme_minimal(base_size = 10) +
        # Add plot title and axis labels
        labs(
          title = paste(input$cluster_method, "Clustering of", input$gene_subset),
          x = "Phenotypic Similarity",
          y = "Genes"
        ) +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
          axis.text.y = element_text(size = 8),                             
          panel.grid.major = element_line(color = "grey90"),                # Light grid lines
          panel.grid.minor = element_blank()                                # Remove minor grid lines
        )
      return(ggplotly(p))
    }
    
    # PCA clustering
    # PCA clustering
    else if (input$cluster_method == "PCA") {
      numeric_data <- data_wide[, sapply(data_wide, is.numeric)]  # Ensure numeric columns
      pca_result <- prcomp(numeric_data, scale. = TRUE)  # Run PCA
      
      # Create PCA data frame
      pca_data <- as.data.frame(pca_result$x[, 1:2])  # Use first two principal components
      pca_data$gene_accession_id <- rownames(data_wide)  # Add gene_accession_id
      
      # Merge gene_symbol
      pca_data <- pca_data %>%
        left_join(gene_symbols, by = "gene_accession_id")  # Add gene_symbol column
      
      # Perform k-means clustering
      km <- kmeans(pca_data[, c("PC1", "PC2")], centers = input$num_clusters)
      pca_data$cluster <- factor(km$cluster)
      
      # Create PCA plot
      plot_title <- paste("PCA Clustering of", input$gene_subset)
      p <- ggplot(pca_data, aes(
        x = PC1, y = PC2, color = cluster, 
        text = paste0("Gene: ", gene_symbol, "<br>Cluster: ", cluster))) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = rainbow(input$num_clusters)) +
        labs(
          title = plot_title,
          x = "Principal Component 1",
          y = "Principal Component 2",
          color = "Cluster"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA, size = 1.5),
          panel.grid.major = element_line(size = 0.5, linetype = "dotted", color = "gray80"),
          panel.grid.minor = element_blank()
        )
      return(ggplotly(p, tooltip = "text"))
    }
    
    # UMAP clustering
    # UMAP clustering
    else if (input$cluster_method == "UMAP") {
      umap_result <- umap(as.matrix(data_wide), n_neighbors = 15, min_dist = 0.1)
      umap_data <- data.frame(
        UMAP1 = umap_result$layout[, 1],
        UMAP2 = umap_result$layout[, 2],
        gene_accession_id = rownames(data_wide)  # Add gene_accession_id
      )
      
      # Merge gene_symbol
      umap_data <- umap_data %>%
        left_join(gene_symbols, by = "gene_accession_id")  # Add gene_symbol column
      
      # Perform k-means clustering
      km <- kmeans(umap_data[, c("UMAP1", "UMAP2")], centers = input$num_clusters)
      umap_data$cluster <- factor(km$cluster)
      
      # Create UMAP plot
      plot_title <- paste("UMAP Clustering of", input$gene_subset)
      p <- ggplot(umap_data, aes(
        x = UMAP1, y = UMAP2, color = cluster, 
        text = paste0("Gene: ", gene_symbol, "<br>Cluster: ", cluster))) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = rainbow(input$num_clusters)) +
        labs(
          title = plot_title,
          x = "UMAP 1",
          y = "UMAP 2",
          color = "Cluster"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA, size = 1.5)
        )
      return(ggplotly(p, tooltip = "text"))
    }
  })
}

