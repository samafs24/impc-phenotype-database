# Install packages
# install.packages(c("shiny", "tidyverse", "plotly", "cluster", "reshape2", "DBI", "RMySQL", "factoextra", "umap", "dendextend", "ggdendro"))

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

# Define UI
ui <- fluidPage(
  titlePanel("Statistical Analysis and Visualisation of Knockout Mice"),
  
  sidebarLayout(
    sidebarPanel(
      
      # Dropdown menu for figure 1 
      conditionalPanel(
        condition = "input.tabs == 'mouse_genotype_tab'",
        selectInput("genotype_mouse", 
                    "Select Knockout Mouse:", 
                    choices = NULL, selected = NULL),  # Dynamically populated
        sliderInput("genotype_threshold", 
                    "Significance Threshold (p-value):", 
                    min = 0, max = 1, value = 0.05, step = 0.01),
        selectInput("genotype_plot_type", 
                    "Select Plot Type:", 
                    choices = c("All Phenotypes", "Top 25 Phenotypes"), selected = "All Phenotypes"),
        selectInput("genotype_mouse_strain", 
                    "Select Mouse Strain:", 
                    choices = NULL, selected = "All"), # Dynamically populated
        selectInput("genotype_life_stage", 
                    "Select Mouse Life Stage:", 
                    choices = NULL, selected = "All")  # Dynamically populated
      ),
      
      # Dropdown menu for figure 2
      conditionalPanel(
        condition = "input.tabs == 'mouse_phenotype_tab'",
        sliderInput("phenotype_threshold", 
                    "Significance Threshold (p-value):", 
                    min = 0, max = 1, value = 0.05, step = 0.01),
        selectInput("phenotype_plot_type", 
                    "Select Plot Type:", 
                    choices = c("All Genotypes", "Top 25 Genotypes"), selected = "All Genotypes"),
        selectInput("procedure", 
                    "Select Procedure:", 
                    choices = NULL, selected = NULL),  # Dynamically populated
        selectInput("phenotype", 
                    "Select Phenotype:", 
                    choices = NULL, selected = NULL),  # Dynamically populated
        textOutput("phenotype_explanation")
      ),
      
      # Dropdown menu for figure 3
      conditionalPanel(
        condition = "input.tabs == 'clustering_tab'",
        # Propose multiple clustering methods 
        selectInput("cluster_method", "Clustering Method:", choices = c("Hierarchical", "PCA", "UMAP"), selected = "Hierarchical"),
        # Only show number of clusters if PCA or UMAP is chosen
        conditionalPanel(
          condition = "input.cluster_method == 'PCA' || input.cluster_method == 'UMAP'",
          numericInput("num_clusters", "Number of Clusters (K-Means):",
                       value = 3, min = 2, max = 50, step = 1)
        ),
        selectInput("gene_subset", "Subset of Genes:",
                    choices = c("All genes", 
                                "Genes with Significant Phenotypes (p<0.05)"),
                    selected = "All genes"),
        )
    ),

    mainPanel(
      tabsetPanel(
        id = "tabs",
        
        # Tab for figure 1
        tabPanel(
          "Phenotype Scores for Selected Knockout Mouse",
          value = "mouse_genotype_tab",
          uiOutput("genotype_plot_container"),  # Placeholder for dynamically rendered plot output
          downloadButton("download_mouse_data", "Download Mouse Data")  # Button for downloading mouse data
        ),
        
        # Tab for figure 2 
        tabPanel(
          "Knockout Gene Scores for Selected Phenotype",
          value = "mouse_phenotype_tab",
          uiOutput("phenotype_plot_container"),  # Placeholder for dynamically rendered plot output
          downloadButton("download_phenotype_data", "Download Phenotype Data")  # Button for downloading phenotype data
        ),
        
        # Tab for figure 3
        tabPanel(
          "Gene Clusters",
          value = "clustering_tab",
          uiOutput("clustering_plot_container"),  # Placeholder for dynamically rendered plot output
          downloadButton("download_cluster_data", "Download Cluster Data")  # Button for downloading cluster data
        )
      )
    )
  )
)
