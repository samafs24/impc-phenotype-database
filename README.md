## Interactive Dashboard and Database for IMPC Phenotypic Data Analysis
### Objective
This repository contains the workflow, database, and interactive visualisation tools developed for the project on International Mouse Phenotyping Consortium (IMPC) data. The project integrates large-scale phenotypic data from knockout mice to explore genotype-phenotype associations and functional insights. It includes a scalable MySQL database and an interactive R Shiny dashboard for visualisation and analysis of phenotypic scores, gene-disease associations, and clustering of similar phenotypes.

The goal of the project is to create a reproducible, user-friendly platform for querying, analysing, and visualising IMPC phenotypic datasets.
### Project Structure
```
├── data/
│   ├── cleaned/            # Cleaned datasets ready for analysis
│   ├── collated/           # Combined datasets
│   ├── database/           # Files for MySQL database population
│   ├── groupings/          # Parameter groupings 
│   ├── metadata/           # Procedures, parameters, and disease metadata
│   └── raw_data/           # Original CSV files from IMPC
│
├── scripts/
│   ├── development_code/   # R and Bash scripts for cleaning, processing, and analysis
│   ├── rshiny_app/         # Shiny dashboard code and supporting scripts
│   └── dcdm_final_quarto.qmd  # Quarto document detailing full reproducible analysis
```
