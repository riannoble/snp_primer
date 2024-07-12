# Load required libraries
library(shiny)
library(tidyr)
library(tibble)
library(rprimer)
#library(Biostrings)
library(htmltools)
#library(kableExtra)
# Data processing
library(DT)
library(dplyr)
library(tidyverse)
library(stringi)
library(stringr)
library(mosaic)
library(purrr)
#graphing
library(ggplot2)
library(hexbin)
library(patchwork)
library(plotly)
# Bioinformatics
library(BiocManager)
library(biomaRt)
library(spgs)
library(primer3)

# Source your background functions
source("ria_fxns_draft.R")

# Define your parameters
primer = "rs9462492, rs58318008, rs1421085"
shift <- 100
desired_tm <- 64
diff <- 3
Heterodimer_tm <- 50
Homodimer <- 45
top <- 2
hairpin <- 45

# Prepare initial data frame
df <- get_primer_candidates(primer, shift)
df <- get_self_filter(df)
df <- get_cross_filter(df)
df <- get_complex_filter(df)
split_dfs <- split(df, df$groupID)
