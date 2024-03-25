#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# install.packages("devtools")
# library(devtools)
# install_github("rianoble/acornfinder")
# library(acornfinder)



#################################################

# Imports

# Probe
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("rprimer")

require("rprimer")
#library(Biostrings)
require(tcltk)

# Data processing
require("DT")
require("dplyr")
require("tidyverse")
require("stringi")
require("stringr")
require("mosaic")
require("purrr")


#graphing
require("ggplot2")
require("hexbin")
require("patchwork")
require("plotly")


# Bioinformatics
require("BiocManager")
require("biomaRt")
require("spgs")
require("primer3")


# Deployment
#require("shinydashboard")
#require("shiny")
#require("shinycssloaders")
#require("shinyWidgets")

#source("functions.R")

#options(repos = BiocManager::repositories())
#options(scipen = 999)

# BACKING FUNCTIONS/////////////

get_strong1 <- function(x, type){
  temp <- ""
  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)

  if (target == "A") {temp <- "G"} else
    if (target == "G") {temp <- "A"} else
      if (target == "C") {temp <- "T"} else
        if (target == "T") {temp <- "C"}

  if(type){
    substring(x, 3, 3) <- temp
  }else
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp

  return(x)
}

get_strong2 <- function(x, type){
  temp <- ""

  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)

  if (target == "T") {
    temp <- "T"
    if(type){
      substring(x, 3, 3) <- temp
    }else
      substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)}
  else
    return("N")
}

get_medium1 <- function(x, type){
  temp <- ""
  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)

  if (target == "A") {temp <- "A"} else
    if (target == "G") {temp <- "G"} else
      if (target == "C") {temp <- "C"} else
        return("N")

  if(type){
    substring(x, 3, 3) <- temp
  }else
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp

  return(x)
}


get_weak1 <- function(x, type){
  temp <- ""
  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)

  if (target == "C") {temp <- "A"} else
    if (target == "A") {temp <- "C"} else
      if (target == "G") {temp <- "T"} else
        if (target == "T") {temp <- "G"}
  if(type){
    substring(x, 3, 3) <- temp
  }else
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
}


## Warngling SNP list to individual rows
list_seq <- function(snp) {
  first_position <- unlist(gregexpr('/', snp))[1]
  number_slash <- str_count(snp, "/")
  block <- str_sub(snp, first_position -1 ,
                   first_position - 2 + number_slash*2 + 1)
  block <- gsub("/", "", block)
  block <- strsplit(block, "")

  for (i in block) {
    k <- paste(str_sub(snp, 1, first_position - 2),
               i,
               str_sub(snp, first_position - 2 + number_slash * 2 + 2, str_length(snp) ),
               sep = "")

  }
  k = gsub("%", "", k)
  k = k[!grepl("W", k)]
  return(k)
}

## Clean out the nested table for the algorithum
get_list <- function(i, j){
  k <- str_flatten(nested_tables[[i]][[j]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  return(k)
}


# Get all endpoints and give parents and children
get_endpoints <- function(lst, current_name = "", parent_names = character()) {
  endpoints <- list()

  if (is.list(lst)) {
    if (length(lst) > 0) {
      for (i in seq_along(lst)) {
        nested_name <- names(lst)[i]
        nested_value <- lst[[i]]

        if (is.list(nested_value)) {
          nested_endpoints <- get_endpoints(
            nested_value,
            paste(current_name, nested_name, sep = "/"),
            c(parent_names, current_name)
          )
          endpoints <- c(endpoints, nested_endpoints)
        } else {
          endpoint <- list(
            endpoint = nested_name,
            parents = c(parent_names, current_name)
          )
          endpoints <- c(endpoints, list(endpoint))
        }
      }
    } else {
      endpoint <- list(
        endpoint = current_name,
        parents = parent_names
      )
      endpoints <- c(endpoints, list(endpoint))
    }
  } else {
    endpoint <- list(
      endpoint = current_name,
      parents = parent_names
    )
    endpoints <- c(endpoints, list(endpoint))
  }
  return(endpoints)
}


# A function tha cleans
clean_endpoints <- function(endpoints){
  for (i in 1:length(endpoints)){
    endpoints[[i]]$parents <- endpoints[[i]]$parents[-1]
    for (j in 1:length(endpoints[[i]]$parents)){

      split_string <- strsplit(endpoints[[i]]$parents[j], "/")[[1]]
      desired_item <- split_string[length(split_string)]
      endpoints[[i]]$parents[j] <- desired_item
    }
  }
  return(endpoints)
}


# find the bad nodes
compute_bad_nodes <- function(endpoints, threshold){
  blacklist <- list()

  for (i in 1:length(endpoints)){
    result = 0
    for (j in 1:length(endpoints[[i]]$parents)){
      result = result + (calculate_dimer(endpoints[[i]]$endpoint,
                                         endpoints[[i]]$parents[j])$temp > threshold)
      # print(calculate_dimer(endpoints[[i]]$endpoint, endpoints[[i]]$parents[j])$temp)
    }
    blacklist <- c(blacklist, result)
  }

  bad_nodes = endpoints[blacklist == 1]
  return(bad_nodes)
}


### Remove unqualafied nods from the OG df
remove_list <- function(lst, path) {
  if (length(path) == 1) {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- NULL
    }
  } else {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- remove_list(lst[[path[[1]]]], path[-1])
      if (is.list(lst[[path[[1]]]]) && length(lst[[path[[1]]]]) == 0 && !any(names(lst[[path[[1]]]]))) {
        lst[[path[[1]]]] <- NULL
      }
    }
  }
  lst
}

## remove the trunk
remove_empty_lists <- function(lst) {
  if (is.list(lst)) {
    lst <- lapply(lst, remove_empty_lists)
    lst <- lst[lengths(lst) > 0]
  }
  lst
}


### Remove based on index
Iterate_remove <- function(level3,bad_nodes){
  for (i in 1:length(bad_nodes)){
    level3 = remove_list(level3, c(bad_nodes[[i]]$parents, bad_nodes[[i]]$endpoint))
  }
  return(level3)
}

# Gather incoming list
incoming_list <- function(arranged_list){
  level4 <- list()
  for (item in arranged_list) {
    # Create a sublist with the name as the item
    sublist <- list(1)
    names(sublist) <- item
    level4[item] <- sublist
  }
  return(level4)
}

# replacing the end nodes so we have a new list at the end
replace_end_nodes <- function(lst, replace_lst) {
  if (is.list(lst)) {
    if (length(lst) == 0) {
      return(replace_lst)
    } else {
      return(lapply(lst, function(x) replace_end_nodes(x, replace_lst)))
    }
  } else {
    return(replace_lst)
  }
}


## get the primers that is around the SNP and apply mismatches rules
extract_substrings <- function(string, center, start_distance , end_distance) {
  # Empty lists to store substrings
  substrings_left <- list()
  substrings_right <- list()

  for (item in string) {

    # right flanking
    for (distance in start_distance:end_distance) {
      sub <- substr(item, center+1, center+1 + distance)
      substrings_right <- c(substrings_right,
                            toupper(reverseComplement(get_strong1(sub,1))),
                            toupper(reverseComplement(get_strong2(sub,1))),
                            toupper(reverseComplement(get_medium1(sub,1))),
                            toupper(reverseComplement(get_weak1(sub,1))))
    }

    # Left flanking
    for (distance in start_distance:end_distance) {
      # print(substr(string, center -distance, center))
      sub <- substr(item, center - distance, center +1)
      substrings_left <- c(substrings_left,
                           get_strong1(sub,0),
                           get_strong2(sub,0),
                           get_medium1(sub,0),
                           get_weak1(sub,0))
    }

    start_distance = start_distance + 1
    end_distance = end_distance + 1

  }

  # Return the extracted substrings
  return(list(left = substrings_left[!substrings_left %in% "N"],
              right = substrings_right[!substrings_right %in% "N"]))
}

## Get the primers that are far away from the SNP
extract_substrings_far <- function(string,
                                   center,
                                   start_distance ,
                                   end_distance,
                                   far,
                                   shift) {
  # Empty lists to store substrings
  substrings_left <- list()
  substrings_right <- list()

  for (i in seq(far, far + shift)){
    # to Right
    for (distance in start_distance:end_distance) {
      # print(substr(string, far + center, far + center + distance))
      sub = substr(string, i + center, i + center + distance)
      substrings_right <- c(substrings_right, toupper(reverseComplement(sub)))
    }

    # Left flanking
    for (distance in start_distance:end_distance) {
      # print(substr(string, center - distance - far, center - far))
      substrings_left <- c(substrings_left,
                           substr(string,
                                  center - distance - i,
                                  center - i))}
  }

  # Return the extracted substrings
  return(list(right = substrings_right, left = substrings_left))
}


## The one that produce all the primers
all_text_wrangling <- function(snp_wrangled,
                               start_distance,
                               end_distance,
                               center,
                               far,
                               shift){

  ## extrac the candidate from the left side (upstream) of the SNP
  ## extrac the candidate from the right side (downstream) of the SNP
  ## We are only getting the primers that are closed to the SNP for now
  grouped_sequences <- snp_wrangled %>%
    group_by(snpID) %>%
    summarize(sequence_list = list(sequence)) %>%
    mutate(substrings = map(sequence_list, ~extract_substrings(.x,
                                                               center,
                                                               start_distance,
                                                               end_distance))) %>%
    unnest(substrings)


  grouped_sequences_far <- snp_wrangled %>%
    group_by(snpID) %>%
    slice(1:1)%>%
    ungroup() %>%
    mutate(substrings = map(sequence,
                            ~extract_substrings_far(.x,
                                                    center,
                                                    start_distance,
                                                    end_distance,
                                                    far,
                                                    shift))) %>% unnest(substrings)

  grouped_sequences$faraway <- grouped_sequences_far$substrings
  grouped_sequences <-  grouped_sequences[, -2]


  grouped_sequences$direction <- duplicated(grouped_sequences[[1]])

  grouped_sequences <- grouped_sequences %>%
    mutate(direction = ifelse(direction == FALSE, "LEFT", "RIGHT"))

  return(grouped_sequences)
}


# Apply all the filters before multiplexing
stage1_filter <- function(df,
                          desired_tm,
                          diff,
                          Homodimer,
                          hairpin){
  df

  # This is the soft filter. We first make sure there is leftafter after the filtering. If not, we keep the best option
  for (i in 1:length(df[[2]])){

    # Homodimer
    k = df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_homodimer)[2,]) < Homodimer]
    if (length(k) > 5) {
      df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_homodimer)[2,]) < Homodimer]
    }else{
      print(paste("Homodimer - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[2]][[i]], calculate_homodimer)
      differences <- abs(unlist(calculated_values[2,]) - Homodimer)
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[2]][[i]] <- df[[2]][[i]][min_diff_indices]
    }

    # Hairpin
    k = df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_hairpin)[2,]) < hairpin]
    if (length(k) > 5) {
      df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_hairpin)[2,]) < hairpin]
    }else{
      print(paste("Hairpin - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[2]][[i]], calculate_hairpin)
      differences <- abs(unlist(calculated_values[2,]) - hairpin)
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[2]][[i]] <- df[[2]][[i]][min_diff_indices]
    }

    # Filter Tm above target
    k = df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_tm)) < desired_tm + diff]
    if (length(k) > 5) {
      df[[2]][[i]] <- k
    }else{
      print(paste("Tm_above - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[2]][[i]], calculate_tm)
      differences <- abs(unlist(calculated_values) - (desired_tm + diff) )
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[2]][[i]] <- df[[2]][[i]][min_diff_indices]
    }

    # df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_tm)) > desired_tm - diff]
    # Filter Tm below target
    k = df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_tm)) > desired_tm - diff]
    if (length(k) > 5) {
      df[[2]][[i]] <- k
    }else{
      print(paste("TM below - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[2]][[i]], calculate_tm)
      differences <- abs(unlist(calculated_values) - (desired_tm - diff) )
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[2]][[i]] <- df[[2]][[i]][min_diff_indices]
    }

  }
  df

  for (i in 1:length(df[[3]])){
    if (length(df[[3]][[i]]) != 0){
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_tm)) > desired_tm - diff]
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_hairpin)[2,]) < hairpin]
    }
  }
  df

  for (i in 1:length(df[[3]])){
    if (length(df[[3]][[i]]) != 0){
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_homodimer)[2,]) < Homodimer]
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_tm)) < desired_tm + diff]
    }
  }

  df
  for (i in length(df[[1]]):1){
    if (length(df[[2]][[i]]) == 0){
      df <- df[-i, ]
    }
  }
  df

  return(df)
}


### select the top n primers for multiplexing
extract_top_n <- function(nested_list, n) {
  modified_list <- lapply(nested_list, function(inner_list) {
    if (length(inner_list) >= n) {
      inner_list[1:n]
    } else {
      inner_list
    }
  })
  return(modified_list)
}

# This handle what part of the tree we want to show
get_display_tree <- function(level3, keep){
  endpoints <- get_endpoints(level3)

  # Endpoints come back a little messy
  endpoints <- clean_endpoints(endpoints)


  display_tree <- list()
  for (i in 1:keep){
    display_tree <- c(display_tree, list(unlist(endpoints[[i]])))
  }

  display_tree <- data.frame(display_tree)

  first_row <- display_tree[1, ]
  display_tree <- display_tree[-1, ]

  display_tree <- rbind(display_tree, first_row)

  colnames(display_tree) <- paste0("Option ", seq(1, keep))


  return(display_tree)
}



soulofmultiplex <- function(df, Heterodimer_tm){
  list_3 <- list()

  for (i in 1:length(df[[1]])){
    list_3 <- c(list_3,
                list(unlist(df[[4]][[i]])),
                list(unlist(df[[5]][[i]])))
  }

  # Arrange the list from small to big
  arranged_list <- list_3

  # Prepare the initial list for multiplexing
  level2 <- list()
  level3 <- list()
  level4 <- list()

  level2 <- incoming_list(arranged_list[[1]])

  level3 <- replace_end_nodes(incoming_list(arranged_list[[1]]),
                              incoming_list(arranged_list[[2]])
  )

  if (length(arranged_list) != 2) {

    level3 <- replace_end_nodes(level3,
                                incoming_list(arranged_list[[3]])
    )

    str(level3)
    # arranged_list
    # Running
    print(length(arranged_list))
    for (i in 4:length(arranged_list)){

      # Start a timer
      start_time <- Sys.time()

      # Get all the end points from the tree
      endpoints <- get_endpoints(level3)

      # Endpoints come back a little messy
      endpoints <- clean_endpoints(endpoints)
      print(paste("Start with ", length(endpoints)))

      # Evalauate all the ned points to its parents
      bad_nodes <- compute_bad_nodes(endpoints, Heterodimer_tm)
      print(paste("We are removing: ", length(bad_nodes)))


      # Remove bad nodes if there are any
      if (length(bad_nodes) != 0){
        level3 <- Iterate_remove(level3,bad_nodes)
        level3 <- remove_empty_lists(level3)
      }


      # If all nodes are bad, return NULL
      if (length(endpoints) == length(bad_nodes)){
        print("All nodes are removed during the process")
        return(NULL)
      }

      print(paste("After trimming: ", length(get_endpoints(level3))))

      # Stop adding list if we are at the last level
      if (1){
        level4 <- incoming_list(arranged_list[[i]])
        print(paste("New list: ", length(level4)))

        level3 <- replace_end_nodes(level3, level4)
        print(paste("level3 + level4: ", length(get_endpoints(level3))))
      }

      # Summarize results for this level
      print(paste("How far are we: ", i))
      print(paste("Time" , round(Sys.time() - start_time, 1)))
      print("--------------------------")
    }
  }

  level5 <- get_display_tree(level3, 3)
  repeated_list <- rep(df[[1]], each = 2)
  suffix <- c("_forward", "_reverse")
  modified_list <- paste0(repeated_list, suffix[rep(1:length(suffix), length.out = length(repeated_list))])
  rownames(level5) <- modified_list

  level5 <- rbind(level5, Tm = c(round(mean(sapply(level5[[1]], calculate_tm)), 2),
                                 round(mean(sapply(level5[[2]], calculate_tm)), 2),
                                 round(mean(sapply(level5[[3]], calculate_tm)),2)))


  return(level5)

}


get_tm_for_all_primers <- function(level5) {


  level5_with_tm_result <- data.frame(matrix(NA, nrow = nrow(level5), ncol = 0))

  # Apply the 'calculate_tm' function to each column of the dataframe
  for (i in seq_along(level5)) {
    # Calculate TM for the column and round the result
    tm_results <- round(calculate_tm(level5[[i]]), 2)

    # Combine the original column with the TM results
    combined <- data.frame(level5[[i]], tm_results)

    # Set the column names for the combined columns
    original_col_name <- names(level5)[i]
    names(combined) <- c(original_col_name, paste0(original_col_name, "_tm"))

    # Bind the new combined columns to the result dataframe
    level5_with_tm_result <- cbind(level5_with_tm_result, combined)
  }

  # Remove the first column if it contains only NA values from the placeholder creation
  level5_with_tm_result <- level5_with_tm_result[, colSums(is.na(level5_with_tm_result)) < nrow(level5_with_tm_result)]
  rownames(level5_with_tm_result) <- rownames(level5)
  level5_with_tm_result <- as.matrix(level5_with_tm_result)
}

# END OF BACKING FUNCTIONS//////////////


# INPUTS FROM SIDEBAR MENU - TO BE IMPLEMENTED INTO FUNCTIONS

#    sidebarMenu(
#      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
# actionButton("run_button", "Run Analysis", icon = icon("play")),
# numericInput(inputId = "shift", label = "Max Length (bp)", value = 400),
# numericInput(inputId = "desired_tm", label = "desired_tm (°C)", value = 60),
# sliderInput("diff", "Max difference in TM", 1, 10, 5),
#      numericInput(inputId = "Heterodimer_tm", label = "Heterodimer (°C)", value = 50),
#      numericInput(inputId = "Homodimer", label = "Homodimer (°C)", value = 30),
#      numericInput(inputId = "top", label = "Top", value = 2),
#      numericInput(inputId = "hairpin", label = "Max Hairpin (°C)", value = 45),
#      div(style = "display: none", downloadButton("downloadData", "Download"))


# These are the parameters used for trouble shooting
#
# primer = "rs53576, rs1815739, rs7412, rs429358, rs6152"
# shift = 100
# desired_tm = 64
# diff = 3
# Heterodimer_tm = 50
# Homodimer <- 45
# top <- 2

#  consoleText <- reactiveVal("")

# Render the console text (from app code)

# FUNCTION - Primer Generator
mart_api <- function(primer,
                     shift){

  # We will start exploring options 800 bp away from the SNP location upstream and downstream
  center <- 800
  hairpin <- 45
  # from that distance of 800, we will search the range from 600 to 1,000. (800+200 and 800-200)
  far <- 200
  start_distance <- 15
  end_distance <- 28

  # Accessing database
  print("Execute MART API")
  snp_list <- strsplit(primer, " ")[[1]]
  upStream <- center
  downStream <- center
  snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp") # possibly establish earlier?
  snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'),
                        filters = c('snp_filter', 'upstream_flank', 'downstream_flank'),
                        checkFilters = FALSE,
                        values = list(snp_list, upStream, downStream),
                        mart = snpmart,
                        bmHeader = TRUE)

  #Create a new data frame
  snp_wrangled <- data.frame(matrix(ncol = 2, nrow = 0))


  # Add each variation as a new string into each row
  for (j in snp_sequence$`Variant name`){
    for (i in list_seq(snp_sequence$`Variant sequences`[snp_sequence$`Variant name`==j])){
      snp_wrangled[nrow(snp_wrangled) + 1,] <- c(j, i)
    }
  }

  # Rename columns and data frame
  colnames(snp_wrangled) = c("snpID", "sequence")


  ### I have a long long string. I want to get the left 18~25 charactors and
  # between 300 ~ 800 units away, I want another 18 ~ 25
  df <- all_text_wrangling(snp_wrangled,
                           start_distance,
                           end_distance,
                           center,
                           far,
                           shift)
  df

  print("Primer generated")
  return(df)
}



get_filter <- function(df, # primer
                       desired_tm,
                       diff, # max diff in tm
                       Heterodimer_tm,
                       Homodimer,
                       hairpin) {

  print("R get filter activated")
  # Applied filters before multiplexing
  df <- stage1_filter(df, desired_tm, diff, Homodimer, hairpin)
  print(df)

  print("Filtered")


  # Count how many candidates there are for each primer group
  df <- df %>%
    mutate(substrings_count = lengths(substrings),
           faraway_count = lengths(faraway)) %>%
    relocate(snpID, substrings_count, faraway_count, everything()) # Moves a block of columns

  # Display the updated nested tibble
  return(df)
}

get_multiplex <- function(df,
                          Heterodimer_tm,
                          top){

  print("Tree search")
  df
  # Keep only certain amount of candidates
  df[[4]] <- extract_top_n(df[[4]], top)
  df[[5]] <- extract_top_n(df[[5]], top)
  # Technical debt
  df <- df[!duplicated(df$snpID), ]



  df <- df %>%
    group_by(snpID) %>%
    filter(substrings_count == max(substrings_count))

  print(df)


  level5 <- soulofmultiplex(df, Heterodimer_tm)
  print(level5)


  level5_with_tm_result <- get_tm_for_all_primers(level5) # What is this fxn??

  return(level5_with_tm_result)
}

# TROUBLESHOOTING
# primer = "rs53576, rs1815739, rs7412, rs429358, rs6152"
# shift = 100
# desired_tm = 64
# diff = 3
# Heterodimer_tm = 50
# Homodimer <- 45
# top <- 2
#hairpin <- 45

findacorn <- function(primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin){
  mart_api(primer, shift)
  get_filter(df, desired_tm, diff, Heterodimer_tm, Homodimer, hairpin) # unexplained error?
  result <- get_multiplex(df, Heterodimer_tm, top)
  tkmessageBox(message = paste("The result is:", result))
}


# This one produces the true table we used
# does get_filter already do this?

# This produced the raw table that has not been filtered
#unfiltered
# does mart_api do this already?


# This produce summary of primer generations
#ProduceGenerationSummary (below)
#  masterTable()[c(1,2,3)] <- what is this?

# This produces the result of multiplexing
#function that uses get_multiplex function to produce a result


# GOAL OF FUNCTION: CREATE SIDE OUTPUT CONTAINING CHARTS AND SNPS THAT CAN BE COPIED/PASTED


# function to download master table as csv to downloads

# function of everything together
