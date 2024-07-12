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


# Deployment
#require("shinydashboard")
#require("shiny")
#require("shinycssloaders")
#require("shinyWidgets")

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
  number_slash <- stringr::str_count(snp, "/")
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

        nested_current_name <- paste(current_name, nested_name, sep = "/")
        nested_parent_names <- c(parent_names, current_name)

        if (is.list(nested_value)) {
          nested_endpoints <- get_endpoints(nested_value, nested_current_name, nested_parent_names)
          endpoints <- c(endpoints, nested_endpoints)
        } else {
          endpoint <- list(endpoint = nested_name, parents = nested_parent_names)
          endpoints <- c(endpoints, list(endpoint))
        }
      }
    } else {
      endpoint <- list(endpoint = current_name, parents = parent_names)
      endpoints <- c(endpoints, list(endpoint))
    }
  } else {
    endpoint <- list(endpoint = current_name, parents = parent_names)
    endpoints <- c(endpoints, list(endpoint))
  }

  return(endpoints)
}



# A function that cleans
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
compute_bad_nodes <- function(endpoints, threshold) {
  bad_nodes <- list()
  for (i in seq_along(endpoints)) {
    endpoint <- endpoints[[i]]
    result <- sum(sapply(endpoint$parents, function(parent) calculate_dimer(endpoint$endpoint, parent)$temp > threshold))
    if (result > 0) {
      bad_nodes <- c(bad_nodes, list(endpoint))
    }
  }
  return(bad_nodes)
}



### Remove unqualafied nods from the OG df
remove_list <- function(lst, path) {
  if (length(path) == 1) {
    lst[[path[1]]] <- NULL
  } else {
    lst[[path[1]]] <- remove_list(lst[[path[1]]], path[-1])
    if (is.list(lst[[path[1]]]) && length(lst[[path[1]]]) == 0) {
      lst[[path[1]]] <- NULL
    }
  }
  return(lst)
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
Iterate_remove <- function(level3, bad_nodes) {
  for (bad_node in bad_nodes) {
    path <- c(bad_node$parents, bad_node$endpoint)
    level3 <- remove_list(level3, path)
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

  ## extract the candidate from the left side (upstream) of the SNP
  ## extract the candidate from the right side (downstream) of the SNP
  ## We are only getting the primers that are closed to the SNP for now
  grouped_sequences <- snp_wrangled %>%
    group_by(snpID, snp_character) %>%
    summarize(sequence_list = list(sequence)) %>%
    mutate(substrings = map(sequence_list, ~extract_substrings(.x,
                                                               center,
                                                               start_distance,
                                                               end_distance))) %>%
    unnest(substrings)


  grouped_sequences_far <- snp_wrangled %>%
    group_by(snpID, snp_character) %>%
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
  grouped_sequences <- grouped_sequences[, -3]

  # what is this code below for?
  grouped_sequences$direction <- duplicated(paste(grouped_sequences$snpID, grouped_sequences$snp_character))

  grouped_sequences <- grouped_sequences %>%
    mutate(direction = ifelse(direction == FALSE, "LEFT", "RIGHT"))

  return(grouped_sequences)
}

all_text_wrangling_reverse <- function(snp_wrangled,
                                       start_distance,
                                       end_distance,
                                       center,
                                       far,
                                       shift){

  ## extract the candidate from the left side (upstream) of the SNP
  ## We are only getting the primers that are close to the SNP for now
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
    slice(1:1) %>%
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

  # This is the soft filter. We first make sure there is left after after the filtering. If not, we keep the best option
  for (i in 1:length(df[[3]])){

    # Homodimer
    k = df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_homodimer)[2,]) < Homodimer]
    if (length(k) > 5) {
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_homodimer)[3,]) < Homodimer]
    }else{
      print(paste("Homodimer - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[3]][[i]], calculate_homodimer)
      differences <- abs(unlist(calculated_values[2,]) - Homodimer)
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[3]][[i]] <- df[[3]][[i]][min_diff_indices]
    }

    # Hairpin
    k = df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_hairpin)[2,]) < hairpin]
    if (length(k) > 5) {
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_hairpin)[2,]) < hairpin]
    }else{
      print(paste("Hairpin - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[3]][[i]], calculate_hairpin)
      differences <- abs(unlist(calculated_values[2,]) - hairpin)
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[3]][[i]] <- df[[3]][[i]][min_diff_indices]
    }

    # Filter Tm above target
    k = df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_tm)) < desired_tm + diff]
    if (length(k) > 5) {
      df[[3]][[i]] <- k
    }else{
      print(paste("Tm_above - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[3]][[i]], calculate_tm)
      differences <- abs(unlist(calculated_values) - (desired_tm + diff) )
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[3]][[i]] <- df[[3]][[i]][min_diff_indices]
    }

    # df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_tm)) > desired_tm - diff]
    # Filter Tm below target
    k = df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_tm)) > desired_tm - diff]
    if (length(k) > 5) {
      df[[3]][[i]] <- k
    } else {
      print(paste("TM below - Bottle neck", df[[1]][[i]]))
      calculated_values <- sapply(df[[3]][[i]], calculate_tm)
      differences <- abs(unlist(calculated_values) - (desired_tm - diff) )
      min_diff_indices <- order(differences)[1:min(5, length(differences))]
      df[[3]][[i]] <- df[[3]][[i]][min_diff_indices]
    }

  }
  df

  for (i in 1:length(df[[4]])){
    if (length(df[[4]][[i]]) != 0){
      df[[4]][[i]] <- df[[4]][[i]][unlist(sapply(df[[4]][[i]], calculate_tm)) > desired_tm - diff]
    }
    if (length(df[[4]][[i]]) != 0){
      df[[4]][[i]] <- df[[4]][[i]][unlist(sapply(df[[4]][[i]], calculate_hairpin)[2,]) < hairpin]
    }
  }
  df

  for (i in 1:length(df[[4]])){
    if (length(df[[4]][[i]]) != 0){
      df[[4]][[i]] <- df[[4]][[i]][unlist(sapply(df[[4]][[i]], calculate_homodimer)[2,]) < Homodimer]
      df[[4]][[i]] <- df[[4]][[i]][unlist(sapply(df[[4]][[i]], calculate_tm)) < desired_tm + diff]
    } else {}
  }

  df
  for (i in length(df[[1]]):1){
    if (length(df[[4]][[i]]) == 0){
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

# This handles what part of the tree we want to show
get_display_tree <- function(level3, keep){
  endpoints <- get_endpoints(level3)

  # Endpoints come back a little messy
  endpoints <- clean_endpoints(endpoints)


  display_tree <- list()
  for (i in 1:keep){
    #if (!is.na(endpoints[[i]])) {
      display_tree <- c(display_tree, list(unlist(endpoints[[i]])))
    #}
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
                list(unlist(df[[5]][[i]])),
                list(unlist(df[[6]][[i]])))
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

      # Evaluate all the ned points to its parents
      bad_nodes <- compute_bad_nodes(endpoints, Heterodimer_tm) # go back to compute_bad_nodes
      print(paste("We are removing: ", length(bad_nodes))) # why aren't we removing anything?
      # to debug, add print statements throughout the for loop for more accurate location of bug


      # Remove bad nodes if there are any
      if (length(bad_nodes) != 0){
        level3 <- Iterate_remove(level3, bad_nodes)
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


#	Check for hairpin loops, self dimers in forward primers
#	Do same for reverse primers (different order?)
#	Check melting temperatures in reverse primers to be in range of forward primer melting temps
#	Check for interactions (hairpin loops/primer-dimers)
#	Best three by closest melting temps lowest number of interactions (DON’T NARROW DOWN TO THREE)
#	Complete the whole process for the opposite strand (don’t need marking mechanism)
#	Do this for all possible snps
#	To narrow down to final three, look at closest melting temps between all, primer-dimer interactions between all, look for max reduction
#	Focus on probes later

# QUESTIONS
# -direction column in df?

# SINGULAR SNP ID

get_primer_candidates <- function(primer,
                                  shift){

  # Explore options 800 bp away from the SNP location upstream and downstream
  center <- 800
  hairpin <- 45
  # Search the range from 600 to 1,000. (800+200 and 800-200)
  far <- 200
  start_distance <- 15
  end_distance <- 28

  # Accessing database
  print("Execute MART API")
  snp_list <- strsplit(primer, " ")[[1]]
  upStream <- center
  downStream <- center

  # Accessing database
  print("Execute MART API")
  snp_list <- strsplit(primer, " ")[[1]]
  print("snp_list generated")
  upStream <- center
  print("upstream")
  downStream <- center
  print("downstream")
  snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp") # possibly establish earlier?
  snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'),
                        filters = c('snp_filter', 'upstream_flank', 'downstream_flank'),
                        checkFilters = FALSE,
                        values = list(snp_list, upStream, downStream),
                        mart = snpmart,
                        bmHeader = TRUE)

  #Create a new data frame
  snp_wrangled <- data.frame(matrix(ncol = 2, nrow = 0))


  # New dataframe contains all snp variation sequences plus a snp id column
  for (j in snp_sequence$`Variant name`){
    for (i in list_seq(snp_sequence$`Variant sequences`[snp_sequence$`Variant name`==j])){
      snp_wrangled[nrow(snp_wrangled) + 1,] <- c(j, i)
    }
  }

  # Rename columns and data frame
  colnames(snp_wrangled) = c("snpID", "sequence")

  # New dataframe contains all snp variation sequences plus a snp id column
  snp_wrangled <- mutate(snp_wrangled, snp_character = substr(snp_wrangled$sequence, 801, 801)) # 804 is the universal snp position for all snps bc of seq setup

  # I have a long long string. I want to get the left 18~25 characters and
  # between 300 ~ 800 units away, I want another 18 ~ 25
  df <- all_text_wrangling(snp_wrangled,
                           start_distance,
                           end_distance,
                           center,
                           far,
                           shift)
  df

  print("Primers generated")
  return(df)
}

get_self_filter <- function(df) {

  print("Reformatting data table")

  # Pivot longer to combine 'substrings' and 'faraway' into 'source' and 'primers'
  self_df <- df %>%
    pivot_longer(
      cols = c("substrings", "faraway"),
      names_to = "source",
      values_to = "primers"
    )

  # Unnesting to handle list-columns if necessary (depending on your data structure)
  self_df <- self_df %>%
    unnest(cols = primers)

  # Collapse list elements in 'primers' column to a single string if necessary
  self_df$primers <- sapply(self_df$primers, function(x) {
    if (is.list(x)) paste(x, collapse = ", ") else x
  })

  print("Calculating hairpins, self dimers, and temps")

  # Calculate and mutate the hairpin, selfdimer, and temp columns
  self_df <- self_df %>%
    mutate(
      hairpin = calculate_hairpin(primers)$dg * 0.004184,
      selfdimer = calculate_homodimer(primers)$dg * 0.004184,
      temp = calculate_tm(primers),
      # Replace NA values with Inf to handle them as large values in sorting
      hairpin = ifelse(is.na(hairpin), Inf, hairpin),
      selfdimer = ifelse(is.na(selfdimer), Inf, selfdimer)
    )

  print("Filtering temps")

  # Filter based on temperature range
  self_df <- self_df %>%
    filter(temp > 45 & temp < 70)

  print("Taking top ten per group")

  # Group by the specified columns, arrange, and select the top 10 rows
  self_df <- self_df %>%
    group_by(snpID, snp_character, direction, source) %>%
    arrange(desc(selfdimer), desc(hairpin)) %>%
    #filter(hairpin > -7,
    #       selfdimer > -7) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(source = recode(source,
                           "substrings" = "forward",
                           "faraway" = "reverse")) %>%
    arrange(snpID, snp_character, direction, source)
  df <- self_df
  df
  return(df)
}

# Example of how to call this function (assuming df is your data frame)
# result_df <- get_self_filter(df)


get_cross_filter <- function(df){

  print("Generating matches")
  forward_df <- filter(df, source == "forward")
  reverse_df <- filter(df, source == "reverse")

  ## THIS IS WHERE DATA RUNS INTO A BOTTLENECK!
  combined_df <- forward_df %>%
    inner_join(reverse_df, by = c("snpID", "snp_character", "direction"), relationship = "many-to-many") %>%
    arrange(snpID, snp_character, direction) %>%
    select(-source.x, -source.y) %>%
    rename(forward_primer = primers.x,
           forward_hairpin = hairpin.x,
           forward_selfdimer = selfdimer.x,
           forward_temp = temp.x,
           reverse_primer = primers.y,
           reverse_hairpin = hairpin.y,
           reverse_selfdimer = selfdimer.y,
           reverse_temp = temp.y)

  print("Filtering matches")
  filtered_combined <- combined_df %>%
    filter(abs(forward_temp - reverse_temp) <= 1) %>%
    mutate(heterodimer = calculate_dimer(forward_primer, reverse_primer)$dg * 0.004184)
    #filter(heterodimer > -7)

  df <- filtered_combined
  # if you want to take top 50
  df <- df %>%
    group_by(snpID, snp_character) %>%
    arrange(desc(heterodimer)) %>%
    slice_head(n = 10)
  df
  return(df)
}

get_complex_filter <- function(df){

  unique_groups <- df %>% distinct(snpID, snp_character)

  #match <- vector("logical", nrow(df))

  print("Filtering for temps between all pairs")

  # add column for tracking matches
  df$match <- ""

  # Iterate over each unique group
  for (first_index in 1:nrow(df)) {

    print(paste0(first_index, " out of ", nrow(df), " rows"))

    for (second_index in 1:nrow(df)) {
      if (df$snpID[first_index] == df$snpID[second_index] & df$snp_character[first_index] == df$snp_character[second_index]) {
        next
        }
      else {
        values <- c(df$forward_temp[first_index], df$reverse_temp[first_index],
                    df$forward_temp[second_index], df$reverse_temp[second_index])
        differences <- abs(diff(values))
        if (max(differences < 1)) {
          df$match[second_index] <- TRUE
        }
        }
      }
  }

  df <- filter(df, match == TRUE)
  df$match <- ""

  for (first_index in 1:nrow(df)) {

    print(paste0(first_index, " out of ", nrow(df), " rows"))

    for (second_index in 1:nrow(df)) {
      if (df$snpID[first_index] == df$snpID[second_index] & df$snp_character[first_index] == df$snp_character[second_index]) {
        next
      }
      else {

        hets <- c((calculate_dimer(df$forward_primer[first_index], df$forward_primer[second_index])$dg * 0.004184), (calculate_dimer(df$forward_primer[first_index], df$reverse_primer[second_index])$dg * 0.004184), (calculate_dimer(df$reverse_primer[first_index], df$forward_primer[second_index])$dg * 0.004184), (calculate_dimer(df$reverse_primer[first_index], df$reverse_primer[second_index])$dg * 0.004184))

        if (min(hets > -10)) {
          df$match[second_index] <- TRUE
        }
      }
    }
  }

  df <- filter(df, match == TRUE)
  df$groupID <- paste(df$snpID, df$snp_character, sep="-")
  df
  return(df)
     #Store the results in the list
    #comparison_results[[paste(current_group$snpID, current_group$snp_character, sep = "_")]] <- mean_diff
}

df <- get_primer_candidates(primer, shift)
df <- get_self_filter(df)
df <- get_cross_filter(df)
df <- get_complex_filter(df)
split_dfs <- split(df, df$groupID)

#generate_final_primers <- function(df){

#  df$groupID <- paste(df$snpID, df$snp_character, sep="-")

  #grouped_data <- split(df, df$groupID)

#  unique_groups <- unique(df$groupID)

  # Function to filter rows for a specific group
#  filter_group <- function(group) {
#    df %>% filter(groupID == group)
#  }

  # Use crossing to generate all combinations
#  combinations <- crossing(
#    lapply(unique_groups, filter_group)
#  )

  # Flatten the list of combinations into a single data frame
#  combined_results <- combinations %>%
#    rowwise() %>%
#    mutate(
#      combined = list(c_across(everything()))
#    ) %>%
#    ungroup() %>%
#    unnest_wider(combined)



#}
