
source("C:/Users/riano/Documents/acornfinder/acornfinder/R/fxns_draft.R")

# TESTING SECTION

primer = "rs9462492, rs58318008"
shift = 100

df <- get_primer_candidates(primer, shift)
df <- get_self_filter(df)
new_df <- get_cross_filter(df)
final_df <- get_final_list(new_df)

# downloading the dataframe
final_df <- final_df[, !sapply(final_df, is.list)]
write.csv(final_df)

