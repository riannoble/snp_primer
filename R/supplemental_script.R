library(biomaRt)
library(tidyverse)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")

# Specify dataset and attributes for human SNPs
dataset <- "hsapiens_snp"
attributes <- c("refsnp_id", "chr_name", "chrom_start", "consequence_type", "allele_string")

# Retrieve SNP data for a specific accession number (e.g., rs1815739)
accession_number <- "rs1815739"

snp_data <- getBM(attributes = attributes,
                  filters = "refsnp_id",
                  values = accession_number,
                  mart = ensembl)

# Print the retrieved SNP data
print(snp_data)

# straightforward
entrez <-
getSequence(id = entrez,
            type="entrezgene",
            seqType="coding_gene_flank",
            upstream=100,
            mart=ensembl)
