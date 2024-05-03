# acornfinder

This is an extension of Archarlie's Shiny App that creates SNP-sensitive primers and probes. This versions attempts to condense it into a package.

### Installation

**Copy and paste the code below:** (Work in Progress)

install.packages("devtools")
library(devtools)
install_github("rianoble/acornfinder")
library(acornfinder)

snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")


## Supplemental Installation

If you cannot use the package correctly from just the code above, paste the code below (contains packages) as well. You may need troubleshooting and installation help from a data scientist or bio/statistician (ex. Ria, Brother Terribilini).

require("rprimer")
require("htmltools")
require("kableExtra")
require("DT")
require("dplyr")
require("tidyverse")
require("stringi")
require("stringr")
require("mosaic")
require("purrr")
require("ggplot2")
require("hexbin")
require("patchwork")
require("plotly")
require("BiocManager")
require("biomaRt")
require("spgs")
require("primer3")


### findacorn function

Parameters: primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin

Result: primers

Example: 

primer = "rs53576, rs1815739, rs7412, rs429358, rs6152"
shift = 100
desired_tm = 64
diff = 3
Heterodimer_tm = 50
Homodimer <- 45
top <- 2
hairpin <- 45
output <- findacorn(primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin)

Attempts to answer all questions in one function. Produces a pop-out output.

-reverse primer has same melting temperature as forward primer
-melting temps should be under 1-2 degrees difference
-paste in entire sequence retrieved by biomart in r into primer blast (~1000 bases of seq, 500 before, 500 after)

## Non-Essential Functions

### mart_api function

Parameters: primer, shift

Result: generates primers

Example: 

Future:
- create a help file?
install.packages("roxygen2")
library(roxygen2)
