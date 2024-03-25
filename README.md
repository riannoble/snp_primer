# acornfinder

This is an extension of Archarlie's Shiny App that creates SNP-sensitive primers and probes. This versions attempts to condense it into a package.

### Installation

**Copy the code below:**

install.packages("devtools")
library(devtools)
install_github("rianoble/acornfinder")
library(acornfinder)
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")



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
findacorn(primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin)

Attempts to answer all questions in one function. Produces a pop-out output.

### mart_api function

Parameters: primer, shift

Result: generates primers

Example: 

