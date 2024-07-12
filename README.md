# acornfinder

This is an extension of Archarlie's Shiny App that creates SNP-sensitive primers and probes (no probes as of July 2024). This versions attempts to condense it into a package.

## FOR FUTURE PROGRAMMERS!

-talk to Bro. Terribilini throughout your 

### Installation

**Copy and paste the code below:** (Work in Progress)

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
# Obesity Primers
primer = "rs9462492, rs58318008, rs1421085"
primer = "rs1550576, rs17025867"
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

## Non-Essential Functions (not complete)

### mart_api function

Parameters: primer, shift

Result: generates primers

Example: 

Future:
- create a help file?
install.packages("roxygen2")
library(roxygen2)
