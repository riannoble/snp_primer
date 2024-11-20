# acornfinder

This is an extension of Archarlie's Shiny App that creates SNP-sensitive primers and probes (no probes as of July 2024).

## FOR FUTURE PROGRAMMERS!

- Main files are ria_fxns_draft and the app
- Talk to Bro. Terribilini throughout your development of the code and app
- Next steps are optimizing the app/code and creating a package
- My multiplex testing primers were held on to by Bro. Terribilini so you can test the specificity of the snp-specific primer design from the code
- You will need to find a better temp calculation method (check accuracy of your method with Oligo Analyzer)
- Better results could be found if you utilize parallelization
- Good luck on your work! If you have any questions for me, my email is rianoble@umich.edu :)

### Installation of Potential Package

**Copy and paste the code below:** (package is a goal for the future)

install.packages("devtools")
library(devtools)
install_github("rianoble/acornfinder")
library(acornfinder)

### findacorn function

Parameters: primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin

Result: table of primers

Example: 

# Obesity Primer Options
primer = "rs53576, rs1815739, rs7412, rs429358, rs6152"
primer = "rs9462492, rs58318008, rs1421085"
primer = "rs1550576, rs17025867"
primer = "rs9462492, rs58318008"
primer = "rs1421085, rs9462492"

shift = 100

The parameters above should be utilized, but not the ones below.

 desired_tm = 64
 diff = 3
 Heterodimer_tm = 50
 Homodimer <- 45
 top <- 2
 hairpin <- 45

The parameters above existed before I took on the project, I am leaving them here just in case someone decides to bring back previous elements of the file. However, they are not needed currently.

output <- findacorn(primer, shift, desired_tm, diff, Heterodimer_tm, Homodimer, top, hairpin)

^ had an idea of compiling everything into one function in the package
Clean run can be utilized through running.R

**Rules**
-reverse primer has same melting temperature as forward primer
-melting temps should be under 1-2 degrees difference
-paste in entire sequence retrieved by biomart in r into primer blast (~1000 bases of seq, 500 before, 500 after)

## Non-Essential Functions (document any (component) functions below)

