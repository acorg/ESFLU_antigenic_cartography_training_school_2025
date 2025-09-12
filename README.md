# ESFLU_antigenic_cartography_training_school_2025
This repository contains the materials for the 2025 ESFLU training school on antigenic cartography.

## Before the workshop

We will be working with free software.

The workshop is run in R. Please make sure you have a working installation of R 
on your computer. Instructions for download and installation can be found 
on the official R website: https://cran.r-project.org

Make sure you install the correct R for your operating system.

To install the `ablandscapes` package, you need to have development tools installed for your operating system. Instructions 
are linked on the R page. Here are the instructions to install the tools for MacOS, XCode and gfortran are needed: https://mac.r-project.org/tools/ 

If installation issues occur, it is often because the R version and the development tools version do not match. We recommend 
installing the most recent version of both and updating R if necessary.

To screenshot antibody landscapes, we will use a package that requires installation of [Google Chrome](https://www.google.com/chrome/) 
on your device.

We also recommend installing RStudio as it provides an accessible user interface: https://posit.co/products/open-source/rstudio/

Once you have everything installed, please open a session in RStudio and run the script `code/00_installation.R` to install all 
packages that we will be using throughout the workshop.

For a python interface to Racmacs, please follow the instructions here: https://github.com/iAvicenna/PyRacmacs 

To download this GitHub repository, click on the green `<> Code` icon on the top right and select `Download ZIP`
in the dropdown menu. You can then create a new `.Rproject` in Rstudio with the unzipped file. To do so, 
open RStudio and select `File > Open Project ... ` and open the unzipped file.

### Coding experience

Coding experience is recommended, but not required. Participants without any coding experience can use the 
`_solutions.Rmd` provided for each exercise. 

Here are links to an [R](https://www.geeksforgeeks.org/r-language/r-tutorial/) and [RStudio tutorial for beginners](https://www.datacamp.com/tutorial/r-studio-tutorial). 
Working through these is *not* a requirement for participating in the workshop.

## During the workshop

The code for each session is in the `code` directory and numbered. We will provide working code for each 
session in the same directory with the appendix `_solution.Rmd` but encourage 
participants to try each task by themselves first. 

