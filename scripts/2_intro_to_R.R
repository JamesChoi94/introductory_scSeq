
####### Interrogating scSeq data 2020/07/07 #######

# The "#" symbol is used to denote the start of a "comment". 
# Comments are text that will not be run as code i.e. they are 
# a handy way to write notes.



# Functions and variables -----------------------------------------------

# Everything you will work with in R is either a variable or a
# function. Variables are data that you can reference with a name.
# Functions are operations (e.g. addition) that can be performed
# on variables.


# "<-" arrow symbol is the assignment operator. It assigns a value
# to a name.


# To run a specific line of code: A) Click cursor on line you want to
# run and press Ctrl + Enter, or B) Copy + paste the line into Console.


# Examples of assigning values to variables

a <- 1
b <- 2.0
c <- 'text'
d <- FALSE


# There are different types of "data structures" that can hold
# multiple data values in an organized way.

x <- c(1,2,3) # numeric vector
y <- c('A','B','C') # character vector
z <- data.frame(x, y) # data.frame


# Use single or double brackets to extract data values from 
# data structures.

x[1] # extract value at position 1
y[c(1,3)] # extract values at positions 1, 2, and 3
z[[1]] # extract column 1


# Functions can be identified by the use of parentheses.

a <- sum(1, 2)
b <- c(1,2,3) # c() itself is a function



# Loading packages --------------------------------------------------------

# Packages are collections of data/files that are pre-packaged into one
# download-able entity. Packages usually contain functions that are tailored
# towards a specific task/analysis. "Seurat" is a package developed/maintained
# by the Satija lab at NYU that is designed for single-cell analyses.


# Load packages into your working session with the library() function.
library('Seurat')
library('dplyr')
library('ggplot2')



# Directories -------------------------------------------------------------

# Working directories
getwd()

# Set your working directory. All saved files/data will automatically be
# saved in this location.
setwd()


# Loading data
macroglia <- readRDS(file = 'data/macroglia.rds')



