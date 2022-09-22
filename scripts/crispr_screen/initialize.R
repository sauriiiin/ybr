##### INITIALIZE PACKAGES AND FUNCTIONS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 03/01/2022 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(ggridges)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(plotly)
library(scales)
library(reshape2)
library(rstatix)
library(gtools)
library(effsize)
library(png)
library(ggstatsplot)
library(cowplot)
library(pzfx)
library(stringi)
library(dplyr)
library(RMariaDB)
library(readxl)
source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")

conn <- initialize.sql("saurin_test")

fig_path <- "~/R/Projects/ybr/output/figs/"

##### FIGURE SIZE (CELL PRESS)
one.c <- 85 #single column
one.5c <- 114 #1.5 column
two.c <- 174 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 8

##### LEGEND GRABBER
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

`%notin%` <- Negate(`%in%`)

