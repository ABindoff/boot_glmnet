library(openxlsx)
library(ggplot2)
library(tidyverse)

d <- read.xlsx('Proteomics raw data.xlsx', colNames = FALSE)
# label Batch, Time, and Sex rows
d[1, 1] <- "Batch"
d[2, 1] <- "Time"
d[3, 1] <- "Sex"

# transpose the entire table
d <- t(d)

# save the top row for column names
cn <- d[1,]
cn[is.na(cn)] <- paste0("FakeID", 1:10)
# look for duplicates
cn[duplicated(cn)]

# delete top row
d <- d[-1,]

# convert to data frame
d <- data.frame(d) 

# name all columns
names(d) <- cn

# remove duplicates
d <- data.frame(d) %>% dplyr::select(-cn[duplicated(cn)])

# make sure everything is a floating point number
d[,5:2373] <- lapply(d[,5:2373], as.numeric)

# scale to unit SD
d[,5:2373] <- lapply(d[,5:2373], scale)

# convert Batch and Sex to 0,1
d$Batch <- ifelse(d$Batch == "Batch 1", 0, 1)
d$Sex <- ifelse(d$Sex == "M", 0, 1)
d$tbi <- ifelse(d$Time == "N", 0, 1)
d$Time <- factor(d$Time, levels = c("N", "28d", "7d", "3d"))

