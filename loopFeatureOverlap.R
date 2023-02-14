#New version of loopFeatureOverlap.py, written in R using Granges instead of Python-based Pyranges version, which had an output issue

#Imports
library('plyr')
library('dplyr')
require('reshape2')
library('purrr')
library('grid')
# library('ChIPpeakAnno')
library('IRanges')
library('GenomicRanges')
library('arrangements')
library('foreach')

#Get args
if(!require(optparse)) {
  stop("Please install the optparse package and try again!", call. = FALSE)
}

library(optparse)

parser <- OptionParser(add_help_option = TRUE)
parser <- add_option(parser, c("-o", "--outdir"), action = "store", type = "character", help = "Directory to output bed files", default = "na")
parser <- add_option(parser, c("-l", "--loops"), action = "store", type = "character", help = "Input bedpe or tsv file containing loops in bedpe format", default = "na")
parser <- add_option(parser, c("-b", "--bed"), action = "store", type = "character", help = "Input bed file or files - if multiple, separate with commas", default = "na")
parser <- add_option(parser, c("-i", "--id"), action = "store", type = "character", help = "Feature names to use for each bed file - will be used to name output files. If multiple, separate with commas", default = "na")
parser <- add_option(parser, c("-e", "--exclusive"), action = "store_true", help = "When set, defines loop anchors with only a single feature - overlap with multiple features is not allowed", default = FALSE)
args <- parse_args(parser)

#########################################################################
#Convert args to variables
path.loops <- args$loops
outdir <- args$outdir
features <- args$bed
ids <- args$id

#########################################################################
#Read in loops
loops <- read.delim(path.loops, header = FALSE, col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2"), sep = "", dec = ".") #opens up the file

#Read in features
featureslist <- unlist(strsplit(features, ","))
idlist <- unlist(strsplit(ids, ","))

#Check features and ids are same length
if (length(featureslist) != length(idlist)) {
  stop("Make sure the same numbers of bed files and IDs are provided")
}

#We only want and know the names of the first three columns, so write a little import function:
import_feature_data <- function(filename) {
  #Get number of columns
  num.cols.to.blank <- max(count.fields(filename, sep = "\t")) - 3
  df <- read.delim(filename, header = FALSE, sep = "", dec = ".", colClasses = c('character', rep('numeric', 2), rep("NULL", num.cols.to.blank)))
  colnames(df) <- c("chr", "start", "end")
  return(df)
}

features.data <- lapply(featureslist, import_feature_data)

#Ensure output directory ends with a /
if(!endsWith(outdir, '/')) {
  outdir <- paste0(outdir, '/')
}

if (outdir == "na") {
  stop("Please provide an option for --out", call. = FALSE)
}

#########################################################################
#Make Granges objects
#For features
features.data.granges <- lapply(features.data, makeGRangesFromDataFrame, keep.extra.columns = FALSE)
#For loops
#First make separate dfs for each anchor
#Make anchors function
make_anchors_separate <- function(loops) {
  #Add loop_id column for merging
  loops$loop_id <- seq.int(nrow(loops))
  #Split loops into two bed-like files
  loops.1 <- data.frame(chr = loops$chr1, start = loops$start1, end = loops$end1, loop_id = loops$loop_id)
  loops.2 <- data.frame(chr = loops$chr2, start = loops$start2, end = loops$end2, loop_id = loops$loop_id)
  #Merge with original to add back the lost info
  return(list(loops.1, loops.2))
}
#Then make Granges
loops.anchors.list <- make_anchors_separate(loops)
#Make separate dfs (maybe don't need to do this?)
loops.anchors.1 <- loops.anchors.list[[1]]
loops.anchors.2 <- loops.anchors.list[[2]]
#Make granges
loops.anchors.list.granges <- lapply(loops.anchors.list, makeGRangesFromDataFrame, keep.extra.columns = FALSE)

#########################################################################
#Compare Granges objects
#Start counter
feature.count <- 1

#Loop through the features and count overlaps
for (feature in features.data.granges) {
    loops.anchors.1[[idlist[feature.count]]] <- countOverlaps(loops.anchors.list.granges[[1]], feature)
    loops.anchors.2[[idlist[feature.count]]] <- countOverlaps(loops.anchors.list.granges[[2]], feature)
    feature.count <- feature.count + 1
}


#Change column names for loops.anchors.2 so they don't match
colnames(loops.anchors.2) <- paste0(colnames(loops.anchors.2), '2')

#Next, merge based on loop_id(2) columns
loops.anchors.remerge <- merge(loops.anchors.1, loops.anchors.2, by.x = "loop_id", by.y = "loop_id2")

########################################################################
#Determine loop classes

#First need to generate loop classes based on ids
#Add null to the class list (for anchors with no features)
idlist.null <- append(idlist, "null")

#Then make all combinations (this is combinations with replacement)
position <- 1
combination.object <- icombinations(idlist.null, k = 2, replace = TRUE)

id.combined.list <- vector("list", length(combination.object$collect())/2)

foreach(x = icombinations(idlist.null, k = 2, replace = TRUE), .combine = c) %do% {
  id.combined.list[[position]] <- paste(idlist.null[x[1]], idlist.null[x[2]], sep = "-")
  position <- position + 1
}

#Classify loops. Which version is run depends on whether inclusive or exclusive loops are desired.
if (args$exclusive) {
  #Make named list of dfs
  output.list <- setNames(replicate(length(id.combined.list), data.frame()), id.combined.list)

  #Exclusive version
  for (i1 in 1:length(idlist)) {
    item1 <- idlist[[i1]]
    #print(item1)
    for (i2 in 1:length(idlist)) {
      item2 <- idlist[[i2]]
      #print(item2)
      if (i1 == i2 | i1 < i2) {
        #Get the loop type from the indices
        looptype <- paste(idlist[[i1]], idlist[[i2]], sep = '-')
        #Get the relevant loops
        temp.df <- loops.anchors.remerge
        #Loop through all ids, check they're 0 except the ones matching the requirements
        for (id in idlist) {
          if (id == idlist[[i1]]) {
            temp.df <- temp.df[which(temp.df[[id]] > 0),]
          } else if (id != idlist[[i1]]) {
            temp.df <- temp.df[which(temp.df[[id]] == 0),]
          }
          if (id == idlist[[i2]]) {
            temp.df <- temp.df[which(temp.df[[paste0(id, '2')]] > 0),]
          } else if (id != idlist[[i2]]) {
            temp.df <- temp.df[which(temp.df[[paste0(id, '2')]] == 0),]
          }
        }
        output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)

      } else if (i1 > i2) {
        #Get the loop type from the indices - need to invert here so that E-P2 and P-E2 loops are both put into
        #the same category
        looptype <- paste(idlist[[i2]], idlist[[i1]], sep = '-')
        #Get the relevant loops
        temp.df <- loops.anchors.remerge

        for (id in idlist) {
          if (id == idlist[[i1]]) {
            temp.df <- temp.df[which(temp.df[[id]] > 0),]
          } else if (id != idlist[[i1]]) {
            temp.df <- temp.df[which(temp.df[[id]] == 0),]
          }
          if (id == idlist[[i2]]) {
            temp.df <- temp.df[which(temp.df[[paste0(id, '2')]] > 0),]
          } else if (id != idlist[[i2]]) {
            temp.df <- temp.df[which(temp.df[[paste0(id, '2')]] == 0),]
          }
        }
        output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)
      }
    }
    looptype <- paste(idlist[[i1]], "null", sep = '-')
    #Generate X-null
    temp.df.1 <- loops.anchors.remerge
    #Loop through the ID columns of the other anchor and select only rows with 0 for each column
    for (id in idlist) {
      if (id == idlist[[i1]]) {
        temp.df.1 <- temp.df.1[which(temp.df.1[[id]] > 0),]
      } else if (id != idlist[[i1]]) {
        temp.df.1 <- temp.df.1[which(temp.df.1[[id]] == 0),]
      }
      temp.df.1 <- temp.df.1[which(temp.df.1[[paste0(id, '2')]] == 0),]
    }
    #Generate null-X
    temp.df.2 <- loops.anchors.remerge[which(loops.anchors.remerge[[paste0(idlist[[i1]], '2')]] > 0),]

    for (id in idlist) {
      temp.df.2 <- temp.df.2[which(temp.df.2[[id]] == 0),]
      if (id == idlist[[i1]]) {
        temp.df.2 <- temp.df.2[which(temp.df.2[[paste0(id, '2')]] > 0),]
      } else if (id != idlist[[i1]]) {
        temp.df.2 <- temp.df.2[which(temp.df.2[[paste0(id, '2')]] == 0),]
      }
    }
    #Combine them
    temp.df <- rbind(temp.df.1, temp.df.2)
    output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)

  }
  looptype <- "null-null"
  temp.df <- loops.anchors.remerge

  for (id in idlist) {
    temp.df <- temp.df[which(temp.df[[id]] == 0 & temp.df[[paste0(id, '2')]] == 0),]
  }
  output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)
} else {
  #Make named list of dfs
  output.list <- setNames(replicate(length(id.combined.list), data.frame()), id.combined.list)

  #Inclusive version
  for (i1 in 1:length(idlist)) {
    item1 <- idlist[[i1]]
    #print(item1)
    for (i2 in 1:length(idlist)) {
      item2 <- idlist[[i2]]
      #print(item2)
      if (i1 == i2 | i1 < i2) {
        #Get the loop type from the indices
        looptype <- paste(idlist[[i1]], idlist[[i2]], sep = '-')
        #Get the relevant loops
        temp.df <- loops.anchors.remerge[which(loops.anchors.remerge[[idlist[[i1]]]] > 0 & loops.anchors.remerge[[paste0(idlist[[i2]], "2")]] > 0), ]
        output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)
      } else if (i1 > i2) {
        #Get the loop type from the indices - need to invert here so that E-P2 and P-E2 loops are both put into
        #the same category
        looptype <- paste(idlist[[i2]], idlist[[i1]], sep = '-')
        #Get the relevant loops
        temp.df <- loops.anchors.remerge[which(loops.anchors.remerge[[idlist[[i1]]]] > 0 & loops.anchors.remerge[[paste0(idlist[[i2]], "2")]] > 0), ]
        output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)
      }
    }
    looptype <- paste(idlist[[i1]], "null", sep = '-')
    #Generate X-null
    temp.df.1 <- loops.anchors.remerge[which(loops.anchors.remerge[[idlist[[i1]]]] > 0),]
    #Loop through the ID columns of the other anchor and select only rows with 0 for each column
    for (id in idlist) {
      temp.df.1 <- temp.df.1[which(temp.df.1[[paste0(id, '2')]] == 0),]
    }
    #Generate null-X
    temp.df.2 <- loops.anchors.remerge[which(loops.anchors.remerge[[paste0(idlist[[i1]], '2')]] > 0),]

    for (id in idlist) {
      temp.df.2 <- temp.df.2[which(temp.df.2[[id]] == 0),]
    }
    #Combine them
    temp.df <- rbind(temp.df.1, temp.df.2)
    output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)

  }
  looptype <- "null-null"
  temp.df <- loops.anchors.remerge

  for (id in idlist) {
    temp.df <- temp.df[which(temp.df[[id]] == 0 & temp.df[[paste0(id, '2')]] == 0),]
  }
  output.list[[looptype]] <- rbind(output.list[[looptype]], temp.df)
}

#Remove any duplicates from the dfs (inclusive calling can result in many)
output.list.nodups <- lapply(output.list, distinct)

#Print the lengs of each df to give the result:
sapply(output.list.nodups, nrow)

#Output the dataframes
for (i in 1:length(output.list.nodups)) {
  temp.df <- output.list.nodups[[i]]
  #Get only the relevant columns
  output.df <- data.frame(chr = temp.df$chr, start = temp.df$start, end = temp.df$end, chr2 = temp.df$chr2, start2 = temp.df$start2, end2 = temp.df$end2)
  write.table(output.df, file = paste0(outdir, names(output.list)[[i]], '.bedpe'), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}
