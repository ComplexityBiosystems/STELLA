#------------------------------------------------------------------------------
#  Input: MACADAM query .tsv files
# Output: pathway-OTU matrix from MACADAM query results (.tsv files)
#         These files are obtained as described in
#         https://github.com/maloleboulch/MACADAM-Database and are saved in the
#         `tmp` directory, since it is used later on and it is not the final
#         result
#------------------------------------------------------------------------------
# Authors: MR Fumagalli mariarita.fumagalli@gmail.com
#          MS Saro      stellamaria.saro@studenti.unimi.it
#------------------------------------------------------------------------------

# Import Libraries
library(dplyr)   # data manipulation
library(tidyr)   # tidy data
library(tibble)  # dataframe
library(stringr) # work with strings
library(htm2txt) # text formatting
#
library(caret)

# More precise errors when using Rscript on the terminal
options(error = quote({
  dump.frames(to.file=TRUE, dumpto="last.dump")
  load("last.dump.rda")
  print(last.dump)
  q()
}))

# Clear workspace, MATLAB style!
rm(list = ls())


#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
# We extract pathways and OTUs from MACADAM query (.tsv) files
files_to_read <- list.files(path="../../data/MACADAM/Results", pattern="\\compact.tsv$", full.names=TRUE)

# create dataframe with 3 columns to store data
all_pathways <- data.frame(Metabolic.Pathway=character(), Target.taxonomies=character(), MetaCyc.URL=character())

for(f in files_to_read){
  if(length(count.fields(f))>1){
    pathway <- read.csv(f, header=TRUE, sep="\t",)
    # select only 3 rows from `pathway`
    pathwayFiltered <- pathway %>% dplyr::select("Metabolic.Pathway", "Target.taxonomies", "MetaCyc.URL")
    # add previous data to dataframe
    all_pathways <- rbind(all_pathways, pathwayFiltered)
  }
}

all_pathways <- distinct(all_pathways) # avoid redundancy
write.csv2(all_pathways, "tmp/all_pathways.csv", row.names=FALSE)

# NB all_pathways.csv has mispelled entries with HTML tags
# `htm2txt` fixes the mispellings (e.g., &beta changes to β)
buff <- all_pathways
buff$Metabolic.Pathway <- htm2txt(buff$Metabolic.Pathway) # proper renaming
all_pathways <- buff
write.csv2(all_pathways, "tmp/all_pathways_corrHTML.csv", row.names=FALSE)


# Files with host-OTUs tables
# Note that these tables may not coincide with MACADAM results, nor # have the
# pathway annotated). Ideally this is a 2-column file, one with short name
# (i.e. either family, genus or species) and one with complete taxonomy
# further columns are not used
path_ms <- "../../data/MACADAM/Filebase/MS/BATTERI_CORRELATI"
path_ms <- "../../data/MACADAM/Filebase/MS/BATTERI_NON_CORRELATI"
files_to_read2 <- list.files(path=path_ms, pattern="\\.csv$", full.names=TRUE)

# List of interesting OTUs
OTUs_ms = "../../data/MACADAM/ListOTUs/Species_Genus_Fam_complete_MS_name_corretto_a_mano.txt"
ListOtu <- read.csv2(OTUs_ms, header=FALSE)

colnames(ListOtu) <- c("short", "complete") # rename dataframe columns
pathwaysfile <- all_pathways

# Missing taxa at next taxonomic level
# This parsing part depends on how the names are written
missingOTUs <- as.character(ListOtu[!ListOtu$short %in% pathwaysfile$Target.taxonomies,]$complete)

# take the last 2 elements (words?) of missingOTUs...
subst <- lapply(missingOTUs, function(x) tail(unlist(strsplit(x, " ")), 2))

# ...and create a new dataframe with these elements only
subst2 <- data.frame(matrix(subst, nrow=length(subst), byrow=TRUE))
names(subst2) <- c("short")
subst2$complete <- missingOTUs
subst2 <- subst2[!grepl("_",subst2$short),] # select elements without "_" in the "short" column
subst2$b1 <- matrix(unlist(subst2$short), nrow=length(subst2$short), byrow=TRUE)[,1]
subst2$b2 <- matrix(unlist(subst2$short), nrow=length(subst2$short), byrow=TRUE)[,2]

ListOtu$shortCorr <- lapply(ListOtu$short, function(x) as.character(subst2$b1[subst2$b2==x][1]))
ListOtu$shortCorr2 <- ifelse(is.na(ListOtu$shortCorr), as.character(ListOtu$short), ListOtu$shortCorr)

# Save space
pathwaysfile <- pathwaysfile[pathwaysfile$Target.taxonomies %in% ListOtu$shortCorr2,]
# Note that now the number of rows of pathwaysfile is almost halved (asd)

# Adjust spelling (spaces, hyphens and dots)
ListOtu$Compl2 <- gsub("[[:space:]][[:space:]]", " ", ListOtu$complete)
ListOtu$Compl2 <- gsub("[[:space:]]", "\\.", ListOtu$Compl2)
ListOtu$Compl2 <- gsub("-", "\\.", ListOtu$Compl2)

# Transform the file into a wide one and save as dataframe with also OTUs
# columns matching with MetaCyc URLs
pathwaysfile$val <- rep(1,length(pathwaysfile$Metabolic.Pathway))
pathwaysfileWide <- pathwaysfile %>%
      	            dplyr::select(-Metabolic.Pathway) %>%
		                pivot_wider(names_from=Target.taxonomies, values_from=val, values_fill=0)

# transform dataframe into a matrix
pathwaysfileWideMatrix <- as.matrix(pathwaysfileWide %>%
			  dplyr::select(-MetaCyc.URL)) # drop column
rownames(pathwaysfileWideMatrix) <- pathwaysfileWide$MetaCyc.URL
write.table(pathwaysfileWideMatrix, "tmp/Matrix_Pathway_OTU.csv")

# for each OTUs table
for(namefile in files_to_read2){
  fileLogabund <- read.csv(namefile, header=TRUE, sep=",") # warning: the separator used is "," (comma)

  # adjust name entries
  colnames(fileLogabund) <- gsub("\\.\\.", "\\.", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("\\.p__", "p__", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("\\.c__", "c__", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("\\.o__", "o__", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("__\\.", "__", colnames(fileLogabund))

  rownames(fileLogabund) <- lapply(seq(1:length(fileLogabund$disease_name)), function(x) paste(fileLogabund$disease_name[[x]], x, sep="_"))
  fileLogabund <- fileLogabund %>%
	          dplyr::select(-disease_name) # drop column

  # scroll through all the patients
  buff <- fileLogabund
  for(i in (1:dim(buff)[[1]])){
    # save an array (vals) inside the last cell of the previous df ?!?!?!
    vals <- as.data.frame(lapply(10^(buff[i, 1:dim(buff)[[2]]])-10^-7, function(x) max(0, x)))
    buff[i, 1:dim(buff)[[2]]] <- vals/sum(vals) # assign to the last element of each row normalized `vals`
  }

  # widen database
  fileabund <- buff
  fileabund$patients <- rownames(fileabund)
  fileabundLong <- pivot_longer(fileabund, !patients)
  fileabundLong$name2 <- lapply(fileabundLong$name, function(x) as.character(ListOtu$shortCorr2[ListOtu$Compl2 %in% x]))
  fileabundLong$name2 <- as.character(fileabundLong$name2)

  # Check if there is no missing value
  fileabundLong2 <- fileabundLong %>%
                    dplyr::select(-name) %>% # drop column
                    group_by(patients, name2) %>%
                    summarise(value=sum(value))

  # The result has to be numeric and not char, otherwise matrix product gives error
  fileabund2 <- pivot_wider(fileabundLong2, names_from=name2, values_fill=0)
  buff2 <- as.matrix(fileabund2[,colnames(fileabund2) %in% ListOtu$shortCorr2])
  rownames(buff2) <- fileabund2$patients
  fileabund2 <- buff2

  # Is the pathway present in the taxonomy?
  XXX <- intersect(colnames(pathwaysfileWideMatrix), colnames(fileabund2))
  pathwaysfileWideMatrix2 <- pathwaysfileWideMatrix[, XXX]
  pathwaysfileWideMatrix2 <- pathwaysfileWideMatrix2[rowSums(pathwaysfileWideMatrix2)>0,]
  bindme <- rbind(pathwaysfileWideMatrix2, fileabund2[, XXX])

  product <- bindme %*% t(bindme)
  write.table(product, paste("product-", tail(unlist(strsplit(namefile, "/")),1), sep=""), dec=".", sep=";", quote=FALSE)
  prod2 <- as.data.frame(product[1:dim(pathwaysfileWideMatrix2)[[1]], -(1:dim(pathwaysfileWideMatrix2)[[1]])])
  prod2 <- prod2[rowSums(prod2)>0,]
  write.table(prod2, paste("tmp/Pathway_Host_", tail(unlist(strsplit(namefile, "/")),1), sep=""), dec=".", sep=";", quote=FALSE)
  write.table(prod2, paste("non-corr-", tail(unlist(strsplit(namefile, "/")),1), sep=""), dec=".", sep=";", quote=FALSE)
}
