## Input: MACADAM query results directory
## that contains a set of files "nametaxa_compact.tsv"
## each file represents the result of one query and its nametaxa is the name of the considered OTU
## Files are obtained as described in https://github.com/maloleboulch/MACADAM-Database
##
## and Host/OTUs abundances tables

# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

#LOAD LIBRARIES
#-----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(htm2txt)
library(klaR)
library(caret)


#LOAD DATA
#-------------------------------------------------------------------------------
### 

##Files with MACADAM QUERY RESULTS
files_to_read <- list.files(path = "your_path_to_MACADAM_query_results", pattern = "\\compact.tsv$", full.names = T)

all_path <- data.frame(Metabolic.Pathway=character(), Target.taxonomies=character(), MetaCyc.URL=character())

for(i in files_to_read){
  print(i)
  print(length(count.fields(i)))
  if(length(count.fields(i))>1){
  filepath <- read.csv(i,  header = TRUE, sep="\t",)
  ##filepath <- filepath[filepath$Median.score >0.5,] #eventually put a cutoff 
  filepath4 <- filepath %>% dplyr::select('Metabolic.Pathway', 'Target.taxonomies', 'MetaCyc.URL')
  all_path <- rbind(all_path,filepath4)
  }
}

all_path <- distinct(all_path)
write.csv2(all_path, "all_path.csv", row.names=FALSE)

buff <- all_path
buff$Metabolic.Pathway <- htm2txt(buff$Metabolic.Pathway)
all_path <- buff
write.csv2(all_path, "all_path_corrHTML.csv", row.names=FALSE)

###################
###
###################
##Files with HOST_OTUs_tables
files_to_read2 <- list.files(path = "your_path_to_HOST_OTUs_tables", pattern = "\\.csv$", full.names = T)
##File with list of all/interesting OTUs expressed in the host (may not coincide with MACADAM reuslts, some genus/species are not in MACADAM or not pathway annotated)
##Ideally, it is a two column file, one with short name (i.e. either falimy,genus,species), one with complete taxonomy
ListOtus <- read.csv2("your_path", header = F)

pathwayfile <- all_path

####################
colnames(ListOtu) <- c("short","complete")


##### Taxa missing taxa at next taxonomic level.
## This parsing part depends on how your name are written
missingOTUs <- as.character(ListOtu[!ListOtu$short %in% pathwayfile$Target.taxonomies,]$complete)
subst <- lapply(missingOTUs, function(x) tail(unlist(strsplit(x, " ")),2))
subst2 <- data.frame(matrix(subst, nrow=length(subst), byrow =TRUE))
names(subst2) <- c("short")
subst2$complete <- missingOTUs
subst2 <- subst2[!grepl("_",subst2$short),]
subst2$b1 <- matrix(unlist(subst2$short), nrow=length(subst2$short),  byrow =TRUE)[,1]
subst2$b2 <- matrix(unlist(subst2$short), nrow=length(subst2$short),  byrow =TRUE)[,2]


ListOtu$shortCorr <- lapply(ListOtu$short, function(x) as.character(subst2$b1[subst2$b2 == x][1]))
ListOtu$shortCorr2 <- ifelse(is.na(ListOtu$shortCorr),as.character(ListOtu$short),ListOtu$shortCorr)

######## End get missing taxa at next level

pathwayfile <-  pathwayfile[pathwayfile$Target.taxonomies %in% ListOtu$shortCorr2,] ##just to save space



ListOtu$Compl2 <- gsub("[[:space:]][[:space:]]", " ", ListOtu$complete)
ListOtu$Compl2 <- gsub("[[:space:]]", "\\.", ListOtu$Compl2)
ListOtu$Compl2 <- gsub("-", "\\.", ListOtu$Compl2)

###transform the file into wide and save as matrix
pathwayfile$val <- rep(1,length(pathwayfile$Metabolic.Pathway))
pathwayfileWide <- pathwayfile %>%  dplyr::select(-Metabolic.Pathway) %>% 
  pivot_wider(names_from=Target.taxonomies, values_from=val, values_fill = 0)

pathwayfileWideMatrix <- as.matrix(pathwayfileWide %>% dplyr::select(-MetaCyc.URL))
rownames(pathwayfileWideMatrix) <- pathwayfileWide$MetaCyc.URL
write.table(pathwayfileWideMatrix, "Matrix_Pathway_OTU.csv")

## for each OTUs table
for(namefile in files_to_read2){
  print(namefile)
  ## these are files with as first coloumn hostID, 
  fileLogabund <- read.csv(namefile,  header = TRUE, sep=",",)
  
  colnames(fileLogabund) <- gsub("\\.\\.", "\\.", colnames(fileLogabund))
  #colnames(fileLogabund) <- gsub("\\.[a-z]__", "\\1""__", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("\\.p__", "p__", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("\\.c__", "c__", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("\\.o__", "o__", colnames(fileLogabund))
  colnames(fileLogabund) <- gsub("__\\.", "__", colnames(fileLogabund))
  rownames(fileLogabund) <- lapply(seq(1:length(fileLogabund$hostID)), function(x) paste(fileLogabund$hostID[[x]],x, sep="_"))
  fileLogabund <- fileLogabund %>% dplyr::select(-hostID)
  
  ### Normalize if needed
  buff <- fileLogabund
  for(i in (1:dim(buff)[[1]])){
    vals <- as.data.frame(lapply(10^(buff[i,1:dim(buff)[[2]]])-10^-7, function(x) max(0,x))) 
    buff[i,1:dim(buff)[[2]]] <- vals/sum(vals)
    total2 <- sum(buff[i,1:dim(buff)[[2]]] )
  }
  
  fileabund <- buff
  fileabund$patients <- rownames(fileabund)
  fileabundLong <- pivot_longer(fileabund, !patients)
  fileabundLong$name2 <- lapply(fileabundLong$name, function(x) as.character(ListOtu$shortCorr2[ListOtu$Compl2 %in% x]))
  fileabundLong$name2 <- as.character(fileabundLong$name2)
  
  fileabundLong2 <- fileabundLong %>% dplyr::select(-name) %>% group_by(patients, name2) %>% summarise(value = sum(value))
  dim(distinct(fileabundLong %>% dplyr::select(-c(value))))
  dim(distinct(fileabundLong %>% dplyr::select(-c(name,value))))
  dim(distinct(fileabundLong2 %>% dplyr::select(-c(value))))
  ##this is to check I am not missing any value
  
  ##the result has to be numeric and not character. Otherwise matrix product gives error
  fileabund2 <- pivot_wider(fileabundLong2, names_from=name2, values_fill = 0) ##values_fill = 0 non dovrebbe servire
  print(colnames(fileabund2)[! colnames(fileabund2) %in% ListOtu$shortCorr2]) 
  buff2 <- as.matrix(fileabund2[,colnames(fileabund2) %in% ListOtu$shortCorr2])
  rownames(buff2) <- fileabund2$patients
  fileabund2 <- buff2

  ##
  XXX <- intersect(colnames(pathwayfileWideMatrix), colnames(fileabund2))
  pathwayfileWideMatrix2 <- pathwayfileWideMatrix[, XXX]
  pathwayfileWideMatrix2 <- pathwayfileWideMatrix2[rowSums(pathwayfileWideMatrix2)>0,]
  bindme <- rbind(pathwayfileWideMatrix2,   fileabund2[, XXX])
  
  
  product <- bindme %*% t(bindme)
  prod2 <- as.data.frame(product[1:dim(pathwayfileWideMatrix2)[[1]],-(1:dim(pathwayfileWideMatrix2)[[1]])])
  prod2 <- prod2[rowSums(prod2)>0,]
  write.table(prod2, paste("Pathway_Host_",tail(unlist(strsplit(namefile, "/")),1), sep=""),dec = ".", sep = ";", quote = FALSE)
} 


#################
## Algorithm and code 
## by Dr. Fumagalli Maria Rita and Stella Maria Saro
## 28th Feb. 2022 
## contact  mariarita.fumagalli@gmail.com
#################
