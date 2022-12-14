# Clear workspace
### Input: from algorithm part 1 and 2

# ------------------------------------------------------------------------------
rm(list = ls())

#LOAD LIBRARIES
#-----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(klaR)
library(caret)

#LOAD DATA
#-------------------------------------------------------------------------------
###

files_to_read <- list.files(path = "your_path", pattern = "^Pathway_Host\\w+\\.csv$", full.names = T)

## the matrix pathways/score 
MetaboliteScore <- read.csv2("your_pathwy/Pathway_reactions_score_table.csv")

##
M2 <-  MetaboliteScore %>% dplyr::select(-c(Pathway, reaction_name)) %>% group_by(ID,reagent)
MetaboliteScoreSum <- MetaboliteScore %>% dplyr::select(-c(Pathway, reaction_name)) %>% group_by(ID,reagent) %>% summarise(score = sum(score))

MetaboliteScoreSumWide <- pivot_wider(MetaboliteScoreSum, names_from=reagent, values_from = score, values_fill = 0, values_fn = list(score = as.numeric))
ggg <- ungroup(MetaboliteScoreSumWide) %>% dplyr::select(-ID)
MetaboliteScoreSumWide <- MetaboliteScoreSumWide[,c(TRUE,colSums(ggg) !=0)]  
rm(ggg)

for(namefile in files_to_read){
  print(namefile)
  filePathOTU <- read.csv(namefile,  header = TRUE, sep=";",)
  dim(filePathOTU)
  
  XXX <- intersect(MetaboliteScoreSumWide$ID, rownames(filePathOTU))
  MetaboliteScoreSumWide2 <- MetaboliteScoreSumWide[MetaboliteScoreSumWide$ID %in% XXX,]
  ggg2 <- ungroup(MetaboliteScoreSumWide2) %>% dplyr::select(-ID)
  MetaboliteScoreSumWide2 <- MetaboliteScoreSumWide2[,c(TRUE,colSums(ggg2) !=0)]
  rm(ggg2)
  MetaboliteScoreSumWideMatrix <- as.matrix(ungroup(MetaboliteScoreSumWide2) %>% dplyr::select(-ID))
  rownames(MetaboliteScoreSumWideMatrix) <- MetaboliteScoreSumWide2$ID
  filePathOTUMatrix <- as.matrix(filePathOTU[XXX,])
  
  bindme <- cbind(MetaboliteScoreSumWideMatrix,   filePathOTUMatrix)
  product <- t(bindme) %*% bindme
  prod2 <- as.data.frame(product[1:dim(MetaboliteScoreSumWideMatrix)[[2]],-(1:dim(MetaboliteScoreSumWideMatrix)[[2]])])
  
  write.table(prod2, paste("Metabolite_Host_",tail(unlist(strsplit(namefile, "Pathway_Host_")),1), sep=""),dec = ".", sep = ";", quote = FALSE)
  
  
  ###Preparing tables for LDA
  
  namecut <- gsub("\\.csv", "_Metabolite", tail(unlist(strsplit(namefile, "Pathway_Host_")),1))
  for_lda <- as.data.frame(product[-(1:dim(MetaboliteScoreSumWideMatrix)[[2]]),1:dim(MetaboliteScoreSumWideMatrix)[[2]]])
  for_lda_corr <- for_lda
  for_lda$disease_name <- gsub("_\\d+","", rownames(for_lda))
  write.table(for_lda, paste(namecut, ".csv",sep=""),dec = ".", sep = ";", quote = FALSE, row.names=FALSE)
  

  ### extract one element for each correlation cluster.
  ### procedure repeated N times in order to compare extactions
  
  for(resortreplica in 1:N){
    listS <- sample(1:ncol(for_lda_corr))
    for_lda_corr <- for_lda_corr[,listS]
    write.table(for_lda, paste("./your_output_path/",namecut, "_resorted_",resortreplica,".csv",sep=""),dec = ".", sep = ";", quote = FALSE, row.names=FALSE)

    corr_m1 <- cor(for_lda_corr)
    
    cccUP <- apply(corr_m1,2,function(x) which(abs(x) >= 0.95))
    names_to_keep <- list()
    allnames <- list()
    doubleval <- list()
    for(i in names(cccUP)){write.table(t(sort(rownames(as.data.frame(cccUP[i])))), 
                                       paste("./your_output_path/Correlated_",namecut, "_resorted_corr095_",resortreplica,".csv",sep=""),
                                       append=TRUE, sep = ";", col.names = FALSE, row.names = FALSE)
      if(! i %in% allnames){
        names_to_keep <- append(names_to_keep,i)
        allnames <- append(allnames,rownames(as.data.frame(cccUP[i])) )
        allnames <- unique(allnames)
      }else{ 
        doubleval <- append(doubleval,i)
        doubleval <- unique(doubleval)
      }
      
      
    }
    if(length(allnames) == dim(for_lda_corr)[[2]]){print("OK")}
    write.table(sort(unlist(names_to_keep)),paste("./your_output_path/NamesCorrelated_",namecut, "_resorted_corr095_",resortreplica,".csv",sep=""),
                append=TRUE, sep = ";", col.names = FALSE, row.names = FALSE)
    
    for_lda_corr2 <- for_lda_corr[,names(for_lda_corr) %in% names_to_keep ]
    the_ones_correlated <- for_lda_corr[,!names(for_lda_corr) %in% names(for_lda_corr2) ]
    
    for_lda_corr2$disease_name <- gsub("_\\d+","", rownames(for_lda_corr2))
    write.table(for_lda_corr2, paste("./your_output_path/",namecut, "_resorted_corr095_",resortreplica,".csv",sep=""),dec = ".", sep = ";", quote = FALSE, row.names=FALSE)
  }
  
  
}

#################
## Algorithm and code 
## by Dr. Fumagalli Maria Rita and Stella Maria Saro
## 28th Feb. 2022
## contact  mariarita.fumagalli@gmail.com
#################