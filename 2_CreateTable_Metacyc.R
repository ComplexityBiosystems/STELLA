## Input: a table with coloumnd c("Pathways","AllReactions","Spec", "ID") 
## that can be easily obtained from metacyc SmartTable selecting
## Pathways;Reaction-List;Species;Object ID
## Data should be downloaded "as common names"
##
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

#############################


# Load libraries
# ------------------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library("openxlsx")
library(htm2txt)


file1 <- read_tsv("./your_path/All-pathways-of-MetaCyc_as_common_names.txt")
colnames(file1) <- c("Pathways","AllReactions","Spec", "ID")   

file2 <- file1 %>% dplyr::select(-Pathways, -AllReactions )
file1 <- file1 %>% dplyr::select(-Spec)

buff <- file1
for(i in (1:dim(buff)[[1]])){
  print(i)
  buff[i,] <- lapply(buff[i,], function(x) htm2txt(x))
  print(buff[i,1])
}
file1 <- buff

 buff <- file2
 for(i in (1:dim(buff)[[1]])){
   print(i)
   buff[i,] <- lapply(buff[i,], function(x) htm2txt(x))
   print(buff[i,1])
 }
 file2 <- buff

react_to_path <- file1 %>%
  separate(
    AllReactions,
    into = as.character(seq(1:300)),
    sep = " // ",
    remove = FALSE,
    convert = FALSE,
    extra = "warn",
    fill = "right"
  )

react_to_path2 <- react_to_path %>% gather(Num,ReactionName, -c(Pathways,AllReactions,ID), na.rm = TRUE) %>% dplyr::select(-Num)

react_to_path2 <- distinct(react_to_path2)

write_csv2(react_to_path2, "All-pathways-of-MetaCyc_as_common_names_reactions_noHTML.txt")



species_to_path <- file2 %>%
  separate(
    Spec,
    into = as.character(seq(1:300)),
    sep = " // ",
    remove = FALSE,
    convert = FALSE,
    extra = "warn",
    fill = "right"
  )

species_to_path2 <- species_to_path %>% gather(Num, SingleSpecies, -c(Spec,ID), na.rm = TRUE) %>% dplyr::select(-Num)

species_to_path3 <- distinct(species_to_path2)

write_csv2(species_to_path3, "All-pathways-of-MetaCyc_as_common_names_species_noHTML.txt")



all_reagents <- react_to_path2 
colnames(all_reagents) <- c("Pathway", "AllReactions", "ID", "reaction_name")

all_reagents2 <- distinct(all_reagents)
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("→", "->", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("←", "<-", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("↔", "<-->", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("ω =", "omega", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("omega =", "omega", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("β", "beta", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("α", "alpha", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("γ", "gamma", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("δ", "delta", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("ω", "omega", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("a \\[", "\\[", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("μ", "mu", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("ε", "epsilon", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("Δ", "Delta", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("λ", "lambda", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("ν", "nu", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("κ", "kappa", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("ι", "iota", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("χ", "chi", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("1<-->1","", x)}))

all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("di-trans", "ditrans", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("octa-cis", "octacis", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[extracellular space\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[periplasm\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[cytosol\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[in\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[out\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[membrane\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[inner membrane\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[unknown space\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[side 2\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[side 1\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[mitochondrial lumen\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[chloroplast stroma\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[chloroplast thylakoid lumen\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[thylakoid membrane\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[chloroplast stroma\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("\\[Golgi lumen\\]", "", x)}))
all_reagents2 <- data.frame(lapply(all_reagents2, function(x) {gsub("  ", " ", x)}))

##additional corretions may be needed depending on the dataset
#

## Specific corrections for reaction with incoherence beteen arrow direction and direction reported in the file
## The list is obtained froma separate code
## we will use arrow direction in the next step to decide reversibility/left-o-rigth or rtl direction
direction_to_change <- read_csv2("your_path/Reactions_to_change_direction.csv") 
direction_to_change <- direction_to_change[direction_to_change$Direction != "REVERSIBLE",]

all_reagents2$reaction_name[all_reagents2$reaction_name %in% direction_to_change$Reaction] <- gsub(" <--> ", " <- ", all_reagents2[all_reagents2$reaction_name %in% direction_to_change$Reaction,]$reaction_name)
all_reagents2$reaction_name[all_reagents2$reaction_name %in% direction_to_change$Reaction] <- gsub(" -> ", " <- ", all_reagents2[all_reagents2$reaction_name %in% direction_to_change$Reaction,]$reaction_name)




react_to_path_reagents <- all_reagents2 %>%
  separate(
    reaction_name,
    into = c("left", "right"),
    sep = " <--> | <- | -> | = ",
    remove = FALSE,
    convert = FALSE,
    extra = "warn"
  )

react_to_path_reagentsX <- react_to_path_reagents[complete.cases(react_to_path_reagents[ , 5:6]),]

react_to_path_reagents2 <- react_to_path_reagentsX %>%
  mutate(left = strsplit(left, "[[:space:]]\\+[[:space:]]"))
react_to_path_reagents2 <- react_to_path_reagents2 %>%
  mutate(right = strsplit(right, "[[:space:]]\\+[[:space:]]"))


### Reversibile here is = left =-1
react_to_path_reagents3 <- react_to_path_reagents2 %>%
  mutate(score_left = case_when(
    str_detect(reaction_name, " ->") ~ -1,
    str_detect(reaction_name, "<- ") ~ 1,
    str_detect(reaction_name, "<-->") ~ -1,   
    str_detect(reaction_name, "=") ~ -1  
  ))

react_to_path_reagents3 <- react_to_path_reagents3 %>% mutate(score_right = (-1) * score_left)


reagents_left <- react_to_path_reagents3 %>%
  mutate(reagent = left, score = score_left) %>%
  dplyr::select(Pathway, ID, reaction_name, reagent, score) 

reagents_left2 <- unnest(reagents_left, cols = c(as.character("reagent")))

reagents_right <- react_to_path_reagents3 %>%
  mutate(reagent = right, score = score_right) %>%
  dplyr::select(Pathway, ID, reaction_name, reagent, score)

reagents_right2 <- unnest(reagents_right, cols = c(reagent))

# Bind the rows of the two tables to have all the reagents and their score
all_reagents3 <- rbind(reagents_left2, reagents_right2)
all_reagents3 <- mutate(all_reagents3, reagent = as.character(reagent))

all_reagents4 <- distinct(all_reagents3)

################################# 
# Consider the rows that have a number followed by a n or a space and multiplying the score for that number
all_reagents4$score[grepl("^[[:digit:]]+[[:space:]]", all_reagents4$reagent)] <- 
  gsub("^([0-9]+).*", "\\1", all_reagents4$reagent[grepl("^[[:digit:]]+[[:space:]]", all_reagents4$reagent)]) %>%
  as.numeric() *             all_reagents4$score  [grepl("^[[:digit:]]+[[:space:]]", all_reagents4$reagent)]

all_reagents4$score[grepl("^[[:digit:]]+n[[:space:]]", all_reagents4$reagent)] <- 
  gsub("^([0-9]+).*", "\\1", all_reagents4$reagent[grepl("^[[:digit:]]+n[[:space:]]", all_reagents4$reagent)]) %>%
  as.numeric() *             all_reagents4$score  [grepl("^[[:digit:]]+n[[:space:]]", all_reagents4$reagent)]

# remove the stoichiometric number from the name of the reagent
all_reagents4$reagent[grepl("^[[:digit:]]+[[:space:]]", all_reagents4$reagent)] <- 
  str_replace_all(all_reagents4$reagent[grepl("^[[:digit:]]+[[:space:]]", all_reagents4$reagent)],"^([0-9]+)", "")

# remove the stoichiometric number from the name of the reagent
all_reagents4$reagent[grepl("^[[:digit:]]+n[[:space:]]", all_reagents4$reagent)] <- 
  str_replace_all(all_reagents4$reagent[grepl("^[[:digit:]]+n[[:space:]]", all_reagents4$reagent)],"^([0-9]+)", "")


all_reagents5 <- data.frame(lapply(all_reagents4, function(x) {gsub("^[[:space:]]+", "", x)}))
all_reagents5 <- data.frame(lapply(all_reagents5, function(x) {gsub("[[:space:]]$", "", x)}))

all_reagents5 <- distinct(all_reagents5)

write_csv2(all_reagents5, "Pathway_reactions_score.csv")

##consider additional substitution e.g.
all_reagents6 <- all_reagents5[!grepl("unknown", all_reagents5$reagent),]
all_reagents6 <- all_reagents6[!grepl("^n an", all_reagents6$reagent),]

all_reagents6$reagent <-  gsub("n-1 1,5-anhydro-D-fructose", "1,5-anhydro-D-fructose", all_reagents6$reagent)
all_reagents6$reagent <-  gsub("n-1 CMP", "CMP", all_reagents6$reagent)
all_reagents6$reagent <-  gsub("n-1 CDP-ribitol", "CDP-ribitol", all_reagents6$reagent)
all_reagents6$reagent <-  gsub("n-1 beta-1,4-D-mannobiose", "beta-1,4-D-mannobiose", all_reagents6$reagent)
all_reagents6$reagent <-  gsub("n-1 ditrans,octacis-undecaprenyl diphosphate", "ditrans,octacis-undecaprenyl diphosphate", all_reagents6$reagent)
all_reagents6$reagent <-  gsub("^n a ", "", all_reagents6$reagent)


all_reagents6$reagent <-  gsub(" \\(P. gingivalis\\)", "", all_reagents6$reagent)                                                                                            
all_reagents6$reagent <-  gsub(" \\(H. pylori\\)", "", all_reagents6$reagent)                                                                                                
all_reagents6$reagent <-  gsub(" \\(P. putida\\)", "", all_reagents6$reagent)                                                                                                
all_reagents6$reagent <-  gsub(" \\(Brucella\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(B. subtilis 168\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(E. coli\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\[lipid IVA\\]", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(S. aureus\\)", "", all_reagents6$reagent) 
all_reagents6$reagent <-  gsub(" \\(Bucella\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(E. faecium, tetrapeptide\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(Salmonella\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(E. coli K-12\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(E. coli K-12 core type\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(P. gingivalis\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(bacteria and plants\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(Enterococcus faecium\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(Chlamydia\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(bacteria\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(fungi\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(mammalian\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(yeast\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(mycobacteria\\)", "", all_reagents6$reagent)
all_reagents6$reagent <-  gsub(" \\(E. faeciums\\)", "", all_reagents6$reagent)


all_reagents6$reaction_name <-  gsub(" \\(P. gingivalis\\)", "", all_reagents6$reaction_name)                                                                                            
all_reagents6$reaction_name <-  gsub(" \\(H. pylori\\)", "", all_reagents6$reaction_name)                                                                                                
all_reagents6$reaction_name <-  gsub(" \\(P. putida\\)", "", all_reagents6$reaction_name)                                                                                                
all_reagents6$reaction_name <-  gsub(" \\(Brucella\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(B. subtilis 168\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(E. coli\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\[lipid IVA\\]", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(S. aureus\\)", "", all_reagents6$reaction_name) 
all_reagents6$reaction_name <-  gsub(" \\(Bucella\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(E. faecium, tetrapeptide\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(Salmonella\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(E. coli K-12\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(E. coli K-12 core type\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(P. gingivalis\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(bacteria and plants\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(Enterococcus faecium\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(Chlamydia\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(bacteria\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(fungi\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(mammalian\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(yeast\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(mycobacteria\\)", "", all_reagents6$reaction_name)
all_reagents6$reaction_name <-  gsub(" \\(E. faeciums\\)", "", all_reagents6$reaction_name)

all_reagents6 <- distinct(all_reagents6)

all_reagents6$reagent <-  gsub("^a ", "", all_reagents6$reagent)
length(unique(all_reagents6$reagent))

write_csv2(all_reagents6, "Pathway_reactions_score_corrected.csv")


#################
## Algorithm and code 
## by Dr. Fumagalli Maria Rita and Stella Maria Saro
## 28th Feb. 2022
## contact  mariarita.fumagalli@gmail.com
#################