#######################################################
#                Ortholog Analysis                    #
#######################################################

#Test Protein: HMGA
#Test Family: Poaceae


#Load important libraries
library("biomaRt")
library("taxize")
library(tidyverse)
library(readr)
library(dplyr)
library(msa)
library(seqinr)
library(ape)

#######################################################
listMarts(host="plants.ensembl.org")
ensemblp <- useMart(biomart="plants_mart",host="plants.ensembl.org")
Pl_Datasets <- listDatasets(ensemblp)
write.csv(Pl_Datasets, "Pl_Datasets.csv", row.names=FALSE) 

#Copy the "Pl_Datasets.csv" to a new file "Pl_taxon.csv"
#Delete the column "version" and manually edit the column 
#"description" to keep the species name only
#rename the "description" column as "species"

#######################################################
#Obtain the taxonomic details of all the datasets
#Now we have two information #1: Datasets, #2: Species
#We need "Family", "genus"
#taxize package wraps APIs for a large suite of taxonomic 
#databases availab on the web.
#make a vector of species
#make a vector of dataset
df <- read.csv("Pl_taxon.csv")
species <- df$species
dataset <-df$dataset
df1 <- tax_name(query = species, get = c("family", "genus"), db = "ncbi")
#Change content of "query" to that of "species" vector
#Change colname "query" to "species"
#Change colname "db" to "database"
df1[["query"]] <- species                             
colnames(df1)[colnames(df1) == "query"] <- "species"  
colnames(df1)[colnames(df1) == "db"] <- "database"    

df1["dataset"] <-dataset      # Add column dataset
write.csv(df1, "Pl_taxon_table.csv", row.names = FALSE)

#Change the order of column
#The "Pl_taxon_table.csv" has following order of columns
# "database", "species", "family", "genus", "dataset"
# Chnage the order of table as "dataset", "species", "family", "genus", "database" 

df2 <- df1[, c(5, 2, 3, 4, 1)] #Note the first comma means keep all the rows
write.csv(df2, "Pl_taxon_table.csv", row.names = FALSE)

#Subset each family       
#The Ensemblp is represented by 24 families

families <- c("Actinidiaceae", "Amborellaceae", "Apiaceae", "Asteraceae", "Bathycoccaceae", "Brassicaceae", 
              "Chenopodiaceae", "Chlamydomonadaceae", "Cyanidiaceae", "Cyanidiaceae", "Dioscoreaceae", 
              "Fabaceae", "Funariaceae", "Gigartinaceae", "Malvaceae", "Musaceae", "Pleosporaceae", 
              "Poaceae", "Ranunculaceae", "Rosaceae", "Salicaceae", "Selaginellaceae", "Solanaceae", "Vitaceae")
for (i in 1:length(families)) {
  tmp <- df2 [grep(families[i], df2$family),]
  out_file <- paste0(families[i], ".csv")
  write.csv(tmp, out_file, row.names=FALSE)
}             
#######################################################
#Test protein is HMGA ("IPR000116", "IPR031059")
#Identify the HMGAs from reference species Rice (as it is most well-studied grass)
#HMGA From Rice  (O. sativa)

valuesIDs = c("IPR000116", "IPR031059")  #Interpro IDs for HMGAs
os = useDataset("osativa_eg_gene", mart = ensemblp)  #HMGA IDs for O sativa japonica

df <- getBM(attributes =c("ensembl_peptide_id", "refseq_peptide", "interpro", "peptide"),
            filters = ("interpro"), 
            values = valuesIDs, 
            mart = os)  

write.csv(df, "o_sativa_HMGAs.csv", row.names=FALSE)

#The order of the columns are as follows: interpro,	peptide,	ensembl_transcript_id, refseq_peptide
## reorder by column name
df <- df[c("ensembl_peptide_id","refseq_peptide", "interpro", "peptide")] 
write.csv(df, "o_sativa_HMGAs.csv", row.names=FALSE)
##Export the table in fasta
fasta <- read.csv("o_sativa_HMGAs.csv")
fasta.fas<-paste0(">",fasta$ensembl_peptide_id, "|", 
                  fasta$refseq_peptide, "|",
                  fasta$interpro, "|HMGA",
                  "\n",fasta$peptide,"\n")
writeLines(fasta.fas,"o_sativa_HMGAs.fasta")
#export the fasta file as a txt file
write.table(fasta.fas, "o_sativa_HMGAs.txt", row.names=FALSE, col.names = FALSE, quote = FALSE)

#######################################################
# Homologs of Rice HMGAs from grasses

#osativa has three peptide IDs for HMGAs [Os05t0597200-01, Os08t0428800-00, Os09t0402100-01]
#Of these Os05t0597200-01 appears to have partial sequence (lacks initial Methionine) and
#also it does not have any corresponding refseq_peptide_id
#therefore, in further analysis Os05t0597200-01 will not be included

ID1 = "Os08t0428800-00"
ID2 = "Os09t0402100-01"

#Vector of datasets in Poaceae family
df <- read.csv("Poaceae.csv")
dataset <-df$dataset
species<-df$species

#Remove "osativa_eg_gene" from vector as I am looking for osativa homolog 
#and therefore, it will not be listed as an attribute in osativa dataset
dataset #"osativa_eg_gene" is listed in 15th position
dataset1 <- dataset[-15]
dataset1
#The string in datset looks like - "zmays_eg_gene"
#The name for homolog looks like - "zmays_eg_homolog_ensembl_peptide" 
#Therefore each elemnet of the vector need to be modified accordingly

dataset2 <-str_remove(dataset1, "gene") #remove "_gene" from the elements of the vector
homolog <- paste(dataset2, "homolog_ensembl_peptide", sep="")
homolog

#There are 9 oryza species in Poaceae family, 
#therefore to avoid any biasness in favor of rice,
#all the oryza species were removed from further analysis

attributes1 <- c("ensembl_peptide_id", "atauschii_eg_homolog_ensembl_peptide", 
                 "bdistachyon_eg_homolog_ensembl_peptide", "hvulgare_eg_homolog_ensembl_peptide",  
                 "lperrieri_eg_homolog_ensembl_peptide") 

attributes2 <- c("ensembl_peptide_id", "sbicolor_eg_homolog_ensembl_peptide", 
                 "sitalica_eg_homolog_ensembl_peptide", "taestivum_eg_homolog_ensembl_peptide", 
                 "zmays_eg_homolog_ensembl_peptide")


IDs <- c("Os08t0428800-00", "Os09t0402100-01")

#for attributes 1

for (i in seq_along(IDs)) {
  df <- getBM(attributes =attributes1,
              filters = ("ensembl_peptide_id"), 
              values = IDs[i], 
              mart = os)
  
  out_file <- paste0("HMGA_Poaceae_1_", IDs[i], ".csv")
  write.csv(df, out_file, row.names=FALSE)
}

#for attributes 2 

for (i in seq_along(IDs)) {
  df <- getBM(attributes =attributes2,
              filters = ("ensembl_peptide_id"), 
              values = IDs[i], 
              mart = os)
  
  out_file <- paste0("HMGA_Poaceae_2_", IDs[i], ".csv")
  write.csv(df, out_file, row.names=FALSE)
}

#Create a df of HMGA_ID1 and attributes 1
#make vector of ensembl_peptide IDs of HMGAs from each species for OS_HMGA_ID1 homolog  

df <-read.csv("HMGA_Poaceae_1_Os08t0428800-00.csv")
HMGA_ID1 <-unname(unlist(df[1,]))  #unname will not include the corresponding column name
HMGA_ID1
species <-colnames(df) #make a vector of columnnames
species <-str_remove(species, "_eg_homolog_ensembl_peptide") #remove "_eg_homolog_ensembl_peptide" from the elements

#the species names looks  like-
species <- c("ensembl_peptide_id", "atauschii", "bdistachyon", "hvulgare", "lperrieri")

#Replace these with actual names
species <- c("Oryza sativa", "Aegilops tauschii", "Brachypodium distachyon", "Hordeum vulgare", "Leersia perrieri")

df <- data.frame(species, HMGA_ID1) #create dataframe
write.csv(df, "HMGA_Poaceae_1_ID1.csv", row.names = FALSE)

#Create a df of HMGA_ID1 and attributes 2

##read data to df
df <-read.csv("HMGA_Poaceae_2_Os08t0428800-00.csv")

# print unique ids for each column.
write.table(list("species","HMGA_ID1"), "Test.csv",sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE)

for (i in 1:ncol(df)){
  eid_list=paste(unique(df[,i]),collapse=",")
  pd=list(gsub("_.*","",colnames(df)[i]),eid_list)
  #print(paste(pd,collapse="\t"), quote=FALSE)
  write.table(pd, "Test.csv",sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, append = T)
}

#Import the newly created csv file into Excel as described in 
#link (https://www.kotive.com/how-to/run-taskflows/import-a-csv-into-microsoft-excel-or-google-sheets-or-apple-numbers/)
#Manually delete first row as it is for osativa
#save it as HMGA_Poaceae_2_ID1.csv

df <-read.csv("HMGA_Poaceae_2_ID1.csv")
# remove row 1 as it is already there in file HMGA_Poaceae_2_ID1.csv
# remove row 3 as sitalica shows NA

df1 <- df[-1, ]
df2 <- df1[-2, ]
df2

species <- c("Sorghum bicolor", "Triticum aestivum", "Zea mays")

df2[, "species"] <- species
write.csv(df2, "HMGA_Poaceae_2_ID1.csv", row.names = FALSE)

#######
##Merge the two csv "HMGA_Poaceae_1_ID1" and "HMGA_Poaceae_2_ID1"

df1 <- read.csv("HMGA_Poaceae_1_ID1.csv")
df2 <- read.csv("HMGA_Poaceae_2_ID1.csv")

df3 = rbind(df1,df2)
write.csv(df3, "HMGA_Poacea_ID1.csv", row.names = FALSE)

#HMGA_ID2    
#Create a df of HMGA_ID2 and attributes 1
#make vector of ensembl_peptide IDs of HMGAs from each species for OS_HMGA_ID1 homolog  

df <-read.csv("HMGA_Poaceae_1_Os09t0402100-01.csv")
HMGA_ID2 <-unname(unlist(df[1,]))  #unname will not include the corresponding column name
HMGA_ID2
species <-colnames(df) #make a vector of columnnames
species <-str_remove(species, "_eg_homolog_ensembl_peptide") #remove "_eg_homolog_ensembl_peptide" from the elements
species
#the species names looks  like-
species <- c("ensembl_peptide_id", "atauschii", "bdistachyon", "hvulgare", "lperrieri")

#Replace these with actual names
species <- c("Oryza sativa", "Aegilops tauschii", "Brachypodium distachyon", "Hordeum vulgare", "Leersia perrieri")

df <- data.frame(species, HMGA_ID2) #create dataframe

write.csv(df, "HMGA_Poaceae_1_ID2.csv", row.names = FALSE)

#Create a df of HMGA_ID2 and attributes 2
##read data to df
df <-read.csv("HMGA_Poaceae_2_Os09t0402100-01.csv")

# print unique ids for each column.
write.table(list("species","HMGA_ID2"), "Test1.csv",sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE)

for (i in 1:ncol(df)){
  eid_list=paste(unique(df[,i]),collapse=",")
  pd=list(gsub("_.*","",colnames(df)[i]),eid_list)
  #print(paste(pd,collapse="\t"), quote=FALSE)
  write.table(pd, "Test1.csv",sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, append = T)
}

#Import the newly created csv file into Excel as described in 
#link (https://www.kotive.com/how-to/run-taskflows/import-a-csv-into-microsoft-excel-or-google-sheets-or-apple-numbers/)
#Manually delete first row as it is for osativa
#save it as HMGA_Poacea_2_ID2.csv

df <-read.csv("HMGA_Poaceae_2_ID2.csv")
df <- df[-1, ]

species <- c("Sorghum bicolor", "Setaria italica", 
             "Triticum aestivum", "Zea mays")

df[, "species"] <- species
write.csv(df, "HMGA_Poaceae_2_ID2.csv", row.names = FALSE)

#Merge the two csv "HMGA_Poaceae_1_ID2" and "HMGA_Poaceae_2_ID2"

df1 <- read.csv("HMGA_Poaceae_1_ID2.csv")
df2 <- read.csv("HMGA_Poaceae_2_ID2.csv")

df3 = rbind(df1,df2)
write.csv(df3, "HMGA_Poaceae_ID2.csv", row.names = FALSE)

#######################################################
#     Obtain sequence from HMGA1_ID1       

#What we need species name, ensembl ID, and dataset name
#the dataset corresponding to each species is given in Pl_taxon_table

df <-read.csv("HMGA_Poaceae_ID1.csv")
#extract species name 
species_poaceae <- df$species
species_poaceae
#extract the dataset names of these species from Pl_taxon_table
df <-read.csv("Pl_taxon_table.csv")

df1 <- which(df$species %in% species_poaceae)
df1
df1 #The rows are 4,  7, 23, 25, 42, 46, 47, 51, 59

#the column we want to keep in df are column 1 (dataset) and column 2 (species) of Pl_taxon_table

df2 <- df[c(4,  7, 23, 25, 41, 46, 51, 59), c(1, 2)]
df2

#The order of sepecies in HMGA_Poaceae_ID1.csv is alphabetical
#So reorder it 

df2 <- df2[order(df2$species),]
write.csv(df2, "Poaceae_spec_mart.csv", row.names = FALSE)

df <- read.csv("Poaceae_spec_mart.csv")
dataset_poaceae <- df$dataset
species_poaceae <- df$species
dataset_poaceae
species_poaceae

#As the "HMGA_Poaceae_ID1.csv" has one entry as 'NA', the corresponding row (row 7) need to be deleted
#so that the for loop in the sequence extraction don't give empty cells after setaria italica

df <- read.csv("HMGA_Poaceae_ID1.csv")

#Order the HMGA IDs based on alphabetical order of species
df1 <- df[order(df$species),]
df1
write.csv(df1, "HMGA_Poaceae_ID1_order.csv", row.names = FALSE)
df <- read.csv("HMGA_Poaceae_ID1_order.csv")
HMGA_ID1 <- df$HMGA_ID1
HMGA_ID1

#Now setup a loop for extracting sequence

for (i in seq_along(dataset_poaceae)){
  ensemblp <- useMart(biomart="plants_mart",host="plants.ensembl.org")
  poa_mart <- useDataset(dataset_poaceae[i], mart = ensemblp)
  df <- getBM(attributes =c("ensembl_peptide_id", "peptide"),
              filters = ("ensembl_peptide_id"), 
              values = HMGA_ID1[i], 
              mart = poa_mart)
  
  df$species <- species_poaceae[i] #to add a column for species name
  species <- tolower(gsub(".","", species_poaceae[i], fixed=TRUE))
  out_file <- paste0("Poaceae_HMGA_ID1_", species, ".csv")
  write.csv(df, out_file, row.names=FALSE)
}

#Obtain sequence from HMGA1_ID2       

df <-read.csv("HMGA_Poaceae_ID2.csv")
#extract species name 
species_poaceae <- df$species
species_poaceae

#extract the dataset names of these species from Pl_taxon_table

df <-read.csv("Pl_taxon_table.csv")

df1 <- which(df$species %in% species_poaceae)
df1
df1 #The rows are 4,  7, 23, 25, 42, 46, 47, 51, 59

#the column we want to keep in df are column 1 (dataset) and column 2 (species) of Pl_taxon_table

df2 <- df[c(4,  7, 23, 25, 41, 46, 47, 51, 59), c(1, 2)]
df2

#The order of sepecies in HMGA_Poaceae_ID1.csv is alphabetical
#So reorder it 

df2 <- df2[order(df2$species),]
write.csv(df2, "Poaceae_spec_mart_2.csv", row.names = FALSE)

df <- read.csv("Poaceae_spec_mart_2.csv")
dataset_poaceae <- df$dataset
species_poaceae <- df$species
dataset_poaceae
species_poaceae

df <- read.csv("HMGA_Poaceae_ID2.csv")

#Order the HMGA IDs based on alphabetical order of species
df1 <- df[order(df$species),]
df1
write.csv(df1, "HMGA_Poaceae_ID2_order.csv", row.names = FALSE)
df <- read.csv("HMGA_Poaceae_ID2_order.csv")
HMGA_ID2 <- df$HMGA_ID2
HMGA_ID2

#Now setup a loop for extracting sequence

for (i in seq_along(dataset_poaceae)){
  ensemblp <- useMart(biomart="plants_mart",host="plants.ensembl.org")
  poa_mart <- useDataset(dataset_poaceae[i], mart = ensemblp)
  df <- getBM(attributes =c("ensembl_peptide_id", "peptide"),
              filters = ("ensembl_peptide_id"), 
              values = HMGA_ID2[i], 
              mart = poa_mart)
  
  df$species <- species_poaceae[i] #to add a column for species name
  species <- tolower(gsub(".","", species_poaceae[i], fixed=TRUE))
  out_file <- paste0("Poaceae_HMGA_ID2_", species, ".csv")
  write.csv(df, out_file, row.names=FALSE)
}

#Merge all HMGA1_ID1 and convert to fasta
#Merge all the plants HMGA csv files into one

getwd() #to check the current directory
setwd('Poaceae_HMGA_ID1/')
list.files()
df <- list.files(full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(df, "Poaceae_HMGA_ID1_All.csv", row.names = FALSE)

fasta <- read.csv("Poaceae_HMGA_ID1_All.csv")
fasta.fas<-paste0(">",fasta$species,  
                  "|HMGA_ID1", "|",
                  fasta$ensembl_peptide_id,"\n",fasta$peptide,"\n")

writeLines(fasta.fas,"Poaceae_HMGA_ID1_All.fasta")

#Merge all HMGA1_ID2 and convert to fasta       

getwd() #to check the current directory
setwd("..") #move one directory up
getwd()
setwd('Poaceae_HMGA_ID2/')
list.files()
df <- list.files(full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 
write.csv(df, "Poaceae_HMGA_ID2_All.csv", row.names = FALSE)

fasta <- read.csv("Poaceae_HMGA_ID2_All.csv")
fasta.fas<-paste0(">",fasta$species,  
                  "|HMGA_ID2", "|",
                  fasta$ensembl_peptide_id,"\n",fasta$peptide,"\n")

writeLines(fasta.fas,"Poaceae_HMGA_ID2_All.fasta")



#######################################################
#Multiple Sequence Alignment
#Upload sequence file for HMGA1
Poaceae_HMGA1_All <- readAAStringSet("Poaceae_HMGA_1_All.fasta")
Poaceae_HMGA1_All

#Now that we have loaded the sequences, we can run the msa() function 
#which, by default, runs ClustalW with default parameters:

Poaceae_HMGA1_ClustalW <- msa(Poaceae_HMGA1_All, "ClustalW")
print(Poaceae_HMGA1_ClustalW, show="complete")

#If we need Muscle alignment instead of ClustalW
Poaceae_HMGA1_Muscle <- msa(Poaceae_HMGA1_All, "Muscle")
Poaceae_HMGA1_Muscle
print(Poaceae_HMGA1_Muscle, show="complete")

#Upload sequence file for HMGA2
Poaceae_HMGA2_All <- readAAStringSet("Poaceae_HMGA2_All.fasta")
Poaceae_HMGA2_ClustalW <- msa(Poaceae_HMGA2_All, "ClustalW")
print(Poaceae_HMGA2_ClustalW, show="complete")

Poaceae_HMGA2_Muscle <- msa(Poaceae_HMGA2_All, "Muscle")
Poaceae_HMGA2_Muscle
print(Poaceae_HMGA2_Muscle, show="complete")

#######################################################
# Computing distance matrix and Constructing phylogenetic tree for HMGA1 

#Now we have two alignment files for HMGA1 (Poaceae_HMGA1_ClustalW, and Poaceae_HMGA1_Muscle)
#Convert the alignment for later processing using seqinr package:

library(seqinr)
HMGA1_ClustalW_Aln2 <- msaConvert(Poaceae_HMGA1_ClustalW, type="seqinr::alignment")
HMGA1_Muscle_Aln2 <- msaConvert(Poaceae_HMGA1_Muscle, type="seqinr::alignment")

#Now we compute a distance matrix using the dist.alignment() function

d <- dist.alignment(HMGA1_ClustalW_Aln2, "identity")
d1 <-dist.alignment(HMGA1_Muscle_Aln2, "identity")

as.matrix(d)[1:9, "O.sativa_Os08t0428800-00", drop=FALSE]
as.matrix(d1)[1:9, "O.sativa_Os08t0428800-00", drop=FALSE]

#Now we can construct a phylogenetic tree with the neighbor joining algorithm 
#using the nj() function from the ape package:
HMGA1_Tree_ClustalW <- nj(d)
HMGA1_Tree_Muscle <- nj(d1)
plot(HMGA1_Tree_ClustalW, main="Phylogenetic Tree of Poaceae HMGA1 Sequences")
plot(HMGA1_Tree_Muscle, main="Phylogenetic Tree of Poaceae HMGA1 Sequences")

# Computing distance matrix and Constructing phylogenetic tree for HMGA2 
HMGA2_ClustalW_Aln2 <- msaConvert(Poaceae_HMGA2_ClustalW, type="seqinr::alignment")
HMGA2_Muscle_Aln2 <- msaConvert(Poaceae_HMGA2_Muscle, type="seqinr::alignment")

d <- dist.alignment(HMGA2_ClustalW_Aln2, "identity")
d1 <-dist.alignment(HMGA2_Muscle_Aln2, "identity")

as.matrix(d)[1:11, "O.sativa", drop=FALSE]
as.matrix(d1)[1:11, "O.sativa", drop=FALSE]

#Now we can construct a phylogenetic tree with the neighbor joining algorithm 
HMGA2_Tree_ClustalW <- nj(d)
HMGA2_Tree_Muscle <- nj(d1)

plot(HMGA2_Tree_ClustalW, main="Phylogenetic Tree of Poaceae HMGA2 Sequences")
plot(HMGA2_Tree_Muscle, main="Phylogenetic Tree of Poaceae HMGA2 Sequences")


################################################################









