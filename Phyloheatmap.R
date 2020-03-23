install.packages("phytools")
library(phytools)
library(ape)
library(ggplot2)
library(seqinr)
library(msa)

#Upload sequence file for HMGA1
Poaceae_HMGA1_All <- readAAStringSet("Poaceae_HMGA_1_All.fasta")
Poaceae_HMGA1_All

#Now that we have loaded the sequences, we can run the msa() function 
#which, by default, runs ClustalW with default parameters:

Poaceae_HMGA1_ClustalW <- msa(Poaceae_HMGA1_All, "ClustalW")
Poaceae_HMGA1_ClustalW
#the default printing function shortens the alignment for the sake of compact output
#print() function provided by the msa package provides some ways for customizing the output,
#such as, showing the entire alignment split over multiple blocks of sub-sequences:

print(Poaceae_HMGA1_ClustalW, show="complete")

Poaceae_HMGA1_Muscle <- msa(Poaceae_HMGA1_All, "Muscle")
Poaceae_HMGA1_Muscle
print(Poaceae_HMGA1_Muscle, show="complete")

############################################################
############################################################

#Upload sequence file for HMGA2
Poaceae_HMGA2_All <- readAAStringSet("Poaceae_HMGA2_All.fasta")
Poaceae_HMGA2_All

#Now that we have loaded the sequences, we can run the msa() function 
#which, by default, runs ClustalW with default parameters:

Poaceae_HMGA2_ClustalW <- msa(Poaceae_HMGA2_All, "ClustalW")
Poaceae_HMGA2_ClustalW
print(Poaceae_HMGA2_ClustalW, show="complete")

Poaceae_HMGA2_Muscle <- msa(Poaceae_HMGA2_All, "Muscle")
Poaceae_HMGA2_Muscle
print(Poaceae_HMGA2_Muscle, show="complete")

#############################################################
#############################################################

### Computing distance matrix and Constructing phylogenetic tree for HMGA1 

#Now we have two aln files for HMGA1 (Poaceae_HMGA1_ClustalW, and Poaceae_HMGA1_Muscle)
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

library(ape)

HMGA1_Tree_ClustalW <- nj(d)
HMGA1_Tree_Muscle <- nj(d1)

plot(HMGA1_Tree_ClustalW, main="Phylogenetic Tree of Poaceae HMGA1 Sequences")
plot(HMGA1_Tree_Muscle, main="Phylogenetic Tree of Poaceae HMGA1 Sequences")

phylo.heatmap(HMGA1_Tree_ClustalW,d)

#############################################################
#############################################################

### Computing distance matrix and Constructing phylogenetic tree for HMGA2 

#Now we have two aln files for HMGA1 (Poaceae_HMGA1_ClustalW, and Poaceae_HMGA1_Muscle)
#Convert the alignment for later processing using seqinr package:

library(seqinr)
HMGA2_ClustalW_Aln2 <- msaConvert(Poaceae_HMGA2_ClustalW, type="seqinr::alignment")
HMGA2_Muscle_Aln2 <- msaConvert(Poaceae_HMGA2_Muscle, type="seqinr::alignment")

#Now we compute a distance matrix using the dist.alignment() function

d <- dist.alignment(HMGA2_ClustalW_Aln2, "identity")
d1 <-dist.alignment(HMGA2_Muscle_Aln2, "identity")

as.matrix(d)[1:11, "O.sativa", drop=FALSE]
as.matrix(d1)[1:11, "O.sativa", drop=FALSE]

#Now we can construct a phylogenetic tree with the neighbor joining algorithm 
#using the nj() function from the ape package:

library(ape)

HMGA2_Tree_ClustalW <- nj(d)
HMGA2_Tree_Muscle <- nj(d1)

plot(HMGA2_Tree_ClustalW, main="Phylogenetic Tree of Poaceae HMGA2 Sequences")
plot(HMGA2_Tree_Muscle, main="Phylogenetic Tree of Poaceae HMGA2 Sequences")

#############################################################################

