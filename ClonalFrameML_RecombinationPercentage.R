#Script takes ClonalFrameML output files and gets a % of each sequence that is recombined sites

library(ape)
library(phangorn)
library(tidyverse)

#Change the working directory to point to where the files are
setwd("~/Desktop/Chimaera/RemoveMultiNov2022/clonalframeml_MixedRemoved")

#Change the prefix for the run here
prefix="MC-ClonalFrame-MixedRemoved"
treefile = paste(prefix,".labelled_tree.newick",sep="")
xreffile = paste(prefix,".position_cross_reference.txt",sep="")
istatefile = paste(prefix,".importation_status.txt",sep="")

#get the total number of sites in the alignment
#if you used an xmfa file you need to add up all the = in the file and multiply by 1000 and add that to the alignment length
#e.g. grep -c "^=$" alignment.xmfa
numOfEquals = 3302 #change this to the number of = that separate the genes
alignmentLength = length(scan(xreffile,sep=","))+(numOfEquals*1000)

#read in the status file
status <- read.table(istatefile,h=T,as.is=T,sep="\t")
#read in the tree
tree = read.tree(treefile)

#add a column to the status file which gets the length of the section (column 3-column 2)
status$Len=status$End-status$Beg

#add up all the length entries for a sequence/node (Node column)
combined<- status %>% 
  group_by(Node) %>% 
  summarise(across(starts_with("Len"), ~sum(., na.rm = TRUE)))

#extract the nodes into a different table from the genomes
nodes <- combined %>% 
  filter(grepl("^NODE", Node))
genomes <- combined %>% 
  filter(!grepl("^NODE", Node))

# get a list of all ancestors for each node in the tree
tree.labels <- c(tree$tip.label, tree$node.label)
tree.ancestors <- Ancestors(tree)
names(tree.ancestors) <- tree.labels
tree.ancestors.labelled <- lapply(tree.ancestors, function(ind.vec) tree.labels[ind.vec])

#Go through the genomes and for each one, get all the internal nodes and add their recombination lengths to the genomes total
for (g in 1:nrow(genomes)) {
  #get the name and list of Nodes that are the ancestor of that leaf
  name <- genomes[g,"Node"]%>% pull()
  g_anc <- tree.ancestors.labelled[name]
  g_anc <- as.data.frame(tree.ancestors.labelled[name],check.names=FALSE)
  
  #variable for extra sites to add to genomes total
  toAdd=0
  
  #go through each node and add the length. 
  for (n in 1:nrow(g_anc)) {
    n_name <-g_anc[n,name]
    n_details <- nodes %>% 
      filter(grepl(n_name, Node))
    #Some nodes have no recombination so wont be in the nodes data frame so these have length 0 and should be skipped
    if (length(n_details$Len) > 0) {
      toAdd = toAdd + n_details$Len
    }  
  }
  
  #Add the length to the genome total
  genomes[g,"Len"] = genomes[g,"Len"]+toAdd
  
  
}

#create a percentage column for the genomes which is the percent of total sites that have recombination
genomes$Percentage=round((genomes$Len/alignmentLength)*100,2)

#write the genomes table out to file
write.table(genomes, file=paste(prefix,"recombinationPercentage.txt", sep = "_"),sep ="\t", row.names = FALSE, quote = FALSE)



