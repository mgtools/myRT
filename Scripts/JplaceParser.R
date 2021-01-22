# Last update Dec 14 2020 by Fatemeh Sharifi

library("treeio")
#library("treeio",lib.loc="~/r-packages")
library("castor")
#library("castor",lib.loc="~/r-packages")
#library("phangorn")
#library("phangorn",lib.loc="~/r-packages")
args = commandArgs(trailingOnly=TRUE)
jp <- read.jplace(args[1])
p <- get.placements(jp, by="all")
# you can use tog (pplacer) to create the placement tree from .jplace, but we directly parse the .jplace placement
#tre <- read.tree(args[2])
sp <- p[order(-p$like_weight_ratio),]
for (n in unique(p$name)) {
entry <- sp[which(sp$name==n)[1],]
#print (entry)
  if (length (which(sp$name==n)) > 1) {
    if ( (sp[which(sp$name==n)[1],]$like_weight_ratio - sp[which(sp$name==n)[2],]$like_weight_ratio) < 0.25 )
    {
      confidence <- 0
    } else {
      confidence <- sp[which(sp$name==n)[1],]$like_weight_ratio
    }
    }
  else {
    confidence <- sp[which(sp$name==n)[1],]$like_weight_ratio
  }
#print (n)
#print(confidence)
node <- entry$node
edgenum <- entry$edge_num
if (node > length(jp@phylo$tip.label)){
  subT <- get_subtree_at_node(jp@phylo,node-length(jp@phylo$tip.label))$subtree
  # Each node is written as a number followed by _ and then RT class, here we are only interested in the class:
  suffixes <- gsub("([0-9]*)_([A-Z|a-z]+[0-9]?[0-9]?[A-Z|a-z]?[0-9]?)$", "\\2", subT$tip.label)
  predict <- names(which.max(table(suffixes)))
  confidence <- 0
  place="non-leaf" 
} else {
  #subT <- get_subtree_at_node(jp@phylo,Ancestors(jp@phylo,node)[1]-length(jp@phylo$tip.label))$subtree
  predict <- gsub("([0-9]*)_([A-Z|a-z]+[0-9]?[0-9]?[A-Z|a-z]?[0-9]?)$", "\\2", jp@phylo$tip.label[node])
  place="leaf"
}
cat(sprintf("%s RVT-%s %2.2f {%d} %s\n", n, predict, confidence, edgenum,place))
#if (length(Ancestors(tre,which(tre$tip.label==node))) > 1 ) {
# plot(get_subtree_at_node(tre,Ancestors(tre,which(tre$tip.label==n))[2]-length(tre$tip.label))$subtree)
#} else {
#  plot(get_subtree_at_node(tre,Ancestors(tre,which(tre$tip.label==n))[1]-length(tre$tip.label))$subtree)
#}
}
