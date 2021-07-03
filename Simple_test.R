library(LineageHomology)
## basic example code

#Loading other packages for simulating data.
library(ape)
library(phytools)
library(phangorn)
library(BactDating)
library(pbapply)

#Simulate data from Norway and rest of the world
set.seed(150)
tree_test = simdatedtree(nsam=12, dateroot=2000)
Q=matrix(c(0.5,0.5,0.5,0.5), nrow=2,ncol=2, byrow=F)
colnames(Q)=c("Norway","RoW")
loc = c("Norway", "Norway","Norway","RoW", "RoW", "RoW", "Norway", "Norway", "Norway", "RoW","RoW", "Norway")
names(loc) = tree_test$tip.label

#Plot
plot.phylo(tree_test,edge.width = 2,label.offset = 0.2, mar=c(0.2,0.2,0.2,0.2))
axisPhylo(root.time=2000, backward=F, lwd=2)
fit1 = ace(x=loc, phy= tree_test, type="discrete", mod="ARD")
nodelabels(pie=fit1$lik.anc,cex=0.7,piecol=c("Red","Blue"))
tips = to.matrix(loc,seq=c("Norway", "RoW"))
tiplabels(pie=tips, cex=0.7,piecol=c("Red","Blue"))

#Reconstruct ancestral states using ace.

Result = LineageHomology(tree_test, ace_nodes=fit1$lik.anc,
                ace_tips = to.matrix(loc, seq=c("Norway", "RoW")), start_time=2000)

Result_multi = pbreplicate(10, LineageHomology_w_uncertainty(tree_test, ace_nodes=fit1$lik.anc,
                         ace_tips = to.matrix(loc, seq=c("Norway", "RoW")), start_time=2000))

Result$`MRCA's`
Result$Halfedge_over_tmrca

Result$Lineage_sizes





#Test function::

mrca = 13
tree = tree_test
to = c("t10","t3","t8")

nodepath_quick = function(tree,taxa) {

  bifurcations = c()
  mrca_node = getMRCA(tree,tip = taxa)
  bifurcations = c(mrca_node)

  tip_nodes = sapply(taxa,function(x,y) which(y==x),y=tree$tip.label)
  new_nodes = tree$edge[,1][which(tree$edge[,2]%in%tip_nodes)]
  bifurcations = unique(c(new_nodes,bifurcations))
  new_nodes = new_nodes[new_nodes!=mrca_node]
  #If any is mrca, remove them from loop.
  while(length(new_nodes)>0) {
    new_nodes = tree$edge[,1][which(tree$edge[,2]%in%new_nodes)]
    new_nodes = new_nodes[new_nodes!=mrca_node]
    bifurcations = unique(c(new_nodes,bifurcations))
  }
  bifurcations #Return all the bifurcation nodes.
}


nodepath_quick_w_mrca = function(tree,taxa,mrca) {

  bifurcations = c()
  mrca_node = mrca
  bifurcations = c(mrca_node)

  tip_nodes = sapply(taxa,function(x,y) which(y==x),y=tree$tip.label)
  new_nodes = tree$edge[,1][which(tree$edge[,2]%in%tip_nodes)]
  bifurcations = unique(c(new_nodes,bifurcations))
  new_nodes = new_nodes[new_nodes!=mrca_node]
  #If any is mrca, remove them from loop.
  while(length(new_nodes)>0) {
    new_nodes = tree$edge[,1][which(tree$edge[,2]%in%new_nodes)]
    new_nodes = new_nodes[new_nodes!=mrca_node]
    bifurcations = unique(c(new_nodes,bifurcations))
  }
  bifurcations #Return all the bifurcation nodes.
}

nodepath_quick_w_mrca(tree_test, c("t8","t3","t10"),17)
