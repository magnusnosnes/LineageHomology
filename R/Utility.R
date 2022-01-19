#   ____________________________________________________________________________
#   Internal functions: nodeheight from phytools and descendants from ips   ####

# convert vector of x to binary matrix
# written by Liam J. Revell 2012
#' Title
#'
#' @param x
#' @param seq
#'
#' @return
#' @export
#'
#' @examples
to.matrix<-function(x,seq){
  X<-matrix(0,length(x),length(seq),dimnames=list(names(x),seq))
  for(i in 1:length(seq)) X[x==seq[i],i]<-1
  return(X)
}




#' Title
#'
#' @param tree
#' @param node
#'
#' @return
#' @export
#'
#' @examples
find_midedge_ancestral = function(tree, node) {

  ind = which(tree$edge[,2]==node)
  if(length(ind)==0) { #Check if node is root.
    #It is the root and we can't shift it back more.
    minus = 0
  }
  else {
    minus = tree$edge.length[ind]/2 #Remove half of the edge.
  }
  minus
}

#' Written by Magnus Nygaard Osnes 2021.
#'
#' @param tree
#' @param taxa
#' @param mrca
#'
#' @return
#' @export
#'
#' @examples
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

#' nodepath quick
#' Written by Magnus Nygaard Osnes 2021.
#'
#' @param tree
#' @param taxa
#' @param mrca
#'
#' @return
#' @export
#'
#' @examples
nodepath_quick_w_mrca_nodenumber = function(tree,taxa,mrca) {

  bifurcations = c()
  mrca_node = mrca
  bifurcations = c(mrca_node)

  tip_nodes = taxa
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


## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-11-07)

#' @param phy
#'
#' @param node
#' @param type
#' @param ignore.tip
#' @param labels
#'
#' @export

descendants <- function(phy, node, type = "t", ignore.tip = TRUE,
                        labels = FALSE){

  # checks and definitions
  # ----------------------
  if ( inherits(phy, "phylo") ){
    edge <- phy$edge
  } else {
    if ( !is.matrix(phy) ) {
      stop("'phy' must be of classes 'phylo' or 'matrix'")
    } else {
      edge <- phy
      labels <- FALSE
    }
  }
  if ( length(node) > 1) stop("'node' must be vector of length 1")
  type <- match.arg(type, c("all", "daughter", "internal", "terminal"))
  tips <- setdiff(edge[, 2], edge[, 1])

  # 'node' is a tip
  # ---------------
  if ( node <= max(tips) ){
    if ( ignore.tip ){
      x <- node
    } else {
      stop("node ", node, " is not an internal node")
    }
  } else {

    # normal procedure when 'node' is internal
    # ----------------------------------------
    x <- edge[edge[,1] == node, 2] # immediate daughter nodes
    if ( type %in% c("internal", "terminal", "all") ){
      repeat{
        xx <- x
        x <- sort(unique(c(x, edge[,2][edge[,1] %in% x])))
        if (identical(x, xx)) break
      }
      if ( type == "internal" ) x <- setdiff(x, tips)
    }
  }
  ## apply 'type' argument:
  ## ----------------------
  if ( type == "terminal" ) {
    x <- intersect(x, tips)
    if (labels) {
      x <- phy$tip.label[x]
    }
  }
  x
}


## function finds the height of a given node
## written by Liam Revell 2014, 2015, 2016
#' Title
#'
#' @param tree
#' @param node
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nodeheight<-function(tree,node,...){
  if(hasArg(root.edge)) root.edge<-list(...)$root.edge
  else root.edge<-FALSE
  if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
  else ROOT<-0
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(node==(Ntip(tree)+1)) h<-0
  else {
    a<-setdiff(c(getAncestors(tree,node),node),Ntip(tree)+1)
    h<-sum(tree$edge.length[sapply(a,function(x,e) which(e==x),e=tree$edge[,2])])
  }
  h+ROOT
}


# returns the heights of each node
# written by Liam J. Revell 2011, 2012, 2013, 2015, 2016
# modified by Klaus Schliep 2017
#' Title
#'
#' @param tree
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nodeHeights<-function(tree,...){
  if(hasArg(root.edge)) root.edge<-list(...)$root.edge
  else root.edge<-FALSE
  if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
  else ROOT<-0
  nHeight <- function(tree){
    tree <- reorder(tree)
    edge <- tree$edge
    el <- tree$edge.length
    res <- numeric(max(tree$edge))
    for(i in seq_len(nrow(edge))) res[edge[i,2]] <- res[edge[i,1]] + el[i]
    res
  }
  nh <- nHeight(tree)
  return(matrix(nh[tree$edge], ncol=2L)+ROOT)
}

## function gets ancestor node numbers, to be used internally by
## written by Liam J. Revell 2014
#' Title
#'
#' @param tree
#' @param node
#' @param type
#'
#' @return
#' @export
#'
#' @examples
getAncestors<-function(tree,node,type=c("all","parent")){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  type<-type[1]
  if(type=="all"){
    aa<-vector()
    rt<-Ntip(tree)+1
    currnode<-node
    while(currnode!=rt){
      currnode<-getAncestors(tree,currnode,"parent")
      aa<-c(aa,currnode)
    }
    return(aa)
  } else if(type=="parent"){
    aa<-tree$edge[which(tree$edge[,2]==node),1]
    return(aa)
  } else stop("do not recognize type")
}

#Written by Magnus Nygaard Osnes 2021.

#' nodepath_quick
#'
#' @param tree
#' @param taxa
#'
#' @return
#' @export
#'
#' @examples
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


#Written by Magnus Nygård Osnes
#' reorder_LH
#' @description
#' Takes the output of LineageHomology and reorders the transmission lineages so that the largest comes first 1, second largest comes 2 etc.
#' @param result_LH #Return from LineageHomology.
#'
#' @return
#' @export
#'
#' @examples
reorder_LH =  function(result_LH) {
    inds = sort.int(result_LH$Lineage_sizes,index.return =T)
    result_LH$Lineage_sizes=result_LH$Lineage_sizes[rev(inds$ix)]
    result_LH$Taxa_names=result_LH$Taxa_names[rev(inds$ix)]
    result_LH$`MRCA's`=result_LH$`MRCA's`[rev(inds$ix)]
    result_LH$lineage_state=result_LH$lineage_state[rev(inds$ix)]
    result_LH$Halfedge_over_tmrca=result_LH$Halfedge_over_tmrca[rev(inds$ix)]
    result_LH
  }
