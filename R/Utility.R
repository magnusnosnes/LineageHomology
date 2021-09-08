

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



#   ____________________________________________________________________________
#   Internal functions: nodeheight from phytools and descendants from ips   ####

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

#' Title
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

#' Title
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

