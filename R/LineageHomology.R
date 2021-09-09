library(ape)
library(lubridate)


#Pseudocode:
#Fint the parent nodes of all tips.
#For each tip; check for transition to its parent node
#If transition, the algorithm has found a singleton.
#If not, move on a branch up and check any new descendant emerging from the parent node.
#If any descendant is found, this is no longer a singleton but represents a local transmission event.
#If no descendants are found, this is just a singleton with a long branch.


#   ____________________________________________________________________________
#   Main functions                                                          ####


#' LineageHomology
#' @description
#' LineageHomology takes the output of an ancestral state reconstruction method with included state probabilities at
#'  each node and counts transmission lineages (TLs). A transmission lineage is defined as a connected group of tips where state transitions between ancestral and descendant nodes
#'  that have a probability lower than 50 percent.
#' The function also counts the number of tips that are not connected to any other tips in this way (singletons).
#' The method is analogous to that introduced by du Plessis et al. (2021) (DOI: 10.1126/science.abf2946).
#' The input format should match the return format from the ape::ace function (see the help files at ?ape::ace).
#' The outputs contain descriptions of the size of TLs, singletons, and other useful summaries.
#' See the full tutorial including visualization methods at: \url{https://github.com/magnusnosnes/LineageHomology}.
#'
#' @param tree Phylogenetic tree
#' @param ace_nodes Node probabilities in ace format.
#' @param ace_tips Tip probabilities in ace format
#' @param give_tips Constrain the computations to a set of tips, e.g. all the names of the tips in a certain state (e.g. one specific geographical location).
#' @param start_time Date of root in the phylogeny.
#'
#' @return
#' @return Import_LocalTrans
#'  The first entry can be
#'  thought of as an estimate of the number of tips that have acquired their state through novel acquisition (a transition in the phylogeny).
#'  Similarly, the second number describes the number of tips that have inherited their state by inheritance (homology).
#'  Since each transmission lineage started is defined by a state transition,
#'  each transmission lineage of size n is composed of 1 novel acquisition and n-1 homologous acquisitions.
#'  Singletons are directly translated to novel acquisitions.
#' @return Lineage_sizes
#'  contains the number of tips that belong to each TL.
#'  This means that all the internal nodes will be mapped to the same state with 50 percent or higher probability for the entire TL.
#' @return Taxa_names is a list of lists where each sub list contains the names of the tips in each transmission lineage.
#' The order is the same as the listed order of the Lineage_sizes.
#' @return MRCA’s
#' contains the time of the most recent common ancestor of the groups.
#' @return Halfedge_over_tmrca
#' is the time point on the middle of the edge that is ancestral to each TL's mrca.
#' This is used in analyses of importation and local transmission.
#' @return Lineage_state
#' contains the state the transmission lineages belong to,
#' the states are treated as numbers in the analyses, so the relative ordering has to be looked up.
#' @export
#'
#' @importFrom ape nodepath
#' @examples
#' \dontrun{
#' set.seed(400)
#' library(BactDating)
#' tree_test = simdatedtree(nsam=300, dateroot=2000) #300 taxa and date of the root at year 2000.
#' tree_test = ladderize(tree_test) #Reorder the tree to make it look nice
#' fit1 = ace(x=trait, phy= tree_test, type="discrete", mod="ARD") #Estimate the ancestral history
#' LineageHomology(tree_test, ace_nodes=fit1$lik.anc,
#' ace_tips = to.matrix(loc, seq=c("Norway", "RoW")), start_time=2000)
#' }

LineageHomology = function(tree,ace_nodes,ace_tips, give_tips=NA,start_time=NA) {
  #The node number of the tips are
  tip_nodes = sapply(tree$tip.label,function(x,y) which(y==x),y=tree$tip.label)

  if(is.na(give_tips)==F){ #If the count is among given tips.
    tip_nodes = sapply(give_tips,function(x,y) which(y==x),y=tree$tip.label)
  }

  tip_states = ace_tips
  node_states = ace_nodes
  node_draws = apply(ace_nodes, 1,FUN=function(x) which.max(x)) #Draw most probable nodes
  rownames(node_states)=(length(tree$tip.label)+1):((length(tree$tip.label))+nrow(ace_nodes)) #Old specification was buggy due to missing nodes.
  tip_label_list = tree$tip.label
  #Parent and child.
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]


  #Size of lineage
  lineage_size = c()
  #Size of lineage
  lineage_state = c()
  #Tips in lineage
  lineage_tips = list()
  #MRCA node list
  mrca.list = c()
  #Initialize importations,
  importations = 0

  #Initialize
  local_transmissions = 0
  #Initialize index counter
  counter=1



  checked_most_basal = F

  while(length(tip_nodes)>0) {
    current_node = tip_nodes[1]
    tip_label_list = c(names(current_node)) #Set up a list containing tips for the current tip.
    tip_nodes = tip_nodes[-1] #Remove tip from search space immediately.
    current_state = which(tip_states[current_node,]==1) #If other columns than this has more than 50%, there is a transition.
    lineage_state = c(lineage_state,current_state)  #Save the state of the lineage.
    goes_all_the_way_back = F
    tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    ancestor_node = parent[which(child==current_node)]
    ancestor_node_row = which(rownames(node_states)==ancestor_node)
    ancestor_state =node_draws[ancestor_node_row]

    if(current_state!=ancestor_state) {  #If more than one mapping on the branch, then it is estimated as an importation.
      tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    }

    else {
      while(current_state==ancestor_state) { #While there is no transition, move one branch up.

        if(((current_node%in%child)==T)) {
          current_node = ancestor_node #If the loop gets all the way through, update ancestor node.
          ancestor_node = parent[which(child==current_node)]
          if((current_node%in%child)==F) { #break loop if we're already at the root node.
            goes_all_the_way_back = T
            break()
          }
          #Note that current state doesn't need to be updated since we're moving up if it is the same.
          ancestor_node_row = which(rownames(node_states)==ancestor_node)
          ancestor_state = node_draws[ancestor_node_row]
        }

      }
      tip_search_subspace = descendants(tree, current_node)
      subspace=tip_search_subspace[tip_search_subspace%in%tip_nodes] #Find all possible tips in searchspace descending from the current node.

      if(length(subspace)!=0) {
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip.
          nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          internal = nodes_to_check[-length(nodes_to_check)] #The minus removes the tip - maybe unessecary.
          nods= sapply(internal, FUN = function(x) node_draws[which(rownames(node_states)==x)]) #node states
          ##The following are unnecessary as tip always have the same state in this example
          tps = which.max(tip_states[nodes_to_check[length(nodes_to_check)],])
          nods = c(nods, tps)

          if(sum(nods==current_state)==length(nods)) { #If all states are equal to the state of the current node proceed.
            tips_in_segment = tips_in_segment+1
            tip_label_list = c(tip_label_list, names(tip_nodes[which(tip_nodes==nodes_to_check[length(nodes_to_check)])])) #Add the tips that will be removed.
            tip_nodes = tip_nodes[-which(tip_nodes==nodes_to_check[length(nodes_to_check)])] #Remove found tip from search space.
          }
        }
      }
    }

    if(tips_in_segment > 1 ) {
      local_transmissions = local_transmissions+(tips_in_segment-1)
      importations = importations+1
    }
    if(tips_in_segment==1) {
      importations = importations+1
    }
    lineage_size = c(lineage_size,tips_in_segment)
    lineage_tips[[counter]]=tip_label_list
    mrca.list[counter] = current_node
    counter = counter+1 #move one up to fill the next index on in the next loop iteration.
    #Debug: print(paste0("imports are:",importations))
    #Debug: print(paste0("Local transmissions are:",local_transmissions))
  }
  #process the mrca nodes to ages.
  names(lineage_tips)=paste0("Lineage no: ", 1:length(lineage_tips))
  lineage_mrca = unlist(lapply(mrca.list, FUN = function(x) nodeheight(tree,node = x)))
  if(start_time){
    lineage_mrca =lineage_mrca+start_time
  }
  #For transition dates it might be more appropriate to use the midpoint of the edge ancestral to the lineage
  subtract_ancestral_edge_time = unlist(lapply(mrca.list, FUN = function(x) find_midedge_ancestral(tree,node = x)))
  halfedge_over_mrca =lineage_mrca - subtract_ancestral_edge_time

  #Returns
  return(list("Import_LocalTrans"=c(importations, local_transmissions),"Lineage_sizes"=lineage_size,"Taxa_names"=lineage_tips, "MRCA's"=lineage_mrca,"lineage_state"=lineage_state,"Halfedge_over_tmrca"=halfedge_over_mrca))

}


#' LineageHomology_w_uncertainty
#' @description
#'  LineageHomology_w_uncertainty returns the same outputs as LineageHomology.
#'  LineageHomology counts transmission lineages according to state transitions between ancestral and descendant nodes that have a probability higher than 50 percent.
#'  LineageHomolog_w_uncertainty instead samples states from the posterior probability of states for each node.
#'  Thus, the transmission lineage division is probabilistic and retains the uncertainty in the posterior distribution if one runs multiple times,
#'   e.g., using base::replicate() (see examples).
#' @param tree Phylogenetic tree
#' @param ace_nodes Node probabilities in ace format
#' @param ace_tips Tip probabilities in ace format
#' @param give_tips Constrain the computations to a set of tips, e.g. all the names of the tips in a certain state (e.g. one specific geographical location).
#' @param start_time Date of root in the phylogeny.
#'
#' @return See LineageHomology for descriptions of the returns.
#' @export
#'
#' @examples
#' \dontrun{
#' multi_counts = pbreplicate(
#'   100,
#'   LineageHomology::LineageHomology_w_uncertainty(
#'     tree,
#'     ace_nodes = test$lik.anc,
#'     give_tips = Norwegian_tips,
#'     ace_tips = to.matrix(Locations, seq = c("Norway", "RoW")),
#'     start_time = start_date
#'   )
#' )
#' }

LineageHomology_w_uncertainty = function(tree,ace_nodes,ace_tips, give_tips=NA,start_time=NA) {
  #The node number of the tips are
  tip_nodes = sapply(tree$tip.label,function(x,y) which(y==x),y=tree$tip.label)

  if(is.na(give_tips)==F){ #If the count is among given tips.
    tip_nodes = sapply(give_tips,function(x,y) which(y==x),y=tree$tip.label)
  }
  #The algorithm proceeds by removing tips from the tip nodes list until there is none left.

  tip_states = ace_tips
  node_states = ace_nodes
  #Sample nodes from the ASR once here, so they don't have to be resampled multiple times.
  node_draws = apply(ace_nodes, 1,FUN=function(x) sample(1:2,1, prob=x))
  rownames(node_states)=(length(tree$tip.label)+1):((length(tree$tip.label))+nrow(ace_nodes)) #Old specification was buggy due to missing nodes.
  tip_label_list = tree$tip.label
  #Parent and child.
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]


  #Size of lineage
  lineage_size = c()
  #Size of lineage
  lineage_state = c()
  #Tips in lineage
  lineage_tips = list()
  #MRCA node list
  mrca.list = c()
  #Initialize importations,
  importations = 0

  #Initialize
  local_transmissions = 0
  #Initialize index counter
  counter=1



  checked_most_basal = F

  while(length(tip_nodes)>0) {
    current_node = tip_nodes[1]
    #print(current_node)
    tip_label_list = c(names(current_node)) #Set up a list containing tips for the current tip.
    tip_nodes = tip_nodes[-1] #Remove tip from search space immediately.
    current_state = which(tip_states[current_node,]==1) #If other columns than this has more than 50%, there is a transition.
    lineage_state = c(lineage_state,current_state)  #Save the state of the lineage.
    goes_all_the_way_back = F
    tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    ancestor_node = parent[which(child==current_node)]
    ancestor_node_row = which(rownames(node_states)==ancestor_node)
    ancestor_state =node_draws[ancestor_node_row]

    if(current_state!=ancestor_state) {  #If more than one mapping on the branch, then it is estimated as an importation.
      tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    }

    else {
      while(current_state==ancestor_state) { #While there is no transition, move one node up.

        if(((current_node%in%child)==T)) {
          current_node = ancestor_node #Move one node up.
          ancestor_node = parent[which(child==current_node)]
          if((current_node%in%child)==F) { #break loop if we're already at the root node.
            goes_all_the_way_back = T
            #print("Broke")
            break()
          }
          #Note that current state doesn't need to be updated since we're moving up if it is the same.
          ancestor_node_row = which(rownames(node_states)==ancestor_node)
          ancestor_state =node_draws[ancestor_node_row]
        }

      }
      tip_search_subspace = descendants(tree, current_node)
      subspace=tip_search_subspace[tip_search_subspace%in%tip_nodes] #Find all from the set that can be included.


      if(length(subspace)!=0) {
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip this - should become buggy with sampling.
          nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          #nodes_to_check = nodepath_quick_w_mrca(tree,subspace[i],mrca = current_node)
          internal = nodes_to_check[-length(nodes_to_check)]
          nods= sapply(internal, FUN = function(x) node_draws[which(rownames(node_states)==x)]) #node states
          ##The following are unnecessary as tip always have the same state in this example
          tps = which.max(tip_states[nodes_to_check[length(nodes_to_check)],])
          nods = c(nods, tps)

          if(sum(nods==current_state)==length(nods)) { #If all states are equal to the state of the current node proceed.
            tips_in_segment = tips_in_segment+1
            tip_label_list = c(tip_label_list, names(tip_nodes[which(tip_nodes==nodes_to_check[length(nodes_to_check)])])) #Add the tips that will be removed.
            tip_nodes = tip_nodes[-which(tip_nodes==nodes_to_check[length(nodes_to_check)])] #Remove found tip from search space.
          }
        }
      }
    }

    if(tips_in_segment > 1 ) {
      local_transmissions = local_transmissions+(tips_in_segment-1)
      importations = importations+1
    }
    if(tips_in_segment==1) {
      importations = importations+1
    }
    lineage_size = c(lineage_size,tips_in_segment)

    lineage_tips[[counter]]=tip_label_list
    mrca.list[counter] = current_node
    counter = counter+1 #move one up to fill the next index on in the next loop iteration.
    #deubg: print(paste0("imports are:",importations))
    #debug: print(paste0("Local transmissions are:",local_transmissions))
  }
  #process the mrca nodes to ages.
  names(lineage_tips)=paste0("Lineage no: ", 1:length(lineage_tips))
  lineage_mrca = unlist(lapply(mrca.list, FUN = function(x) nodeheight(tree,node = x)))
  if(start_time){
    lineage_mrca =lineage_mrca+start_time
  }

  #For transition dates it might be more appropriate to use the midpoint of the edge ancestral to the lineage
  subtract_ancestral_edge_time = unlist(lapply(mrca.list, FUN = function(x) find_midedge_ancestral(tree,node = x)))
  halfedge_over_mrca =lineage_mrca - subtract_ancestral_edge_time

  return(list("Import_LocalTransmission"=c(importations, local_transmissions),"Lineage_sizes"=lineage_size,"Taxa_names"=lineage_tips, "MRCA's"=lineage_mrca,"lineage_state"=lineage_state,"Halfedge_over_tmrca"=halfedge_over_mrca))

}


#' LineageHomology_w_uncertainty_v2
#' @description
#'  LineageHomology_w_uncertainty_v2 returns the same outputs as LineageHomology.
#'  This is the same function as LineageHomology_w_uncertainty but faster.
#'  LineageHomology counts transmission lineages according to state transitions between ancestral and descendant nodes that have a probability higher than 50 percent.
#'  LineageHomolog_w_uncertainty instead samples states from the posterior probability of states (that is usually included in a phylogeographical consensus tree) for each node. Thus, the transmission lineage division is probabilistic and retains the uncertainty in the posterior distribution if the function is run multiple times, e.g., using base::replicate().
#'
#' @param tree Phylogenetic tree
#' @param ace_nodes Node probabilities in ace format
#' @param ace_tips Tip probabilities in ace format
#' @param give_tips Constrain the computations to a set of tips, e.g. all the names of the tips in a certain state (e.g. one specific geographical location).
#' @param start_time Date of root in the phylogeny.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' multi_counts = pbreplicate(
#'   100,
#'   LineageHomology::LineageHomology_w_uncertainty_v2(
#'     tree,
#'     ace_nodes = test$lik.anc,
#'     give_tips = Norwegian_tips,
#'     ace_tips = to.matrix(Locations, seq = c("Norway", "RoW")),
#'     start_time = start_date
#'   )
#' )
#' }
#'
LineageHomology_w_uncertainty_v2 = function(tree,ace_nodes,ace_tips, give_tips=NA,start_time=NA) {

  #The node number of the tips are
  tip_nodes = sapply(tree$tip.label,function(x,y) which(y==x),y=tree$tip.label)

  if(is.na(give_tips)==F){ #If the count is among given tips.
    tip_nodes = sapply(give_tips,function(x,y) which(y==x),y=tree$tip.label)
  }
  #The algorithm proceeds by removing tips from the tip nodes list until there is none left.

  tip_states = ace_tips
  node_states = ace_nodes
  #Sample nodes from the ASR once here, so they don't have to be resampled multiple times.
  node_draws = apply(ace_nodes, 1,FUN=function(x) sample(1:2,1, prob=x))
  rownames(node_states)=(length(tree$tip.label)+1):((length(tree$tip.label))+nrow(ace_nodes)) #Old specification was buggy due to missing nodes.
  tip_label_list = tree$tip.label
  #Parent and child.
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]


  #Size of lineage
  lineage_size = c()
  #Size of lineage
  lineage_state = c()
  #Tips in lineage
  lineage_tips = list()
  #MRCA node list
  mrca.list = c()
  #Initialize importations,
  importations = 0

  #Initialize
  local_transmissions = 0
  #Initialize index counter
  counter=1



  checked_most_basal = F

  while(length(tip_nodes)>0) {
    current_node = tip_nodes[1]
    #print(current_node)
    tip_label_list = c(names(current_node)) #Set up a list containing tips for the current tip.
    tip_nodes = tip_nodes[-1] #Remove tip from search space immediately.
    current_state = which(tip_states[current_node,]==1) #If other columns than this has more than 50%, there is a transition.
    lineage_state = c(lineage_state,current_state)  #Save the state of the lineage.
    goes_all_the_way_back = F
    tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    ancestor_node = parent[which(child==current_node)]
    ancestor_node_row = which(rownames(node_states)==ancestor_node)
    ancestor_state =node_draws[ancestor_node_row]

    if(current_state!=ancestor_state) {  #If more than one mapping on the branch, then it is estimated as an importation.
      tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    }

    else {
      while(current_state==ancestor_state) { #While there is no transition, move one node up.

        if(((current_node%in%child)==T)) {
          current_node = ancestor_node #Move one node up.
          ancestor_node = parent[which(child==current_node)]
          if((current_node%in%child)==F) { #break loop if we're already at the root node.
            goes_all_the_way_back = T
            #print("Broke")
            break()
          }
          #Note that current state doesn't need to be updated since we're moving up if it is the same.
          ancestor_node_row = which(rownames(node_states)==ancestor_node)
          ancestor_state =node_draws[ancestor_node_row]
        }

      }
      tip_search_subspace = descendants(tree, current_node)
      subspace=tip_search_subspace[tip_search_subspace%in%tip_nodes] #Find all from the set that can be included.

      ####  Must add one section for the nodepath we just went through sampling up to the ancestor.###
      # Nodepath from current_node to tip node - all the internal will be the same state since we were allowed to move.
      # Add a starting node in the script, then find the nodepath to this one. Setup a list with these nodes, and replace the
      # whichmax nodes - soon to be sampled nodes, with these.

      if(length(subspace)!=0) {
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip this - should become buggy with sampling.
          #nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          nodes_to_check = head(nodepath_quick_w_mrca_nodenumber(tree = tree,taxa = subspace[i],mrca = current_node),-1)
          nods= sapply(nodes_to_check, FUN = function(x) node_draws[which(rownames(node_states)==x)]) #node states
          tps = which.max(tip_states[subspace[i],]) #Check state of tip
          nods = c(nods, tps)

          if(sum(nods==current_state)==length(nods)) { #If all states are equal to the state of the current node proceed.
            tips_in_segment = tips_in_segment+1
            tip_label_list = c(tip_label_list, names(tip_nodes[which(tip_nodes==subspace[i])])) #Add the tips that will be removed.
            tip_nodes = tip_nodes[-which(tip_nodes==subspace[i])] #Remove found tip from search space.
          }
        }
      }
    }

    if(tips_in_segment > 1 ) {
      local_transmissions = local_transmissions+(tips_in_segment-1)
      importations = importations+1
    }
    if(tips_in_segment==1) {
      importations = importations+1
    }
    lineage_size = c(lineage_size,tips_in_segment)

    lineage_tips[[counter]]=tip_label_list
    mrca.list[counter] = current_node
    counter = counter+1 #move one up to fill the next index on in the next loop iteration.
    #debug: print(paste0("imports are:",importations))
    #debug: print(paste0("Local transmissions are:",local_transmissions))
  }
  #process the mrca nodes to ages.
  names(lineage_tips)=paste0("Lineage no: ", 1:length(lineage_tips))
  lineage_mrca = unlist(lapply(mrca.list, FUN = function(x) nodeheight(tree,node = x)))
  if(start_time){
    lineage_mrca =lineage_mrca+start_time
  }

  #For transition dates it might be more appropriate to use the midpoint of the edge ancestral to the lineage
  subtract_ancestral_edge_time = unlist(lapply(mrca.list, FUN = function(x) find_midedge_ancestral(tree,node = x)))
  halfedge_over_mrca =lineage_mrca - subtract_ancestral_edge_time

  return(list("Import_LocalTransmission"=c(importations, local_transmissions),"Lineage_sizes"=lineage_size,"Taxa_names"=lineage_tips, "MRCA's"=lineage_mrca,"lineage_state"=lineage_state,"Halfedge_over_tmrca"=halfedge_over_mrca))

}

