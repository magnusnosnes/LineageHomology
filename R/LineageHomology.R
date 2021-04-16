library(phytools)
library(ape)
library(phangorn)
library(TreeTools)
library(BactDating)
library(ips)


#Algorithm:
#Fint the parent nodes of all tips.
#For each tip; check for transition to its parent node
#If transition, the algorithm has found a singleton.
#If not, move on branch up, and check any new descendant emerging from the parent node.
#If any descendant is found this is no longer a singleton, but represents a local tranmission event.
#If no descendants is found, this is just a singleton with a long branch.

#For testing
#consensus_mapping = consensus_Q

#For testing

#' Title
#'
#' @param tree #Phylogenetic tree
#' @param ace_nodes #Node probabilities in ace format
#' @param ace_tips #Tip probabilities in ace format
#' @param give_tips #Taxa to do the computations for
#' @param start_time Date of root in the phylogeny.
#'
#' @return
#' @export
#'
#' @importFrom phytools nodeHeights nodeheight
#' @importFrom ips descendants
#' @importFrom ape nodepath
#' @importFrom ape nodepath
#' @examples
#' \dontrun{
#' #LineageHomology(tree_test, ace_nodes=fit1$lik.anc,
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
  rownames(node_states)=(length(tree$tip.label)+1):(2*length(tree$tip.label)-1)
  tip_label_list = tree$tip.label
  #Parent and child.
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]


  #Size of lineage
  lineage_size = c()
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
          #tps = which.max(tip_states[nodes_to_check[length(nodes_to_check)],])
          #both = c(nods, tps)

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
    #print(paste0("imports are:",importations))
    #print(paste0("Local transmissions are:",local_transmissions))
  }
  #process the mrca nodes to ages.
  names(lineage_tips)=paste0("Lineage no: ", 1:length(lineage_tips))
  lineage_mrca = unlist(lapply(mrca.list, FUN = function(x) nodeheight(tree,node = x)))
  if(start_time){
    lineage_mrca =lineage_mrca+start_time
  }
  return(list("Import_LocalTrans"=c(importations, local_transmissions),"Lineage_sizes"=lineage_size,"Taxa_names"=lineage_tips, "MRCA's"=lineage_mrca))

}

#' Title
#'
#' @param tree #Phylogenetic tree
#' @param ace_nodes #Node probabilities in ace format
#' @param ace_tips #Tip probabilities in ace format
#' @param give_tips #Taxa to do the computations for
#' @param start_time Date of root in the phylogeny.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' #LineageHomology_w_uncertainty(tree_test, ace_nodes=fit1$lik.anc,
#' ace_tips = to.matrix(loc, seq=c("Norway", "RoW")), start_time=2000)
#' }

LineageHomology_w_uncertainty = function(tree,ace_nodes,ace_tips, give_tips=NA,start_time=NA) {

  #NB. There is some bug in here. Must be fixed before the Uk report can be updated.
  #Line for debugging using the UK script
  #tree = tree; ace_nodes = castor$ancestral_likelihoods; give_tips = Norwegian_tips; ace_tips = to.matrix(Locations,seq=c("Norway","RoW"));start_time = start_date

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
  rownames(node_states)=(length(tree$tip.label)+1):(2*length(tree$tip.label)-1)
  tip_label_list = tree$tip.label
  #Parent and child.
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]


  #Size of lineage
  lineage_size = c()
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
          nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          internal = nodes_to_check[-length(nodes_to_check)]
          nods= sapply(internal, FUN = function(x) node_draws[which(rownames(node_states)==x)]) #node states
          ##The following are unnecessary as tip always have the same state in this example
          #tps = which.max(tip_states[nodes_to_check[length(nodes_to_check)],])
          #both = c(nods, tps)

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
    #print(paste0("imports are:",importations))
    #print(paste0("Local transmissions are:",local_transmissions))
  }
  #process the mrca nodes to ages.
  names(lineage_tips)=paste0("Lineage no: ", 1:length(lineage_tips))
  lineage_mrca = unlist(lapply(mrca.list, FUN = function(x) nodeheight(tree,node = x)))
  if(start_time){
    lineage_mrca =lineage_mrca+start_time
  }
  return(list("Import_LocalTransmission"=c(importations, local_transmissions),"Lineage_sizes"=lineage_size,"Taxa_names"=lineage_tips, "MRCA's"=lineage_mrca))

}



