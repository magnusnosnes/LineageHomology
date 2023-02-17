library(ape)
library(lubridate)

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

  if(length(give_tips)>1){ #If the count is among given tips.
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


#toolfunction for the function below.



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

  if(length(give_tips)>1){ #If the count is among given tips.
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

      #current node is now the most ancestral node in the TL
      #Must find all descendent taxa that were not included in the path
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
      #Add nodes that are in the TL but do not lead to a observed sequence.


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


#   ____________________________________________________________________________
#   Utility functions for LineageHomology_w_uncertainty_v2                  ####

#After the most ancestral node has been found, check for descendant nodes that
#do not lead to sequences. (These are technically also in the TL and signify local transmission
#in that location)
find_remaining_nodes_tl = function(found_nodes, current_state, node_draws, parent, child) {

  # Pesudocode
  #While possible children is not zero search for more.
  while(1==1) {

    # Find children not in found nodes
    possible_children_index = which(parent %in% found_nodes)
    #Check that those children is not in the found node set.
    possible_children=child[possible_children_index][((child[possible_children_index] %in% found_nodes) == F)]
    # Check that the (probability drawn) state of those children matches the current state.
    possible_children_state = node_draws[names(node_draws)%in%possible_children]

    new_nodes = as.integer(names(possible_children_state[which(possible_children_state == current_state)]))

    #Break the while loop if there are no new nodes to add.
    if((length(new_nodes)>0)==FALSE) {
      break()
    }
    found_nodes = c(found_nodes, new_nodes)
  }
  found_nodes
}

# find_observed_TLs defines TLs by finding connected segments of tip nodes,
# until the number of tip_nodes remaining is zero.
find_observed_TLs = function(tip_nodes, tip_states, node_draws, parent, child){

  remaining_nodes = tip_nodes
  current_node_set = tip_nodes[1]
  mrca_node = current_node_set

  #Size of lineage
  lineage_size = c()
  #Size of lineage
  lineage_state = c()
  #Tips in lineage
  lineage_tips = list()
  #Tree nodes in lineage.
  lineage_nodes=list() #Including tip nodes
  #MRCA node list
  mrca_list = c()
  #Initialize importations,
  importations = 0
  #Initialize local transmissions
  local_transmissions = 0
  #Initialize index counter
  counter=1


  #Pseudocode:
  #1. Define the current state, a nodeset, mrca_node, remaining_nodes, current_node_set, counter for indexes to fill.
  #2. while lenght(remaining_nodes)
  #2. Check for parents and children not in the nodeset.
  #3. Grow the nodeset if any of the possible new nodes matches the current state.
  #4. Update mrca node if a new parent node is found.
  #5. Check if nodeset includes any tip nodes, if so remove them from the remaining nodes.
  #6.

  while(length(remaining_nodes)>0) {
  counter = counter+1
  }

    return(list("tl_size"=lineage_size,
              "tl_states"=lineage_state,
              "tl_tips" = lineage_tips,
              "tl_nodes"=lineage_nodes,
              "tl_mrcas"=mrca_list,
              "importations"=1,
              "local_transmissions"=1))


#   ____________________________________________________________________________
#   Old code to be replaced by this function                                ####

  while(length(tip_nodes)>0) {
    current_node = tip_nodes[1]
    #print(current_node)
    tip_label_list = c(names(current_node)) #Set up a list containing taxa names for the lineage.
    node_label_list = c(current_node)  #Set up a list node names for the lineage.
    tip_nodes = tip_nodes[-1] #Remove tip from search space immediately.
    current_state = which(tip_states[current_node,]==1) #If other columns than this has more than 50%, there is a transition.
    lineage_state = c(lineage_state,current_state)  #Save the state of the lineage.
    tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    ancestor_node = parent[which(child==current_node)]
    ancestor_node_row = which(rownames(node_states)==ancestor_node)
    ancestor_state =node_draws[ancestor_node_row]

    if(current_state!=ancestor_state) {  #If more than one mapping on the branch, then it is estimated as an importation.
      tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    }

    else {
      while(current_state==ancestor_state) { #While there is no state transition, move one node up.

        if(((current_node%in%child)==T)) {
          current_node = ancestor_node #Move one node up.
          node_label_list =c(node_label_list, current_node) #Save new node number to node list.
          ancestor_node = parent[which(child==current_node)]

          if((current_node%in%child)==F) { #break loop if we're already at the root node.
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

      if(length(subspace)!=0) { #Check if there are any descendant tips that haven´t been removed from the search space yet.
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip this
          #nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          nodes_to_check = head(nodepath_quick_w_mrca_nodenumber(tree = tree,taxa = subspace[i],mrca = current_node),-1) #head, -1 removes the last entry
          nods= sapply(nodes_to_check, FUN = function(x) node_draws[which(rownames(node_states)==x)]) #node states
          tps = which.max(tip_states[subspace[i],]) #Check state of tip
          nods = c(nods, tps)


          if(sum(nods==current_state)==length(nods)) { #If all states are equal to the state of the current node proceed.
            tips_in_segment = tips_in_segment+1
            tip_label_list = c(tip_label_list, names(tip_nodes[which(tip_nodes==subspace[i])])) #Add the tips that will be removed.
            node_path= c(nodes_to_check,subspace[i]) #Node path to tips including the tip node.
            nodes_not_added_indexes = which((node_path %in% node_label_list)==F) #Which of the nodes are not already in the node list.
            nodes_to_add = node_path[nodes_not_added_indexes] #Select the corresponding nodes.
            node_label_list =c(node_label_list,nodes_to_add) #Add them.
            tip_nodes = tip_nodes[-which(tip_nodes==subspace[i])] #Remove found tip from search space.
          }
        }
      }
      #Add nodes that are in the TL but do not lead to a observed sequence.
      node_label_list = find_remaining_nodes_tl(found_nodes = node_label_list,
                                                current_state = current_state,
                                                node_draws =  node_draws,
                                                parent = parent,
                                                child = child) #Function defined above this one.
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
    lineage_nodes[[counter]] = node_label_list
    mrca_list[counter] = current_node

    #move one up to fill the next index on in the next loop iteration.
    counter = counter+1

    #debug: print(paste0("imports are:",importations))
    #debug: print(paste0("Local transmissions are:",local_transmissions))

  }

}

#Utility function that finds nodes unassigned to TLs and returns the nodes in addition to their mapped state.
find_unassigned_nodes = function(tree, lineage_nodes) {
  all_nodes = unique(as.vector(tree$edge))
  assigned_nodes = unlist(lineage_nodes)
  #Any tip will by definition already be assigned.
  #We only have to search for internal nodes.
  unassigned_nodes = all_nodes[(all_nodes %in% assigned_nodes ) == F]
  unassigned_nodes
}

#Utility function that finds TLs form unassigned nodes.
#Must add conditions for these not to execute.
find_unobserved_TLs = function(unassigned_nodes, node_draws,name_of_states, parent, child){ #Parent and child are tree$edge[,1] and tree$edge[,2] respectively.

  #   ____________________________________________________________________________
  #   Debug                                                                   ####
  #unassigned_nodes = unassigned_nodes;  node_draws = node_draws;   name_of_states = colnames(node_states);   parent = parent; child=child;

  #   ____________________________________________________________________________
  #   End debug                                                               ####
  #initialize unobserved lineages.
  unobserved_lineages = list()
  #State of the unobserved lineages
  unobserved_lineages_state = c()
  # Set up a list with nodes to remove
  remaining_nodes = unassigned_nodes
  # While the list is not empty
  counter = 1

  while(length(remaining_nodes)>0) {


    #   ____________________________________________________________________________
    #   Debug                                                                   ####

    #if(counter==11) break

    #   ____________________________________________________________________________
    #   end debug                                                               ####

    #Set up a node set and expand this by parents and children satisfying criteria each time.
    current_node_set = remaining_nodes[1]
    #Remove this node from the remaining nodes
    remaining_nodes = remaining_nodes[-1]
    current_node_set_state = node_draws[names(node_draws)==current_node_set]
    unobserved_lineages_state=c(unobserved_lineages_state,current_node_set_state) #State is determined on the line above. Store it here.
    #While possible children and possible parents is not zero.

    condition_to_continue = T



    while(condition_to_continue == T) {

      #.Search for possible children nodes.
      possible_children = child[which(parent %in% current_node_set)]
      possible_children = possible_children[possible_children %in% remaining_nodes]
      #Check state of the possible_children
      correct_state_child_ind = node_draws[as.character(possible_children)]==current_node_set_state #Index the children states by
      correct_child = possible_children[correct_state_child_ind]

      #Search for possible parent nodes in the same fashion as for the children.
      # Check that they all have the same state.
      possible_parents = parent[which(child %in% current_node_set)]
      possible_parents = possible_parents[possible_parents %in% remaining_nodes]
      #Check state of the possible_parents
      correct_state_parent_ind = node_draws[as.character(possible_parents)]==current_node_set_state #Index the children states by
      correct_parent = possible_parents[correct_state_parent_ind]

      condition_to_continue=(length(c(correct_parent, correct_child))>0) #Update loop condition.

      if(condition_to_continue) {
        #Current_node_set
        current_node_set = c(current_node_set,correct_parent,correct_child)
        #Update remaining nodes
        remaining_nodes = remaining_nodes[-c((which(remaining_nodes %in% current_node_set)))]
      }
      }

    # Save TL
    unobserved_lineages[[counter]]=current_node_set
    # remove TL from remaining nodes.

    #Update counter
    counter = counter + 1

    # Return

  }
  #if condition needed because we cannot set attributes of null´s.

  if(length(unobserved_lineages_state)>0){
    names(unobserved_lineages_state) = name_of_states[unobserved_lineages_state] #Name the states with the corresponding state.
  }

  return(list("unobserved_tls"=unobserved_lineages,"unobserved_tl_state"=unobserved_lineages_state))

}

#Utility function that assigns dates to the unassigned nodes.
assign_date_unobserved_TLs = function(tree, unobserved_tl, root_time) {



  root_node = unique(as.vector(tree$edge))[(unique(as.vector(tree$edge)) %in% tree$edge[,2])==F]

  if(length(unobserved_tl$unobserved_tls)>0) {

    halfedge_over_mrca = c()
    for(i in 1:length(unobserved_tl$unobserved_tls)) {


      nodes_to_consider = unobserved_tl$unobserved_tls[i][[1]]

      if(root_node %in% nodes_to_consider) {
        new_height = root_time
        halfedge_over_mrca = c(halfedge_over_mrca,new_height)
      }
      else {
        parental_edges = which(tree$edge[,2] %in% nodes_to_consider) #indexes of parental branches
        nodes_with_ancestral_branch = tree$edge[,2][parental_edges] #Nodes with parental branches
        heights = unlist(lapply(nodes_with_ancestral_branch, FUN = function(x) nodeheight(tree,node = x))) #Times of nodes with parental branches
        heights = heights - (tree$edge.length[parental_edges]/2) #Subtract half of the ancestral branch
        new_height = heights[which.min(heights)]+root_time #Add root time
        halfedge_over_mrca = c(halfedge_over_mrca,new_height) #Add the
      }

    }
    return(halfedge_over_mrca)

  }

  else {
    return(NULL)
  }


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


  #   ____________________________________________________________________________
  #   Debug                                                                   ####

      #give_tips=NA;     tree=tree_test;    ace_nodes=fit1$lik.anc;    ace_tips = to.matrix(loc, seq=c("Norway", "RoW"));    start_time = 2000;

    #For importation runs:
    #tree=tree;   ace_nodes=ace_NOR$lik.anc;   ace_tips = to.matrix(locations_NOR, seq=c("Norway", "ROW"));   start_time = 1750; give_tips = NA

  #   ____________________________________________________________________________
  #   End debug                                                               ####

  #The node number of the tips are
  tip_nodes = sapply(tree$tip.label,function(x,y) which(y==x),y=tree$tip.label)

  if(length(give_tips)>1){ #If the count is among given tips.
    tip_nodes = sapply(give_tips,function(x,y) which(y==x),y=tree$tip.label)
  }
  #The algorithm proceeds by removing tips from the tip nodes list until there is none left.

  tip_states = ace_tips
  node_states = ace_nodes
  #Sample nodes from the ASR once here, so they don't have to be resampled multiple times.
  node_draws = apply(ace_nodes, 1,FUN=function(x) sample(1:2,1, prob=x))
  #Set rownames for the matrix
  rownames(node_states)=(length(tree$tip.label)+1):((length(tree$tip.label))+nrow(ace_nodes))
  #Set names for the vector of samples states.
  names(node_draws)=(length(tree$tip.label)+1):((length(tree$tip.label))+nrow(ace_nodes))
  #Set of a list for the names of the sequences.
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
  #Tree nodes in lineage.
  lineage_nodes=list() #Including tip nodes
  #MRCA node list
  mrca_list = c()
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
    tip_label_list = c(names(current_node)) #Set up a list containing taxa names for the lineage.
    node_label_list = c(current_node)  #Set up a list node names for the lineage.
    tip_nodes = tip_nodes[-1] #Remove tip from search space immediately.
    current_state = which(tip_states[current_node,]==1) #If other columns than this has more than 50%, there is a transition.
    lineage_state = c(lineage_state,current_state)  #Save the state of the lineage.
    tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    ancestor_node = parent[which(child==current_node)]
    ancestor_node_row = which(rownames(node_states)==ancestor_node)
    ancestor_state =node_draws[ancestor_node_row]

    if(current_state!=ancestor_state) {  #If more than one mapping on the branch, then it is estimated as an importation.
      tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    }

    else {
      while(current_state==ancestor_state) { #While there is no state transition, move one node up.

        if(((current_node%in%child)==T)) {
          current_node = ancestor_node #Move one node up.
          node_label_list =c(node_label_list, current_node) #Save new node number to node list.
          ancestor_node = parent[which(child==current_node)]

          if((current_node%in%child)==F) { #break loop if we're already at the root node.
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

      if(length(subspace)!=0) { #Check if there are any descendant tips that haven´t been removed from the search space yet.
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip this
          #nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          nodes_to_check = head(nodepath_quick_w_mrca_nodenumber(tree = tree,taxa = subspace[i],mrca = current_node),-1) #head, -1 removes the last entry
          nods= sapply(nodes_to_check, FUN = function(x) node_draws[which(rownames(node_states)==x)]) #node states
          tps = which.max(tip_states[subspace[i],]) #Check state of tip
          nods = c(nods, tps)


          if(sum(nods==current_state)==length(nods)) { #If all states are equal to the state of the current node proceed.
            tips_in_segment = tips_in_segment+1
            tip_label_list = c(tip_label_list, names(tip_nodes[which(tip_nodes==subspace[i])])) #Add the tips that will be removed.
            node_path= c(nodes_to_check,subspace[i]) #Node path to tips including the tip node.
            nodes_not_added_indexes = which((node_path %in% node_label_list)==F) #Which of the nodes are not already in the node list.
            nodes_to_add = node_path[nodes_not_added_indexes] #Select the corresponding nodes.
            node_label_list =c(node_label_list,nodes_to_add) #Add them.
            tip_nodes = tip_nodes[-which(tip_nodes==subspace[i])] #Remove found tip from search space.
          }
        }
      }
      #Add nodes that are in the TL but do not lead to a observed sequence.
      node_label_list = find_remaining_nodes_tl(found_nodes = node_label_list,
                                                current_state = current_state,
                                                node_draws =  node_draws,
                                                parent = parent,
                                                child = child) #Function defined above this one.
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
    lineage_nodes[[counter]] = node_label_list
    mrca_list[counter] = current_node
    counter = counter+1 #move one up to fill the next index on in the next loop iteration.
    #debug: print(paste0("imports are:",importations))
    #debug: print(paste0("Local transmissions are:",local_transmissions))
    #

  }


  #1. Find unobserved nodes
  unassigned_nodes = find_unassigned_nodes(tree,
                                           lineage_nodes=lineage_nodes)
  #unassigned_nodes; length(unassigned_nodes);
  #2. Find connected state groups of parent and children nodes among the unobserved nodes.

  unobserved_TLs = find_unobserved_TLs(unassigned_nodes = unassigned_nodes,
                                       node_draws = node_draws,
                                       name_of_states = colnames(node_states),
                                       parent = parent, child=child)
#   ____________________________________________________________________________
#   Debug                                                                   ####

  #print(unobserved_TLs)
  #print(paste0("Reason for bug:",unobserved_TLs$unobserved_tls))
  #Assign import dates to these
  # print(paste0("Hello", unobserved_TLs$unobserved_tls[[1]]))

#   ____________________________________________________________________________
#   End debug                                                                   ####

  unobserved_TLs_halfedge_above_mrca = assign_date_unobserved_TLs(tree=tree,
                                                                  unobserved_tl = unobserved_TLs,
                                                                  root_time = start_time)
  #Process the mrca nodes to ages.
  names(lineage_tips)=paste0("Lineage no: ", 1:length(lineage_tips))
  lineage_mrca = unlist(lapply(mrca_list, FUN = function(x) nodeheight(tree,node = x)))
  if(start_time){
    lineage_mrca =lineage_mrca+start_time
  }


  #For transition dates it might be more appropriate to use the midpoint of the edge ancestral to the lineage
  subtract_ancestral_edge_time = unlist(lapply(mrca_list, FUN = function(x) find_midedge_ancestral(tree,node = x)))
  halfedge_over_mrca = lineage_mrca - subtract_ancestral_edge_time

#   ____________________________________________________________________________
#   Debug                                                                   ####

  #unassigned_nodes ;  unobserved_TLs;   unobserved_TLs_halfedge_above_mrca

#   ____________________________________________________________________________
#   End Debug                                                               ####

   return(list("Import_LocalTransmission" = c(importations, local_transmissions),
              "Lineage_sizes" = lineage_size,
              "Taxa_names" = lineage_tips,
              "MRCA's" = lineage_mrca,
              "lineage_state" = lineage_state,
              "Halfedge_over_tmrca" = halfedge_over_mrca,
              "lineage_nodes" = lineage_nodes,
              "unobserved_tl_nodes"=unobserved_TLs$unobserved_tls,
              "unobserved_tl_states"=unobserved_TLs$unobserved_tl_state,
              "unobserved_tl_halfedge_above_mrca"=unobserved_TLs_halfedge_above_mrca)
         )

}




#   ____________________________________________________________________________
#   Utility function that calculates exports and export times from the nodes in a TL####

find_exports_from_nodes = function(tree,tl_nodes, tip_nodes,start_time) {

  ## DEBUG
  #tree = tree; tip_nodes = unchanged_tip_nodes; start_time = start_time; tl_nodes = lineage_nodes;
  ## END DEBUG

  #For each transmission lineage.
  tl_exports = c()
  tl_export_times = list()

  for(i in 1:length(tl_nodes)) {
    group_nodes = tl_nodes[[i]]
    group_tip_nodes = group_nodes[group_nodes %in% tip_nodes]
    condition1 = (tree$edge[,1] %in% group_nodes) #Parent node is in TL
    condition2 = (tree$edge[,2] %in% c(group_nodes,group_tip_nodes))==F #Child node is not in TL group nodes, Child node is not among TL taxa.
    exportation_edges = which(condition1 & condition2) # Potential bug here with the array notation

    if(length(exportation_edges)>0) { #if the TL produced exports
      ancestral_edge_exportation = tree$edge[exportation_edges,1] #Find parent edge
      time_ancestral_edge = start_time + unlist(lapply(ancestral_edge_exportation, FUN = function(x) nodeheight(tree, x))) #This is slow, could speed everything up by replacing
      half_edge_to_child = tree$edge.length[exportation_edges]/2 #Midpoint of leading to the exportation
      export_times = time_ancestral_edge+half_edge_to_child
      n_exports = length(export_times)
      tl_exports = c(tl_exports, n_exports)
      tl_export_times[[i]] = export_times
    }
    else {
      tl_exports = c(tl_exports, 0)
      tl_export_times[[i]] = NA
    }

  }
  exports = sum(tl_exports)
  return(list("Exports"=exports, "Exports_in_TL"=tl_exports, "Export_times_from_TL"=tl_export_times))
}



#' LineageHomology_v2

#' @description
#'  LineageHomology_v2 is faster version of LineageHomology.
#'  LineageHomology counts transmission lineages according to state transitions between ancestral and descendant nodes that have a probability higher than 50 percent.
#'  LineageHomolog_w_uncertainty instead samples states from the posterior probability of states (that is usually included in a phylogeographical consensus tree) for each node. Thus, the transmission lineage division is probabilistic and retains the uncertainty in the posterior distribution if the function is run multiple times, e.g., using base::replicate().
#'
#' @param tree Phylogenetic tree
#' @param ace_nodes Node probabilities in ace format
#' @param ace_tips Tip probabilities in ace format
#' @param give_tips Constrain the computations to a set of tips, e.g. all the names of the tips in a certain state (e.g. one specific geographical location).
#' @param start_time Date of root in the phylogeny.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Lineage_sizes}{A vector containing the number of tips that belong to each transmission lineage.}
#'   \item{Import_Export_LocalTransmission}{A list containing the estimated number of importation, local transmission and export events.}
#'   \item{Import_times}{A vector containing the times at the middle of the ancestral edge of the TLs.}
#'   \item{Exports_from_TL}{A vector containing the number of exports from each transmission lineage.}
#'   \item{Exports_times_from_TL}{A list of vectors containing the times of exports from each transmission lineage.}
#'   \item{Genome_identifiers_in_TLs}{A list of lists where each sub list contains the names of the tips in each transmission lineage. The order is the same as the listed order of the Lineage_sizes.}
#'   \item{TMRCA}{A vector containing the time of the most recent common ancestor of each transmission lineage.}
#'   \item{lineage_state}{A vector containing the state each transmission lineage belongs to. The states are treated as numbers in the analyses, so the relative ordering has to be looked up.}
#'   \item{lineage_nodes}{A list of lists where each sub list contains the node IDs of the nodes that belong to each transmission lineage. The order is the same as the listed order of the Lineage_sizes.}
#'   \item{unobserved_tl_nodes}{A list containing the node IDs of the unobserved transmission lineages.}
#'   \item{unobserved_tl_states}{A vector containing the state of each unobserved transmission lineage.}
#'   \item{unobserved_tl_halfedge_above_mrca}{A vector containing the time points on the middle of the edges ancestral to each unobserved transmission lineage's MRCA. Used in analyses of importation and local transmission.}
#' }
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
LineageHomology_v2 = function(tree,ace_nodes,ace_tips, give_tips=NA,start_time=NA) {


  # Pseudocode.
  # The algorithm proceeds by removing tips from the tip nodes list until there is none left.

  #   ____________________________________________________________________________
  #   Debug                                                                   ####

  #give_tips=NA;     tree=tree_test;    ace_nodes=fit1$lik.anc;    ace_tips = to.matrix(loc, seq=c("Norway", "RoW"));    start_time = 2000;

  #For importation runs:
  #tree=tree;   ace_nodes=ace_NOR$lik.anc;   ace_tips = to.matrix(locations_NOR, seq=c("Norway", "ROW"));   start_time = 1750; give_tips = NA

  #   ____________________________________________________________________________
  #   End debug                                                               ####

  #The node number of the tips are
  tip_nodes = sapply(tree$tip.label,function(x,y) which(y==x),y=tree$tip.label)
  unchanged_tip_nodes = tip_nodes

  if(length(give_tips)>1){ #If the count is among given tips.
    tip_nodes = sapply(give_tips,function(x,y) which(y==x),y=tree$tip.label)
  }

  tip_states = ace_tips
  node_states = ace_nodes
  #Sample nodes from the ASR once here, so they don't have to be resampled multiple times.
  node_draws = apply(ace_nodes, 1,FUN=function(x) which.max(x))
  #Set rownames for the matrix
  rownames(node_states)=(length(tree$tip.label)+1):((length(tree$tip.label))+nrow(ace_nodes))
  #Set names for the vector of samples states.
  names(node_draws)=(length(tree$tip.label)+1):((length(tree$tip.label))+nrow(ace_nodes))
  #Set of a list for the names of the sequences.
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
  #Tree nodes in lineage.
  lineage_nodes=list() #Including tip nodes
  #MRCA node list
  mrca_list = c()


  #Initialize importations,
  importations = 0
  #Initialize
  local_transmissions = 0
  #Initialize index counter
  exports = 0


  # Lineage counter used to fill lists.
  counter=1

  checked_most_basal = F

  while(length(tip_nodes)>0) {
    current_node = tip_nodes[1]
    #print(current_node)
    tip_label_list = c(names(current_node)) #Set up a list containing taxa names for the lineage.
    node_label_list = c(current_node)  #Set up a list node names for the lineage.
    tip_nodes = tip_nodes[-1] #Remove tip from search space immediately.
    current_state = which(tip_states[current_node,]==1) #If other columns than this has more than 50%, there is a transition.
    lineage_state = c(lineage_state,current_state)  #Save the state of the lineage.
    tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    ancestor_node = parent[which(child==current_node)]
    ancestor_node_row = which(rownames(node_states)==ancestor_node)
    ancestor_state =node_draws[ancestor_node_row]

    if(current_state!=ancestor_state) {  #If more than one mapping on the branch, then it is estimated as an importation.
      tips_in_segment = 1 #Count how many tips are considered for each mapping segment.
    }

    else {
      while(current_state==ancestor_state) { #While there is no state transition, move one node up.

        if(((current_node%in%child)==T)) {
          current_node = ancestor_node #Move one node up.
          node_label_list =c(node_label_list, current_node) #Save new node number to node list.
          ancestor_node = parent[which(child==current_node)]

          if((current_node%in%child)==F) { #break loop if we're already at the root node.
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


      if(length(subspace)!=0) { #Check if there are any descendant tips that haven´t been removed from the search space yet.
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip this
          #nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          nodes_to_check = head(nodepath_quick_w_mrca_nodenumber(tree = tree,taxa = subspace[i],mrca = current_node),-1) #head, -1 removes the last entry
          nods= sapply(nodes_to_check, FUN = function(x) node_draws[which(rownames(node_states)==x)]) #node states
          tps = which.max(tip_states[subspace[i],]) #Check state of tip
          nods = c(nods, tps)


          if(sum(nods==current_state)==length(nods)) { #If all states are equal to the state of the current node proceed.
            tips_in_segment = tips_in_segment+1
            tip_label_list = c(tip_label_list, names(tip_nodes[which(tip_nodes==subspace[i])])) #Add the tips that will be removed.
            node_path= c(nodes_to_check,subspace[i]) #Node path to tips including the tip node.
            nodes_not_added_indexes = which((node_path %in% node_label_list)==F) #Which of the nodes are not already in the node list.
            nodes_to_add = node_path[nodes_not_added_indexes] #Select the corresponding nodes.
            node_label_list =c(node_label_list,nodes_to_add) #Add them.
            tip_nodes = tip_nodes[-which(tip_nodes==subspace[i])] #Remove found tip from search space.
          }
        }
      }
      #Add nodes that are in the TL but do not lead to a observed sequence.
      node_label_list = find_remaining_nodes_tl(found_nodes = node_label_list,
                                                current_state = current_state,
                                                node_draws =  node_draws,
                                                parent = parent,
                                                child = child) #Function defined above this one.
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
    lineage_nodes[[counter]] = node_label_list
    mrca_list[counter] = current_node
    counter = counter+1 #move one up to fill the next index on in the next loop iteration.
    #debug: print(paste0("imports are:",importations))
    #debug: print(paste0("Local transmissions are:",local_transmissions))
    #

  }

  #Calculate Exports
  export_list = find_exports_from_nodes(tree = tree,tl_nodes = lineage_nodes,tip_nodes = unchanged_tip_nodes,start_time = start_time)


  #1. Find unobserved nodes
  unassigned_nodes = find_unassigned_nodes(tree,
                                           lineage_nodes=lineage_nodes)

  #2. Find connected state groups of parent and children nodes among the unobserved nodes.
  unobserved_TLs = find_unobserved_TLs(unassigned_nodes = unassigned_nodes,
                                       node_draws = node_draws,
                                       name_of_states = colnames(node_states),
                                       parent = parent, child=child)

  #3. Assign date at the midpoint of the ancestral edge.
  unobserved_TLs_halfedge_above_mrca = assign_date_unobserved_TLs(tree=tree,
                                                                  unobserved_tl = unobserved_TLs,
                                                                  root_time = start_time)
  #Process the mrca nodes to ages.
  names(lineage_tips)=paste0("Lineage no: ", 1:length(lineage_tips))
  lineage_mrca = unlist(lapply(mrca_list, FUN = function(x) nodeheight(tree,node = x)))
  if(start_time){
    lineage_mrca =lineage_mrca+start_time
  }

  #For transition dates it might be more appropriate to use the midpoint of the edge ancestral to the lineage
  subtract_ancestral_edge_time = unlist(lapply(mrca_list, FUN = function(x) find_midedge_ancestral(tree,node = x)))
  halfedge_over_mrca = lineage_mrca - subtract_ancestral_edge_time


  return(list("Lineage_sizes" = lineage_size,
              "Import_Export_LocalTransmission" = c(importations, local_transmissions, export_list$Exports),
              "Import_times" = halfedge_over_mrca,
              "Exports_from_TL" = export_list$Exports_in_TL,
              "Exports_times_from_TL" = export_list$Export_times_from_TL,
              "Genome_identifiers_in_TLs" = lineage_tips,
              "TMRCA" = lineage_mrca,
              "lineage_state" = lineage_state,
              "lineage_nodes" = lineage_nodes,
              "unobserved_tl_nodes"=unobserved_TLs$unobserved_tls,
              "unobserved_tl_states"=unobserved_TLs$unobserved_tl_state,
              "unobserved_tl_halfedge_above_mrca"=unobserved_TLs_halfedge_above_mrca)
  )

}



