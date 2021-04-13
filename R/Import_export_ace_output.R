library(phytools)
library(ape)
library(phangorn)
library(TreeTools)
library(BactDating)
library(ips)

#Simulate
# set.seed(350)
# tree_test = simdatedtree(nsam=10, dateroot=2000)
# tree_test = ladderize(tree_test)
# Q=matrix(c(0.5,0.5,0.5,0.5), nrow=2,ncol=2, byrow=F)
# loc = c("Norway", "Norway","Norway","RoW", "RoW", "Norway", "Norway", "RoW", "RoW", "RoW")
# names(loc) = tree_test$tip.label
# mapped_Q = sim.history(tree_test, Q,nsim=10)
# colnames(Q)=c("Norway","RoW")
# Q
# mapped_Q = make.simmap(tree_test, x = loc, nsim=10, mod="ER")
# plot(mapped_Q[[2]])
# consensus_Q = summary(mapped_Q)
# plot(consensus_Q)
# nodelabels()
# 
# fit1 = ace(x=loc, phy= tree_test, type="discrete", mod="ARD")
# tree_test$tip.label
# fit1
# plot.phylo(tree_test)
# nodelabels(pie=fit1$lik.anc)
# tips = to.matrix(loc,seq=c("Norway", "RoW"))
# tiplabels(pie=tips, cex=0.1)
# fit1$lik.anc

#Algorithm:
#Fint the parent nodes of all tips.
#For each tip, check for transition to its parent node
#If transition, the algorithm has found a singleton.
#If not, move on branch up, and check any new descendant emerging from the parent node.
#If any descendant is found this is no longer a singleton, but represents a local tranmission event.
#If no descendants is found, this is just a singleton with a long branch.

#For testing
#consensus_mapping = consensus_Q

#For testing

count_import_export_ace = function(tree,ace_nodes,ace_tips, give_tips=NA,start_time=NA) {
  
  #The node number of the tips are
  tip_nodes = sapply(tree$tip.label,function(x,y) which(y==x),y=tree$tip.label)
  
  if(is.na(give_tips)==F){ #If the count is among given tips.
    tip_nodes = sapply(give_tips,function(x,y) which(y==x),y=tree$tip.label)
  }
  
  tip_states = ace_tips
  node_states = ace_nodes
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
    ancestor_state =which.max(node_states[ancestor_node_row,])
    
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
          ancestor_state = which.max(node_states[ancestor_node_row,])       
        }
        
      }
      tip_search_subspace = descendants(tree, current_node)
      subspace=tip_search_subspace[tip_search_subspace%in%tip_nodes] #Find all possible tips descending from the current node.
      
      if(length(subspace)!=0) { 
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip.
          nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          internal = nodes_to_check[-length(nodes_to_check)]
          nods= sapply(internal, FUN = function(x) which.max(node_states[which(rownames(node_states)==x),])) #node states
          tps = which.max(tip_states[nodes_to_check[length(nodes_to_check)],]) #Tip state
          both = c(nods, tps)
          
          if(sum(both==current_state)==length(both)) { #If all states are equal to the state of the current node proceed.
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
  return(list("Import_Export"=c(importations, local_transmissions),"Lineage_sizes"=lineage_size,"Taxa_names"=lineage_tips, "MRCA's"=lineage_mrca))  
  
}

count_import_export_ace_uncertainty = function(tree,ace_nodes,ace_tips, give_tips=NA,start_time=NA) {
  
  #The node number of the tips are
  tip_nodes = sapply(tree$tip.label,function(x,y) which(y==x),y=tree$tip.label)
  
  if(is.na(give_tips)==F){ #If the count is among given tips.
    tip_nodes = sapply(give_tips,function(x,y) which(y==x),y=tree$tip.label)
  }
  
  tip_states = ace_tips
  node_states = ace_nodes
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
    ancestor_state =sample(1:2,1,prob=node_states[ancestor_node_row,])
    
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
          ancestor_state = sample(1:2,1,prob=node_states[ancestor_node_row,])       
        }
        
      }
      tip_search_subspace = descendants(tree, current_node)
      subspace=tip_search_subspace[tip_search_subspace%in%tip_nodes] #Find all possible tips descending from the current node.
      
      if(length(subspace)!=0) { 
        for(i in 1:length(subspace)) { #Check for transitions on nodepath to each tip.
          nodes_to_check = nodepath(tree, from = current_node, to=subspace[i])[-1] #path to the tips
          internal = nodes_to_check[-length(nodes_to_check)]
          nods= sapply(internal, FUN = function(x) which.max(node_states[which(rownames(node_states)==x),])) #node states
          tps = which.max(tip_states[nodes_to_check[length(nodes_to_check)],]) #Tip state
          both = c(nods, tps)
          
          if(sum(both==current_state)==length(both)) { #If all states are equal to the state of the current node proceed.
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
  return(list("Import_Export"=c(importations, local_transmissions),"Lineage_sizes"=lineage_size,"Taxa_names"=lineage_tips, "MRCA's"=lineage_mrca))  
  
}
#count_import_export_ace(tree_test, ace_nodes=fit1$lik.anc,ace_tips = to.matrix(loc, seq=c("Norway", "RoW")), start_time=2000)
#Test
# Test = count_import_export_consensus(consensus_Q)
# Test[1]; Test[2]
# 
# #Test on full tree
# q2 = summary(mapped_cladeA)
# 
# #Find out why it doesn't work on the phylogeny - bug in the loop when it goes all the way back. Debug that aspect.
# #DEBUG SECTION
# q1 = summary(mapped_cladeA)
# count_import_export_consensus(q1)[2]
# #END OF DEBUG SECTION
# 
# 
# q1 = summary(mapped_cladeA)
# q2 = summary(mapped_cladeB)
# q3 = summary(mapped_cladeC)
# q4 = summary(mapped_cladeD)
# a_lineages =count_import_export_consensus(q1, give_tips=tipsa)[[2]]
# b_lineages =count_import_export_consensus(q2, give_tips=tipsb)[[2]]
# c_lineages =count_import_export_consensus(q3, give_tips=tipsc)[[2]]
# d_lineages =count_import_export_consensus(q4, give_tips=tipsd)[[2]]
# all_lineages = c(a_lineages,b_lineages,c_lineages,d_lineages)
# 
# paste0("G:",1:length(all_lineages)," S: ", all_lineages)
# atreemap = data.frame(group=paste0("G:",1:length(all_lineages)," S: ", all_lineages),Value=all_lineages)
# pdf("~/Dropbox/10000_Lovers/Junk_used_for_exploratory/Lineage_sizes_norway_consensus_Ngon.pdf",width=10,height=10)
# treemap(atreemap,index=("group"),vSize="Value", type="index", palette="Paired")
# dev.off()
