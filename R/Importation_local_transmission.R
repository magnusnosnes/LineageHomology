#' import_local_transmission
#' @description
#' import_local_transmission takes the output of LineageHomology_w_uncertainty_v2 and produces estimates of
#' the number of importation and local transmission events.
#' @param tree
#' @param LineageHomology_replicates
#' @param start_time The numeric date of the root of the tree.
#' @param time_interval_size The length of the time intervals to aggregate, e.g. 1/52 corresponds to 1 week.
#' @details
#' Singletons are translated directly to importation events.
#' The time of the importation for singletons is set to the midpoint of the edge downstream of the
#' geographical transition leading to the singleton. For transmission lineages (TLs),
#' the function assumes that a TL starts with an importation event.
#' The timing of this importation event is set to the midpoint of the
#' edge on the ancestral branch of the most recent common ancestor of the TL.
#' In each TL, the function assumes that the branching points (nodes) define local transmission events.
#' Thus the procedure does not consider within-host diversity,
#' and the timing of events should be interpreted with caution.
#' @return
#' The function returns matrices with counts of local transmission ("LC") and importation events ("Import") on the time intervals, where each row contains the results from one replicate. The function also returns the times of the midpoints of the time intervals ("week_time")
#' @export

import_local_transmission = function(tree,LineageHomology_replicates,start_time,time_interval_size=1/52) {

  Result_multi = LineageHomology_replicates

  #Dates must have names that match the tip labels.
  #The matrices are way to big, but that doesn't matter in this application.
  size = 2*length(tree$tip.label)
  LC = matrix(0,nrow=ncol(Result_multi),ncol=size)
  Import = matrix(0,nrow=ncol(Result_multi), ncol=size)
  time_end = start_time+max(nodeHeights(tree))
  date_matrix = matrix(time_end+1,nrow=ncol(Result_multi), ncol=size)  ## NB. Must be generalized to take any other dataset.


  #pb <- txtProgressBar(min = 0, max =ncol(Result_multi), style = 3)
  for(i in 1:ncol(Result_multi)) {
    #setTxtProgressBar(pb, i)
    counter = 1
    for(j in 1:length(Result_multi[,i]$Lineage_sizes)) {
      #cat(paste0(" ",j))
      group_size = Result_multi[,i]$Lineage_sizes[j]
      condition = (group_size>1) #Check if group size is larger than one.

      if(condition==F) { #For groups less than one.
        Import[i,counter]=1
        date_matrix[i,counter] = Result_multi[,i]$Halfedge_over_tmrca[j] #Add the date for the event.
        counter=counter+1  #Update counter

      }

      else { #For groups larger than 1.

        #THe importation part
        Import[i,counter]=1
        date_matrix[i,counter] = Result_multi[,i]$Halfedge_over_tmrca[j] #Add dates for the events.
        counter = counter+1 #Update counter


        #The local transmission part. Time of branching points is used as dates of local transmission events.
        mrca_node = ape::getMRCA(tree, tip=c(Result_multi[,i]$Taxa_names[[j]]))
        group_taxa_names = Result_multi[,i]$Taxa_names[[j]]
        group_nodes = nodepath_quick(tree = tree, taxa = group_taxa_names) #Function is in Utiliy.R
        LTT = start_time + unlist(lapply(group_nodes, FUN = function(x) nodeheight(tree, x)))
        nltt = length(LTT)
        LC[i,counter:(counter+(nltt-1))]=(group_size-1)/group_size

        date_matrix[i,counter:(counter+(nltt-1))] = LTT #Add dates for the events.
        counter = counter+nltt #Update counter.


      }
    }
  }



  #A: Should add an extra loop to add the halfedge_above_groups larger than 1.
  #Set windows for counting importations and local transmission in.

  weeks = seq(start_time, time_end, by = time_interval_size)
  week_time=weeks[1:(length(weeks)-1)]+time_interval_size/2 #This is only for returning results at midpoint of weeks.

  date_indexes = matrix(nrow=nrow(date_matrix), ncol=ncol(date_matrix))

  for(i in 1:ncol(Result_multi)) {
    date_indexes[i,] = sort.int(date_matrix[i,],decreasing = F,index.return = T)$ix
    Import[i,] = Import[i,date_indexes[i,]] #Sort imports after dates of observed samples.
    LC[i,] = LC[i,date_indexes[i,]] #Sort local transmissions after dates of observed samples.

  }

  #Aggregate weekly importations
  t1 = table(cut(date_matrix[1,date_indexes[1,]], breaks=weeks))

  #Set up weekly matrices
  weekly_importations = matrix(0,nrow=ncol(Result_multi),ncol=length(t1))
  weekly_LC = matrix(0,nrow=ncol(Result_multi),ncol=length(t1))


  #print("Aggregating weeks ")
  for(k in 1:ncol(Result_multi)){

    #print(k/ncol(Result_multi))
    t1 = table(cut(date_matrix[k,date_indexes[k,]], breaks=weeks)) #Table that sorts observed tip dates
    j = 1
    for(i in 1:length(as.numeric(t1))){
      n = as.numeric(t1)[i]
      if(n == 0) {
        weekly_importations[k,i] = 0
        weekly_LC[k,i] = 0
        next
      }
      weekly_importations[k,i] = sum(Import[k,j:(j+n-1)])  #These are sorted, so we can just sum indexes
      weekly_LC[k,i] = sum(LC[k,j:(j+n-1)])      #These are sorted, so we can just sum indexes
      j = j+n
    }
  }


  list("Import"=weekly_importations,"LC"=weekly_LC, "Week_time"=week_time)

}



#   ____________________________________________________________________________
#   Plotting methods and summary functions                                  ####


#' Summarize_import_local_transmission
#' @description
#' This function summarizes the estimates of importation and local transmission events from
#' replicated runs of LineageHomology (e.g. by using LineageHomoolog_w_uncertainty_v2).
#' @param LineageHomology_replicates Results from running LineageHomology_w_uncertainty_v2
#'
#' @return
#' The function returns quantiles of the 2.5, 50, and 97.5 percent quantiles of the observed number of importations and local transmissions in the replicated runs.
#' @export
#'
#' @examples
Summarize_import_local_transmission = function(LineageHomology_replicates){
  imports_temp = unlist(LineageHomology_replicates[1,])[seq(1,length(LineageHomology_replicates[1,])*2,by=2)]
  local_temp = unlist(LineageHomology_replicates[1,])[seq(2,length(LineageHomology_replicates[1,])*2,by=2)]
  relative_temp = imports_temp/(imports_temp+local_temp)
  c_one = quantile(imports_temp,c(0.025,0.50,0.975))
  c_two = quantile(local_temp,c(0.025,0.50,0.975))
  c_three = quantile(relative_temp,c(0.025,0.50,0.975))
  imp_loc = rbind(c_one,c_two, c_three)
  rownames(imp_loc)=c("Import", "Local transmission", "Import / Total")
  imp_loc
}

#' plot_importation_local_transmission
#' @description
#' plot_importation_local_transmission plots the estimates of importation and local
#' transmission events that are obtained from running the import_local_transmission function.
#' @param tree tree for which the geographical estimation has been done.
#' @param result_import_local_transmission results from running import_local_transmission
#' @param time_interval_size size of the time intervals used in import_local_transmission
#' @param start_time date of the root of the phylogeny
#' @param date_breaks The interval between dates plotted on the x-axis. Takes arguments such as "1 week", "2 months","1 year" etc..
#' @param importation_or_local
#' Whether to plot importation, local transmission or both. The variable takes the arguments: "importation","local" or "both".
#' For plotting both the function requires the grid r-package.
#'
#' @return

#' @export
#'
#' @examples
plot_importation_local_transmission = function(tree,result_import_local_transmission, time_interval_size=1/52,start_time, date_breaks="1 month", importation_or_local="both")
{

  nrep= nrow(result_import_local_transmission$Import)
  time_end = start_time+max(nodeHeights(tree))
  time = seq(start_time, time_end, by = time_interval_size)
  mid_time=time[1:(length(time)-1)]+(time_interval_size/2) #This is only for returning results at midpoint of weeks.
  interval_importations = apply(result_import_local_transmission$Import,2, FUN="mean")
  interval_local_tranmissions = apply(result_import_local_transmission$LC,2, FUN="mean")
  sum_cases = interval_importations+interval_local_tranmissions

  #Plotting method
  imports = t(apply(result_import_local_transmission$Import, 2,quantile,c(0.025,0.5,0.975)))
  local_transmission = t(apply(result_import_local_transmission$LC, 2,quantile,c(0.025,0.5,0.975)))
  #Set up data
  dat=data.frame(cbind(imports,local_transmission),dates=as.Date(date_decimal(mid_time)))
  colnames(dat)=c("i2.5","i50","i97.5","e2.5","e50","e97.5","dates")

  if(importation_or_local=="both") {
    library(grid) #Requires the grid package for plotting
    g1_1 = ggplot(dat, aes(x=dates,y=i50))+geom_line(col="Blue")+theme_minimal(base_size=20)+
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
      geom_ribbon(aes(ymin=i2.5, ymax=i97.5,alpha=0.2), fill="blue")+
      scale_x_date(date_breaks=date_breaks,limits = as.Date(date_decimal(c(start_time,time_end))),labels=date_format("%d %b %Y"), expand=c(0,0))+
      ylab("Import")+
      theme(legend.position = "none")
    #g1_1


    g2_2 = ggplot(dat, aes(x=dates,y=e50))+geom_line(col="red")+theme_minimal(base_size=20)+
      theme(axis.title.x = element_blank())+
      geom_ribbon(aes(ymin=e2.5, ymax=e97.5,alpha=0.2), fill="red")+
      scale_x_date(date_breaks=date_breaks,limits = as.Date(date_decimal(c(start_time,time_end))),labels=date_format("%d %b %Y"), expand=c(0,0))+
      ylab("Local transmissions")+
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      theme(legend.position = "none")
    #g2_2
    grid.draw(rbind(ggplotGrob(g1_1), ggplotGrob(g2_2), size = "last"))
  }

  else if(importation_or_local=="importation")
  {g1_1 = ggplot(dat, aes(x=dates,y=i50))+geom_line(col="Blue")+theme_minimal(base_size=20)+
    theme(axis.title.x = element_blank())+
    geom_ribbon(aes(ymin=i2.5, ymax=i97.5,alpha=0.2), fill="blue")+
    scale_x_date(date_breaks=date_breaks,limits = as.Date(date_decimal(c(start_time,time_end))),labels=date_format("%d %b %Y"), expand=c(0,0))+
    ylab("Import")+
    theme(axis.text.x=element_text(angle=60, hjust=1))+
    theme(legend.position = "none")
  g1_1}

  else if(importation_or_local=="local") {
    g2_2 = ggplot(dat, aes(x=dates,y=e50))+geom_line(col="red")+theme_minimal(base_size=20)+
      theme(axis.title.x = element_blank())+
      geom_ribbon(aes(ymin=e2.5, ymax=e97.5,alpha=0.2), fill="red")+
      scale_x_date(date_breaks=date_breaks,limits = as.Date(date_decimal(c(start_time,time_end))),labels=date_format("%d %b %Y"), expand=c(0,0))+
      ylab("Local transmissions")+
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      theme(legend.position = "none")
    g2_2
  }

  }
