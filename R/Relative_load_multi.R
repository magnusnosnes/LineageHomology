#Given a result output, whether if it belongs to a group size that is larger than one.
#Naively calculate the importation weight by the local group size
#E.G: If the group size is 2, then it is +0.5 importations, and +0.5 local transmissions.
#E.G: IF the group size is 30 then it is +1/30 importation and + 29/30 local transmissions

# library(data.table)

#Load metadata Southafrican: Only for testing
# metadata = fread("~/Dropbox/Covid/Southafrican_03_24/nextstrain_groups_niph_ncov_southafrican-2021-04-07_metadata.tsv")
# sum(tree$tip.label %in% metadata$Strain)
# length(tree$tip.label)
# tip_list = Norwegian_tips
# dates = decimal_date(metadat$`Collection Data`)
# names(dates)=metadata$Strain

#Load metadata for the UK-variant.

# metadata = fread("~/Dropbox/Covid/UK_03_24/nextstrain_groups_niph_ncov_UK-centric-2021-04-07_metadata.tsv")
# sum(tree$tip.label %in% metadata$Strain)
# length(tree$tip.label)
# tip_list = Norwegian_tips
# #Setup dates.
# dates = decimal_date(metadata$`Collection Data`)
# names(dates)=metadata$Strain


library(phytools)
#' Title
#'
#' @param tree Takes the phylogenetic tree used for ancestral state reconstruction.
#' @param multicount Takes the Result_multi output from and replicated count script on castor (which runs on ace, simmap etc.)
#' @param dates Takes argument dates, should have names matching tip labels from tree.
#' @param start_time Takes the date of the root of the phylogenetic tree.
#' @param tip_list Takes a list of taxa which to compute the importation load for: e.g. all tips in Norway
#'
#' @return
#' @export
#' @importFrom phytools nodeHeights
#' @examples
Relative_load_import_multi = function(tree,multicount, dates,start_time,tip_list) {
  Result_multi = multicount

  #Dates must have names that match the tip labels.
  LC = matrix(nrow=ncol(Result_multi),ncol=length(tip_list))
  Import = matrix(nrow=ncol(Result_multi), ncol=length(tip_list))
  date_matrix = matrix(nrow=ncol(Result_multi), ncol=length(tip_list))

  print("Collecting import dates ")
  for(i in 1:ncol(Result_multi)) {
    print(i)/ncol(Result_multi)
    for(j in 1:length(tip_list)){
      print(j)
      tip = tip_list[j]
      date_new = dates[which(names(dates)==tip)]
      if(length(date_new)==0) next
      group_no= which(lapply(Result_multi[,i]$Taxa_names, FUN=function(x) tip%in%x)==T)
      #group_no = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))
      group_size = Result_multi[,i]$Lineage_sizes[group_no]
      Import[i,j] = 1/group_size
      LC[i,j] =  1-1/group_size
      date_matrix[i,j] = date_new
    }
  }

  #Set windows for counting importations and local transmission in.
  start_time
  time_end = start_time+max(nodeHeights(tree))
  weeks = seq(start_time, time_end, by = 1/52)
  week_time=weeks[1:(length(weeks)-1)]+1/(52*2) #This is only for returning results at midpoint of weeks.

  date_indexes = matrix(nrow=ncol(Result_multi), ncol=length(tip_list))

  for(i in 1:ncol(Result_multi)) {
    date_indexes[i,] = sort.int(date_matrix[i,],decreasing = F,index.return = T)$ix
    Import[i,] = Import[i,date_indexes[i,]] #Sort imports after dates of observed samples.
    LC[i,] = LC[i,date_indexes[i,]] #Sort local transmissions after dates of observed samples.
  }

  #Aggregate weekly importations
  t1 = table(cut(date_matrix[1,date_indexes[1,]], breaks=weeks))
  #Set up weekly matrices
  weekly_importations = matrix(nrow=ncol(Result_multi),ncol=length(t1))
  weekly_LC = matrix(nrow=ncol(Result_multi),ncol=length(t1))


  print("Aggregating weeks ")
  for(k in 1:ncol(Result_multi)){

    print(k)/ncol(Result_multi)
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


#
# #Example of usage:
# Result_with_uncertainty = Relative_load_import_multi(tree, multi_counts,dates, start_date)
