#Given a result output, whether if it belongs to a group size that is larger than one.
#Naively calculate the importation weight by the local group size
#E.G: If the group size is 2, then it is +0.5 importations, and +0.5 local transmissions.
#E.G: IF the group size is 30 then it is +1/30 importation and + 29/30 local transmissions


#Load metadata Southafrican: Only for testing
# library(data.table)
# metadata = fread("~/Dropbox/Covid/Southafrican_03_24/nextstrain_groups_niph_ncov_southafrican-2021-04-07_metadata.tsv")
# sum(tree$tip.label %in% metadata$Strain)
# length(tree$tip.label)
# tip_list = Norwegian_tips
# dates = decimal_date(metadat$`Collection Data`)
# names(dates)=metadata$Strain

#Load metadata for the UK-variant.
# library(data.table)
# metadata = fread("~/Dropbox/Covid/UK_03_24/nextstrain_groups_niph_ncov_UK-centric-2021-04-07_metadata.tsv")
# sum(tree$tip.label %in% metadata$Strain)
# length(tree$tip.label)
# tip_list = Norwegian_tips
# #Setup dates.
# dates = decimal_date(metadata$`Collection Data`)
# names(dates)=metadata$Strain



# @ Takes argument dates, should have names matching tip labels from tree.
# @ Takes the phylogenetic tree used for ancestral state reconstruction.
# @ Takes the Result_multi output from and replicated count script on castor (which runs on ace, simmap etc.)
# @ Takes a list of taxa which to compute the importation load for: e.g. all tips in Norway 

Relative_load_import_multi = function(tree,multicount, dates,start_time) {
  Result_multi = multicount
  
  #Dates must have names that match the tip labels.
  LC = matrix(nrow=ncol(Result_multi),ncol=length(tip_list))
  Import = matrix(nrow=ncol(Result_multi), ncol=length(tip_list))
  date_matrix = matrix(nrow=ncol(Result_multi), ncol=length(tip_list))
  
  print("Collecting import dates ")
  for(i in 1:ncol(Result_multi)) {
    print(i)/ncol(Result_multi)
    for(j in 1:length(tip_list)){
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
#  
# start_time = start_date
# #Set windows for counting importations and local transmission in.
# start_time 
# time_end = start_time+max(nodeHeights(tree))
# weeks = seq(start_time, time_end, by = 1/52)
# week_time=weeks[1:(length(weeks)-1)]+1/(52*2) #This is only for returning results at midpoint of weeks.
# weekly_importations = apply(Result_with_uncertainty$Import,2, FUN="mean")
# weekly_local_tranmissions = apply(Result_with_uncertainty$LC,2, FUN="mean")
# sum_cases = weekly_importations+weekly_local_tranmissions
# 
# 
# #This uses the mean and discards uncertainty.
# ggdata = data.frame(dates = as.Date(date_decimal(week_time)),
#                   fraction=c(weekly_importations/sum_cases,weekly_local_tranmissions/sum_cases),
#                   group=c(rep("Import",length(week_time)),
#                           rep("Local transmission",length(week_time))))
# 
# pdf(paste0(path_to_results,"relative_load.pdf"), width=15,height=7.5)
# g1 = ggplot(data=ggdata, aes(x=dates, y=fraction, fill=group)) + 
#   geom_area(alpha=0.6 , size=1, colour="black")+
#   theme_minimal(base_size=20)+scale_fill_brewer(palette="Set1")+
#   scale_x_date(date_breaks="1 week")+
#   theme(axis.text.x=element_text(angle=60, hjust=1))+
#   ylab("")+xlab("")+
#   annotate("rect",xmin=as.Date(decimal2Date(time_end-1/12)),xmax=as.Date(dat$dates[length(dat$dates)]),ymin=-Inf,ymax=Inf,alpha=0.2)+
#   theme(legend.position = c(0.55,0.15), legend.title=element_blank())
# g1
# dev.off()
# 
# #Plotting method
# imports = t(apply(Result_with_uncertainty$Import, 2,quantile,c(0.025,0.5,0.975)))
# exports = t(apply(Result_with_uncertainty$LC, 2,quantile,c(0.025,0.5,0.975)))
# dat=data.frame(cbind(imports,exports),dates=decimal2Date(week_time))
# colnames(dat)=c("i2.5","i50","i97.5","e2.5","e50","e97.5","dates")
# 
# 
# #Line for UK
# #scale_x_date(date_breaks="1 week",limits =c(as.Date("2020-11-30"), as.Date("2021-03-15")))+ # limits=c(as.Date("2020-06-30"), as.Date("2021-03-29"))
# 
# svglite(paste0(path_to_results,"Imports.svg"), width=10,height=5)
# g1 = ggplot(dat, aes(x=dates,y=i50))+geom_line(col="Red")+theme_minimal(base_size=20)+
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
#   ylim(0,5)+
#   geom_ribbon(aes(ymin=i2.5, ymax=i97.5,alpha=0.2), fill="lightblue")+
#   scale_x_date(date_breaks="1 week",limits =c(as.Date("2020-11-30"), as.Date("2021-03-28")),expand=c(0,0))+ # limits=c(as.Date("2020-06-30"), as.Date("2021-03-29"))
#   #theme(axis.text.x=element_text(angle=60, hjust=1))+scale_fill_continuous("Grey")+xlab("")+
#   ylab("Import")+
#   theme(legend.position = "none")+annotate("rect",xmin=as.Date(decimal2Date(time_end-1/12)),xmax=as.Date(dat$dates[length(dat$dates)]),ymin=-Inf,ymax=Inf,alpha=0.2)
# g1
# dev.off()
# 
# 
# #Line for UK
# #scale_x_date(date_breaks="1 week",limits =c(as.Date("2020-11-30"), as.Date("2021-03-15")))+ #limits =c(as.Date("2021-01-18"), as.Date("2021-03-15"))
# 
# svglite(paste0(path_to_results,"Local_transmissions.svg"), width=10,height=5)
# g2 = ggplot(dat, aes(x=dates,y=e50))+geom_line(col="Darkgreen")+theme_minimal(base_size=20)+ylim(0,45)+
#   scale_x_date(date_breaks="1 week",limits =c(as.Date("2020-11-30"), as.Date("2021-03-28")), expand=c(0,0))+ #limits =c(as.Date("2021-01-18"), as.Date("2021-03-15"))
#   geom_ribbon(aes(ymin=e2.5, ymax=e97.5,alpha=0.2), fill="lightblue")+
#   theme(axis.text.x=element_text(angle=60, hjust=1))+scale_fill_continuous("Grey")+ylab("Local transmission")+xlab("")+
#   theme(legend.position = "none")+annotate("rect",xmin=as.Date(decimal2Date(time_end-1/12)),xmax=as.Date(dat$dates[length(dat$dates)]),ymin=-Inf,ymax=Inf,alpha=0.2)
# g2
# dev.off()
# 
# 
# graphics.off()
# pdf(paste0(path_to_results,"Import_export.pdf"),width=7.5,height=7.5)
# #grid.arrange(grobs=list(g1,g2),layout_matrix=rbind(c(1,1),c(2,2)))
# #grid.arrange(grobs=list(g1,g2),ncol=1)
# grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
# dev.off()
