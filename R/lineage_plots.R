library(lubridate)
library(ggridges)
library(forcats)
library(scales)

#' Title
#'
#' @param Result_Lineage_homology The output from linege homology.
#' @param name_date #Data frame with taxa names in one column and numeric dates of sampling in another. The column names should be c("name", "date")
#'
#' @return
#' @export
#'
#' @examples
lineage_info = function(Result_Lineage_homology, name_date) {

  result_matrix = matrix(nrow=sum(Result_Lineage_homology$Lineage_sizes)+length(Result_Lineage_homology$Lineage_sizes), ncol=5)
  counter = 1 #Used for indexing correctly.
  for(i in 1:length(Result_Lineage_homology$Taxa_names)) {
    numeric_dates = c(name_date$date[which(name_date$name%in%unlist(Result_Lineage_homology$Taxa_names[i]))], Result_Lineage_homology$`MRCA's`[i])
    group_identity = rep(i, length(numeric_dates))
    group_size = Result_Lineage_homology$Lineage_sizes[i]
    range_group = rep(range(numeric_dates)[2]-range(numeric_dates)[1],length(numeric_dates))
    lineage_state=rep(Result_Lineage_homology$lineage_state[i],length(numeric_dates))
    this_lineage = cbind(numeric_dates,group_identity,group_size,range_group,lineage_state) #The dates has one more entry than group size because of the mrca.
    result_matrix[counter:(counter+length(unlist(Result_Lineage_homology$Taxa_names[i]))),]=this_lineage
    counter = counter+length(unlist(Result_Lineage_homology$Taxa_names[i]))+1
  }

  colnames(result_matrix)=c("dates", "group_no","group_size","range_group", "lineage_state")
  result_matrix=as.data.frame(result_matrix)
  result_matrix$group_size = as.numeric(result_matrix$group_size)
  result_matrix$group_no = as.factor(result_matrix$group_no)
  result_matrix$range_group = as.numeric(result_matrix$range_group)
  result_matrix$dates = as.numeric(result_matrix$dates)
  result_matrix$lineage_state = as.factor(result_matrix$lineage_state)

  result_matrix
}

#' Title
#'
#' @param Result_LineageHomology Result from the function LineageHomology
#'
#' @return
#' @export
#'
#' @examples
treemap_lineagehomology = function(Result_LineageHomology) {
  atreemap = data.frame(group=paste0("G:",1:length(Result$Lineage_sizes)," S: ", Result_LineageHomology$Lineage_sizes, "\n TMRCA: ", round(Result_LineageHomology$`MRCA's`,2)),Value=Result_LineageHomology$Lineage_sizes)
  treemap::treemap(atreemap,index=("group"),vSize="Value", type="index", palette="Paired")
}


#' Title
#'
#' @param Result_lineage_info
#' @param groups_larger_than
#' @param datelims
#'
#' @return
#' @export
#'
#' @examples
ridgeplot_lineagedensities = function(Result_lineage_info, datelims,groups_larger_than=4,color_by_state=FALSE) {
  #Requires dplyr, scales
  dateupplow = as.POSIXct(strptime(c(paste0(datelims[1]," 03:00"),paste0(datelims[2]," 16:00")),format = "%Y-%m-%d %H:%M"))
  Result_lineage_info= Result_lineage_info %>% dplyr::mutate(group_no=forcats::fct_reorder(group_no, group_size)) #Reorder group name factor levels by group size.
  Result_lineage_info=Result_lineage_info[as.numeric(Result_lineage_info$group_size)>groups_larger_than,]
  Result_lineage_info$dates=date_decimal(Result_lineage_info$dates)
  if(color_by_state == TRUE) {
    ggplot(Result_lineage_info,aes(x=dates, y=group_no,fill=lineage_state))+
      ggridges::geom_density_ridges(aes(point_color = lineage_state, point_fill = lineage_state),alpha=0.4,point_alpha=0.4,scale = 2.5, size = 0.25, rel_min_height = 0.03,jittered_points = TRUE)+
      #ggridges::geom_density_ridges(alpha=0.4,scale = 2.5, size = 0.25, rel_min_height = 0.03,point_shape = "|", point_size = 1,jittered_points = TRUE,position = ggridges::position_points_jitter(height = 0))+
      ggridges::theme_ridges(font_size = 10)+scale_y_discrete(labels=sort(aggregate(group_size~group_no, data=Result_lineage_info, FUN=function(x) c(mean=mean(x), count=length(x)))[[2]][,1]))+
      scale_x_datetime(date_breaks = datelims[3],limits =  dateupplow,labels=date_format("%d %b %Y"))+ylab("Group size")+
      scale_fill_manual(values = c("red","blue"))+theme(legend.position ="none")+
      ggplot2:::manual_scale("point_color", values=c("red","blue"), guide = "none")
  }
  else{
    ggplot(Result_lineage_info,aes(x=dates, y=group_no))+
      ggridges::geom_density_ridges(alpha=0.4,point_alpha=0.4,scale = 2.5, size = 0.25, rel_min_height = 0.03,jittered_points = TRUE)+
      ggridges::theme_ridges(font_size = 10,)+scale_y_discrete(labels=sort(aggregate(group_size~group_no, data=Result_lineage_info, FUN=function(x) c(mean=mean(x), count=length(x)))[[2]][,1]))+
      scale_x_datetime(date_breaks = datelims[3],limits =  dateupplow,labels=date_format("%d %b %Y"))+ylab("Group size")
  }

}



#' Title
#'
#' @param Result_lineage_info The output from lineage info.
#' @param datelims Limits of the dates in format c("yyyy-mm-dd","yyyy-mm-dd", ""). The last argument can be e.g. "1 week"/"1 year"/"1 month"
#' @param color_by_state Color by the state of the lineage.
#'
#' @return
#' @export
#'
#' @examples
lineage_growth_cumulative = function(Result_lineage_info,datelims,color_by_state=FALSE) {
  dateupplow = as.POSIXct(strptime(c(paste0(datelims[1]," 03:00"),paste0(datelims[2]," 16:00")),format = "%Y-%m-%d %H:%M"))
  Result_lineage_info$group_cumsum = ave(Result_lineage_info$dates,Result_lineage_info$group_no, FUN=function(x) rank(x,ties.method = "first"))-1 #Use the rank function to obtain cumulative counts. The -1 removes the count of the mrca.

  if(color_by_state == TRUE) {
    g1 = ggplot(Result_lineage_info,aes(x=date_decimal(dates), y=group_cumsum, group=group_no,color=lineage_state))+geom_line(alpha=0.4)+theme_bw()+
      scale_x_datetime(date_breaks = datelims[3],limits =  dateupplow,labels=date_format("%d %b %Y"))+ylab("Group size")+xlab(element_blank())+scale_color_manual(values = c("red","blue"))+
      theme(legend.position ="none")
  }
  else {
    g1 = ggplot(Result_lineage_info,aes(x=date_decimal(dates), y=group_cumsum,group=group_no))+geom_line(alpha=0.4)+theme_bw()+
    scale_x_datetime(date_breaks = datelims[3],limits =  dateupplow,labels=date_format("%d %b %Y"))+ylab("Group size")+xlab(element_blank())
    }
  g1
}



#Example of usage:
# name_date_category = data.frame(metadat$Strain, decimal_date(metadat$`Collection Data`),metadat$`Pangolin clade`)
# colnames(name_date_category)=c("name", "date", "category")
# Result_Lineage_Homology = Result
# dummy_matrix = lineage_densities_B117(Result, name_date_category)
# ridgeplot_lineagedensities(dummy_matrix,4)
