library(magrittr)
library(rms)
options(stringsAsFactors = FALSE)

split_values <- function(data1,n){
  
  out_data <- data1
  median_values <- data.frame(cancer_types=NA,features=NA,medians=NA)
  
  for (j in n:ncol(data1)) {
    median(data1[,j],na.rm = T) -> Values
    median_values[(j-1),] <- c(NA,colnames(data1)[j],Values)
    for (i in 1:nrow(data1)) {
      if(!is.na(data1[i,j])){
        if(data1[i,j] > Values){
          out_data[i,j] <- "high"
        }else{
          out_data[i,j] <- "low"
        }
      }
    }
  }
  return(list(out_data,median_values))
}


# fn_splite_by_TCGA <- function(data1,n){
#   
#   out_data <- matrix(data = NA,ncol = ncol(data1),nrow = nrow(data1))%>%as.data.frame()
#   colnames(out_data)<- colnames(data1)
#   
#   for (i in n:ncol(data1)) {
#     which(cancer_median$features[] == colnames(data1)[i])->ID
#     
#     if(length(ID)>0){
#       for (j in 1:nrow(data1)) {
#         if(data1[j,i] > cancer_median$medians[ID]){
#           out_data[j,i] <- "high"
#         }else{
#           out_data[j,i] <- "low"
#         }
#       }
#     }
#     
#   }
#   
#   out_data[,1:(n-1)] <- data1[,1:(n-1)]
#   out_data
# }


readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/muli_uni_survival/OS_multi_res.tsv")%>%
  dplyr::filter(coxp < 0.1)%>%
  dplyr::filter(cancer_types == "ACC")%>%
  dplyr::select(Features)%>%
  as.matrix()%>%as.character()->features
gsub("_high","",features)-> features

readr::read_rds("/home/liull/TCGA_nomograph/Filtered_GEO_clinical_TIL_expr.rds.gz")%>%
  dplyr::filter(Project_ID == "GSE33371")-> cancer_rds
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/OS_Median_values.txt")%>%
  dplyr::filter(cancer_types == "ACC")->cancer_median

cancer_rds[[3]][[1]]%>%
  as.data.frame()%>%
  dplyr::select(Sample_ID,OS_days,OS_status)->survivals

cancer_rds[[4]][[1]]%>%
  as.data.frame()%>%
  dplyr::select(Sample_ID,Tumor_stage)->stages

cancer_rds[[5]][[1]]%>%
  as.data.frame()%>%
  dplyr::select(Sample_ID,intersect(features,colnames(cancer_rds[[5]][[1]])))->Infiltrations

gsub("-", ".",cancer_rds[[6]][[1]]$Symbol,fixed = TRUE) ->cancer_rds[[6]][[1]]$Symbol
cancer_rds[[6]][[1]]%>%
  as.data.frame()%>%
  dplyr::filter(Symbol %in% intersect(features,cancer_rds[[6]][[1]]$Symbol))->exprs
rownames(exprs) <- exprs$Symbol
exprs[,-1]%>%t()%>%as.data.frame()%>%tibble::rownames_to_column()->exprs

merge(survivals,stages)%>%
  merge(Infiltrations)%>%
  merge(exprs,by.x = "Sample_ID",by.y = "rowname")->all_infor

# gsub("^1$","I",all_infor$Tumor_stage) %>%
#   gsub("^2$","II",.) %>% gsub("^3","III",.)%>% gsub("^4$","IV",.) -> all_infor$Tumor_stage

#splite by GEO's own median
split_values(all_infor,5) ->split_res
split_res[[1]]->all_infor
split_res[[2]]->median_values
#splite by TCGA's median
fn_splite_by_TCGA(all_infor,5)->all_infor


all_infor[,-1] ->all_infor
colnames(all_infor)[1:3] <- c("OS", "Status","Stage")
all_infor$OS <- as.numeric(all_infor$OS)
all_infor$Status <- as.numeric(all_infor$Status)
#cbind(all_infor,HLA.DRB1 = rep("low",nrow(all_infor)))->all_infor
ddist <- datadist(all_infor)
options(datadist='ddist')
f_cph_0 <- cph(Surv(OS,Status) ~., data=all_infor,surv=TRUE,x=TRUE, y=TRUE,time.inc=3650)
1-rcorrcens(Surv(OS,Status) ~ predict(f_cph_0), data =  all_infor)[[1]]

