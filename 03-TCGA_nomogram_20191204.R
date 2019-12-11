#####get the best ICP genes' ssGSEA score and TIL features(self adjusted and pure tumor/fixed tumor)
##do not define high/low,use the number to do the survival analysis
##use the uni survival result

rm(list=ls())
library(magrittr)
options(stringsAsFactors = FALSE)

#get the best features in TIL\gene_sets seperatly---------------------------------------------------------------
fn_nomogram_model <- function(data){
  
  if(max(data$OS,na.rm = TRUE) >= 1825){
    f_cph <- cph(Surv(OS,Status) ~., data=data,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825,singular.ok = TRUE)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  data)[[1]]%>%signif(3)->C_index
    
  }else{
    f_cph <- cph(Surv(OS,Status) ~., data=data,surv=TRUE,x=TRUE, y=TRUE,time.inc=1095,singular.ok = TRUE)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  data)[[1]] ->C_index
    
  }
  
  return(list(f_cph,C_index))
}

readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/muli_uni_survival_20191204/OS_pure_uni_res.sig.tsv")->OS_pure_uni_sig

TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_fixed_TIL_ssGSEA_data_20191204.rds.gz")%>%
  dplyr::filter(cancer_types %in% unique(OS_pure_uni_sig$cancer_types))

#TNM stage
single_C_index <- data.frame(cancer_types=NA,TNM_C_index=NA)
for (i in 1:nrow(TCGA_all)) {
  
  TCGA_all$OS_stage[[i]]%>%as.data.frame()%>%
    dplyr::select(-Age)%>%
    dplyr::filter(Stage != "NA")%>%
    dplyr::select(-barcode)-> OS_stage
  OS_stage$OS=as.numeric(OS_stage$OS)
  OS_stage$Status=as.numeric(OS_stage$Status)
  
  if(nrow(OS_stage) > 0){
    ddist <- datadist(OS_stage)
    options(datadist ='ddist')
    single_C_index[i,]<- c(TCGA_all$cancer_types[i],fn_nomogram_model(OS_stage)[[2]])
  }else{
    single_C_index[i,]<- c(TCGA_all$cancer_types[i],"NA")
  }
  
}
readr::write_tsv(single_C_index,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model_20191204/TNM_C_index.tsv")



#gene set
single_C_index <- data.frame(cancer_types=NA,all_gene_set_features=NA,all_C_index=NA,new_gene_set_features=NA,new_gene_set_C_index=NA)
for (i in 1:nrow(TCGA_all)) {
  
  dplyr::filter(OS_pure_uni_sig,cancer_types == TCGA_all$cancer_types[i])%>%
    dplyr::select(features)%>%as.matrix()%>%as.character()%>%unique()->i_all_features
  
  TCGA_all$OS_stage[[i]]%>%as.data.frame()%>%
    dplyr::select(-Age,-Stage)-> OS_stage
  OS_stage$OS=as.numeric(OS_stage$OS)
  OS_stage$Status=as.numeric(OS_stage$Status)
  
  TCGA_all$ssGSEA[[i]]%>%as.data.frame()-> ssGSEA
  dplyr::select(ssGSEA,barcode,intersect(i_all_features,colnames(ssGSEA)))-> ssGSEA
  
  dplyr::inner_join(OS_stage,ssGSEA,by="barcode")%>%dplyr::select(-barcode)->all_ssGSEA
  
  if(ncol(all_ssGSEA) > 3){
    ddist <- datadist(all_ssGSEA)
    options(datadist ='ddist')
    fn_nomogram_model(all_ssGSEA)[[2]]->C_index_all
    
    Max_C_index <- C_index_all
    Max_rm_ID <- character()
    for (j in 1:(ncol(all_ssGSEA)-3)) {
      combn(3:ncol(all_ssGSEA),j) ->rm_ID
      for (m in 1:ncol(rm_ID)) {
        all_ssGSEA[,-rm_ID[,m]]->ssGSEA_m
        fn_nomogram_model(ssGSEA_m)[[2]]->C_index_m
        if(C_index_m > Max_C_index){
          Max_C_index <- C_index_m
          Max_rm_ID <- rm_ID[,m]
        }
        
      }
      
    }
    if(length(Max_rm_ID) != 0){
      all_ssGSEA[,-Max_rm_ID]->ssGSEA_Max
      fn_nomogram_model(ssGSEA_Max)[[2]]->C_index_Max
    }else{
      all_ssGSEA->ssGSEA_Max
      C_index_all->C_index_Max
    }
    
    single_C_index[i,] <- c(TCGA_all$cancer_types[i],paste(colnames(all_ssGSEA)[-c(1,2)],collapse = " + "),C_index_all,
                            paste(colnames(ssGSEA_Max)[-c(1,2)],collapse = " + "),C_index_Max)
  }else if(ncol(all_ssGSEA) == 3){
    ddist <- datadist(all_ssGSEA)
    options(datadist ='ddist')
    fn_nomogram_model(all_ssGSEA)[[2]]->C_index_all
    single_C_index[i,] <- c(TCGA_all$cancer_types[i],paste(colnames(all_ssGSEA)[-c(1,2)],collapse = " + "),C_index_all,
                            paste(colnames(all_ssGSEA)[-c(1,2)],collapse = " + "),C_index_all)
  }else{
    single_C_index[i,] <- c(TCGA_all$cancer_types[i],"NA","NA","NA","NA")
  }
  
}
readr::write_tsv(single_C_index,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model_20191204/ssGSEA_C_index.tsv")

#TIL
single_C_index <- data.frame(cancer_types=NA,all_TIL_features=NA,all_C_index=NA,new_TIL_features=NA,new_TIL_C_index=NA)
for (i in 1:nrow(TCGA_all)) {
  
  dplyr::filter(OS_pure_uni_sig,cancer_types == TCGA_all$cancer_types[i])%>%
    dplyr::select(features)%>%as.matrix()%>%as.character()%>%unique()->i_all_features
  
  TCGA_all$OS_stage[[i]]%>%as.data.frame()%>%
    dplyr::select(-Age,-Stage)-> OS_stage
  OS_stage$OS=as.numeric(OS_stage$OS)
  OS_stage$Status=as.numeric(OS_stage$Status)
  
  TCGA_all$pure_tumor_Infiltration[[i]]%>%as.data.frame()-> pure_TIL
  dplyr::select(pure_TIL,barcode,intersect(i_all_features,colnames(pure_TIL)))-> pure_TIL
  
  dplyr::inner_join(OS_stage,pure_TIL,by="barcode")%>%dplyr::select(-barcode)->all_pure_TIL
  
  if(ncol(all_pure_TIL) > 3){
    ddist <- datadist(all_pure_TIL)
    options(datadist ='ddist')
    fn_nomogram_model(all_pure_TIL)[[2]]->C_index_all
    
    Max_C_index <- C_index_all
    Max_rm_ID <- character()
    for (j in 1:(ncol(all_pure_TIL)-3)) {
      combn(3:ncol(all_pure_TIL),j) ->rm_ID
      for (m in 1:ncol(rm_ID)) {
        all_pure_TIL[,-rm_ID[,m]]->TIL_m
        fn_nomogram_model(TIL_m)[[2]]->C_index_m
        if(C_index_m > Max_C_index){
          Max_C_index <- C_index_m
          Max_rm_ID <- rm_ID[,m]
        }
        
      }
      
    }
    if(length(Max_rm_ID) != 0){
      all_pure_TIL[,-Max_rm_ID]->TIL_Max
      fn_nomogram_model(TIL_Max)[[2]]->C_index_Max
    }else{
      all_pure_TIL->TIL_Max
      C_index_all->C_index_Max
    }
    
    single_C_index[i,] <- c(TCGA_all$cancer_types[i],paste(colnames(all_pure_TIL)[-c(1,2)],collapse = " + "),C_index_all,
                            paste(colnames(TIL_Max)[-c(1,2)],collapse = " + "),C_index_Max)
  }else if(ncol(all_pure_TIL) == 3){
    ddist <- datadist(all_pure_TIL)
    options(datadist ='ddist')
    fn_nomogram_model(all_pure_TIL)[[2]]->C_index_all
    single_C_index[i,] <- c(TCGA_all$cancer_types[i],paste(colnames(all_pure_TIL)[-c(1,2)],collapse = " + "),C_index_all,
                            paste(colnames(all_pure_TIL)[-c(1,2)],collapse = " + "),C_index_all)
  }else{
    single_C_index[i,] <- c(TCGA_all$cancer_types[i],"NA","NA","NA","NA")
  }
  
}
readr::write_tsv(single_C_index,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model_20191204/pure_TIL_C_index.tsv")


#compare effect TNM~gene_set~TIL~gene_set+TIL.. use the best features---------------------------------------------------------
