options(stringsAsFactors = FALSE)


#all GEO expr,TIL,survival,stage---------------------------------------------------------------------------
readxl::read_xlsx("/home/liull/TCGA_nomograph/GEO_cancer_metadata.xlsx",sheet = "Sheet1")%>%
  dplyr::select(-Cancer_type,-Data_type,-Sex,-TNM_stage_T,-TNM_stage_N,-TNM_stage_M)->metadata

#tidyr::nest expr and TIL
readxl::read_xlsx("/home/liull/reference/All_EntrezID_Symbl_NCBI.xlsx")%>%as.data.frame()->relationship
metadata$TCGA_cancer_type %>% unique()->all_Cancers
for (j in 1:length(all_Cancers)) {
  dplyr::filter(metadata,TCGA_cancer_type == all_Cancers[j])%>%
    dplyr::select(Project_ID)%>%
    unique()%>%as.matrix()%>%as.character()->j_Projects
  
  for (i in 1:length(j_Projects)) {
    list.files("/home/liull/TCGA_nomograph/GEO_expr_data",pattern = j_Projects[i],full.names = TRUE) ->i_expr_path
    read.table(i_expr_path,header = T,sep = "\t",quote = "")%>%
      tibble::rownames_to_column()%>%
      merge(relationship,.,by.x = "Symbol",by.y = "rowname")->i_expr
    cbind(Project_ID= rep(j_Projects[i],nrow(i_expr)),i_expr)%>%
      cbind(TCGA_cancer_type= rep(all_Cancers[j],nrow(i_expr)),.)->i_expr
    tidyr::gather(i_expr,"Samples","Values",colnames(i_expr)[-c(1:4)])->i_expr
    
      if(i == 1 & j==1){
        all_exprs <- i_expr
      }else{
        rbind(all_exprs,i_expr)->all_exprs
      }
    
  }
  
  
  for (m in 1:length(j_Projects)) {
    list.files("/home/liull/TCGA_nomograph/GEO_TIL_data",pattern = j_Projects[m],full.names = TRUE) ->m_TIL_path
    read.table(m_TIL_path,header = T,sep = "\t",quote = "")%>%
      tibble::rownames_to_column()->m_TIL
    colnames(m_TIL)[1] <- "Sample_ID"
    cbind(Project_ID= rep(j_Projects[m],nrow(m_TIL)),m_TIL)%>%
      cbind(TCGA_cancer_type= rep(all_Cancers[j],nrow(m_TIL)),.)->m_TIL
    
    if(m == 1 & j==1){
      all_TILs <- m_TIL
    }else{
      rbind(all_TILs,m_TIL)->all_TILs
    }
    
  }
  
  
}
  
tidyr::nest(all_exprs,-TCGA_cancer_type,-Project_ID)->Exprs
dplyr::mutate(Exprs,expr = purrr::map2(data,Project_ID,.f=function(.x,.y){ 
  print(.y) 
  .x %>% tidyr::spread(Samples,Values)
}))%>%
  dplyr::select(-data)->Exprs  
colnames(Exprs)[3] <- "Exprs"


tidyr::nest(all_TILs,-TCGA_cancer_type,-Project_ID)->TILs
colnames(TILs)[3] <- "TILs"

  
#tidyr::nest survival & stage
dplyr::select(metadata,-Tumor_stage) %>% 
  tidyr::nest(-TCGA_cancer_type, -Project_ID)->Survivals
colnames(Survivals)[3] <- "Survivals" 
dplyr::select(metadata,TCGA_cancer_type,Project_ID,Sample_ID,Tumor_stage) %>% 
  tidyr::nest(-TCGA_cancer_type, -Project_ID)->Stages
colnames(Stages)[3] <- "Stages"


merge(Survivals,Stages, by = c("TCGA_cancer_type","Project_ID"))%>%
  merge(TILs,by = c("TCGA_cancer_type","Project_ID"))%>%
  merge(Exprs,by = c("TCGA_cancer_type","Project_ID"))->clinical_info 
 
clinical_info %>% 
  readr::write_rds("/home/liull/TCGA_nomograph/GEO_clinical_TIL_expr.rds.gz")



# inhibit_cells <- c("Exhausted","iTreg","Neutrophil","Monocyte","nTreg","Tr1")
# non_inhibit_TIL_markers <- readr::read_tsv("/home/liull/TCGA_nomograph//markers_used_to_predict_TIL_miao_TCAP.txt")%>% 
#   tidyr::gather(key="Cell_type",value="markers")%>% 
#   dplyr::filter(!is.na(markers)) %>% 
#   dplyr::filter(! Cell_type %in% inhibit_cells)
# 
# checkpoints_as_TILMarker <- readr::read_tsv("/home/liull/TCGA_nomograph/ICPs_all_info_class.tsv") %>% 
#   dplyr::filter(! symbol %in% non_inhibit_TIL_markers$markers) 



#  searched gene list related immuncheckpoint
checkpoints <- readr::read_tsv("/home/liull/TCGA_nomograph/ICPs_all_info_class.tsv") 

# expression data 
genelist_exp <- readr::read_rds("/home/liull/TCGA_nomograph/GEO_clinical_TIL_expr.rds.gz") %>% 
  dplyr::mutate(exp_filter = purrr::map2(Exprs,Project_ID,.f=function(.x,.y){
    
    print(.y) 
    .x %>% 
      dplyr::filter(Symbol %in% checkpoints$symbol) 
    
  })) %>%
  dplyr::select(-Exprs)  

genelist_exp %>%
  readr::write_rds("/home/liull/TCGA_nomograph/Filtered_GEO_clinical_TIL_expr.rds.gz")





#GEO data with TNM stage-------------------------------------------------------------------
readxl::read_xlsx("/home/liull/TCGA_nomograph/GEO_TNM_metadata.xlsx",sheet = "Sheet1")%>%
  dplyr::select(-Cancer_type,-Data_type,-Sex,-TNM_stage_T,-TNM_stage_N,-TNM_stage_M)%>%
  dplyr::filter(Tumor_stage != "NA")%>%
  dplyr::filter(OS_status != "NA")->metadata

#tidyr::nest expr and TIL
readxl::read_xlsx("/home/liull/reference/All_EntrezID_Symbl_NCBI.xlsx")%>%as.data.frame()->relationship
metadata$TCGA_cancer_type %>% unique()->all_Cancers
for (j in 1:length(all_Cancers)) {
  dplyr::filter(metadata,TCGA_cancer_type == all_Cancers[j])%>%
    dplyr::select(Project_ID)%>%
    unique()%>%as.matrix()%>%as.character()->j_Projects
  
  for (i in 1:length(j_Projects)) {
    metadata %>%
      dplyr::filter(Project_ID == j_Projects[i])%>%
      dplyr::select(Sample_ID)%>%
      unique()%>%as.matrix()%>%as.character()->i_samples
    
    list.files("/home/liull/TCGA_nomograph/GEO_expr_data",pattern = j_Projects[i],full.names = TRUE) ->i_expr_path
    read.table(i_expr_path,header = T,sep = "\t",quote = "")%>%
      tibble::rownames_to_column()%>%
      merge(relationship,.,by.x = "Symbol",by.y = "rowname")%>%
      dplyr::select(Symbol,GeneID,i_samples)->i_expr
    cbind(Project_ID= rep(j_Projects[i],nrow(i_expr)),i_expr)%>%
      cbind(TCGA_cancer_type= rep(all_Cancers[j],nrow(i_expr)),.)->i_expr
    tidyr::gather(i_expr,"Samples","Values",colnames(i_expr)[-c(1:4)])->i_expr

    
    
    list.files("/home/liull/TCGA_nomograph/GEO_TIL_data",pattern = j_Projects[i],full.names = TRUE) ->i_TIL_path
    read.table(i_TIL_path,header = T,sep = "\t",quote = "")%>%
      tibble::rownames_to_column()%>%
      dplyr::filter(rowname %in% i_samples)->i_TIL
    colnames(i_TIL)[1] <- "Sample_ID"
    cbind(Project_ID= rep(j_Projects[i],nrow(i_TIL)),i_TIL)%>%
      cbind(TCGA_cancer_type= rep(all_Cancers[j],nrow(i_TIL)),.)->i_TIL
    
    if(i == 1 & j==1){
      
      all_exprs <- i_expr
      all_TILs <- i_TIL
    }else{
      
      rbind(all_exprs,i_expr)->all_exprs
      rbind(all_TILs,i_TIL)->all_TILs
    }
    
  }
  
  
  
}

tidyr::nest(all_exprs,-TCGA_cancer_type,-Project_ID)->Exprs
dplyr::mutate(Exprs,expr = purrr::map2(data,Project_ID,.f=function(.x,.y){ 
  print(.y) 
  .x %>% tidyr::spread(Samples,Values)
}))%>%
  dplyr::select(-data)->Exprs  
colnames(Exprs)[3] <- "Exprs"


tidyr::nest(all_TILs,-TCGA_cancer_type,-Project_ID)->TILs
colnames(TILs)[3] <- "TILs"

dplyr::select(metadata,-Tumor_stage) %>% 
  tidyr::nest(-TCGA_cancer_type, -Project_ID)->Survivals
colnames(Survivals)[3] <- "Survivals" 
dplyr::select(metadata,TCGA_cancer_type,Project_ID,Sample_ID,Tumor_stage) %>% 
  tidyr::nest(-TCGA_cancer_type, -Project_ID)->Stages
colnames(Stages)[3] <- "Stages"


merge(Survivals,Stages, by = c("TCGA_cancer_type","Project_ID"))%>%
  merge(TILs,by = c("TCGA_cancer_type","Project_ID"))%>%
  merge(Exprs,by = c("TCGA_cancer_type","Project_ID"))->clinical_info 

clinical_info %>% 
  readr::write_rds("/home/liull/TCGA_nomograph/GEO_TNM_clinical_TIL_expr.rds.gz")

#select related genes
checkpoints <- readr::read_tsv("/home/liull/TCGA_nomograph/ICPs_all_info_class.tsv") 

genelist_exp <- readr::read_rds("/home/liull/TCGA_nomograph/GEO_TNM_clinical_TIL_expr.rds.gz") %>% 
  dplyr::mutate(exp_filter = purrr::map2(Exprs,Project_ID,.f=function(.x,.y){
    
    print(.y) 
    .x %>% 
      dplyr::filter(Symbol %in% checkpoints$symbol) 
    
  })) %>%
  dplyr::select(-Exprs)  

genelist_exp %>%
  readr::write_rds("/home/liull/TCGA_nomograph/Filtered_GEO_TNM_clinical_TIL_expr.rds.gz")

#death ratio
GEO_TNM <- readr::read_rds("/home/liull/TCGA_nomograph/Filtered_GEO_TNM_clinical_TIL_expr.rds.gz")
Death_ratio <- data.frame(cancer_types=NA,OS_death_ratio_1y=NA,OS_death_ratio_3y=NA,OS_death_ratio_5y=NA)

for (i in 1:nrow(GEO_TNM)) {
  GEO_TNM[[3]][[i]]%>%as.data.frame()%>%dplyr::select(OS_days,OS_status) -> OS_infor
  as.numeric(OS_infor$OS_days)->OS_infor$OS_days
  as.numeric(OS_infor$OS_status)->OS_infor$OS_status
  
  
  dplyr::filter(OS_infor,OS_status == 1)%>%
    dplyr::filter(OS_days < 365)%>% nrow()->OS_death_num_1y
  dplyr::filter(OS_infor,OS_status == 1)%>%
    dplyr::filter(OS_days < 1096)%>% nrow()->OS_death_num_3y
  dplyr::filter(OS_infor,OS_status == 1)%>%
    dplyr::filter(OS_days < 1825)%>% nrow()->OS_death_num_5y
  
  
  Death_ratio[i,] <- c(GEO_TNM$TCGA_cancer_type[i],signif(OS_death_num_1y/nrow(OS_infor),3),
                       signif(OS_death_num_3y/nrow(OS_infor),3),signif(OS_death_num_5y/nrow(OS_infor),3))
}

write_tsv(Death_ratio,"/home/liull/TCGA_nomograph/GEO_TNM_death_ratio_summary.tsv")
