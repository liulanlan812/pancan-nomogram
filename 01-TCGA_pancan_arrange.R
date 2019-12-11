# #TCGA's rds file,remove normal samples.
# options(stringsAsFactors = FALSE)
# 
# all_exp <- readr::read_rds("/home/liull/TCGA_nomograph/pancan33_expr.rds.gz")
# 
# all_exp %>% 
#   dplyr::mutate(cancer_expr = purrr::map2(expr,cancer_types,.f=function(.x,.y){
#     
#     print(.y)
#     .x %>% 
#       as.data.frame()%>%
#       dplyr::filter(symbol != "?") ->.x
#     .x[which(.x$entrez_id == "728661"),1] <- "SLC35E2B"
#     .x[which(.x$entrez_id == "9906"),1] <- "SLC35E2A"
#      rownames(.x) <- .x$symbol
#      
#      which(substr(colnames(.x),14,14) == "0")->cancer_id
#      
#     .x %>%
#       dplyr::select(cancer_id) 
#     
#   })) %>%
#   dplyr::select(-expr) ->Cancer_expr
# Cancer_expr %>%
#   readr::write_rds("/home/liull/TCGA_nomograph/TCGA_result/pancan33_cancer_expr.rds.gz")
#ImmuCellAI.R


#load("/home/liull/TCGA_nomograph/TCGA_result/TCGA_infil_ls.Rdata")#all TCGA's cancer sample's TIL


#summary all TCGA-pancancer's survival situation------------------------------------------------------------
TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_fixed_TIL_ssGSEA_data_20191123.rds.gz")

Death_ratio <- data.frame(cancer_types=NA,cancer_Num=NA,paired_normal_Num=NA,OS_death_ratio_1y=NA,OS_death_ratio_3y=NA,OS_death_ratio_5y=NA,
                          PFS_death_ratio_1y=NA,PFS_death_ratio_3y=NA,PFS_death_ratio_5y=NA)

for (i in 1:nrow(TCGA_all)) {
  
  
  TCGA_all[[2]][[i]]%>%as.data.frame()%>%dplyr::select(OS,Status) -> OS_infor
  as.numeric(OS_infor$OS)->OS_infor$OS
  as.numeric(OS_infor$Status)->OS_infor$Status
  
  TCGA_all[[3]][[i]]%>%as.data.frame()%>%dplyr::select(PFS.time,PFS) -> PFS_infor
  as.numeric(PFS_infor$PFS.time)->PFS_infor$PFS.time
  as.numeric(PFS_infor$PFS)->PFS_infor$PFS
  
  dplyr::filter(OS_infor,Status == 1)%>%
    dplyr::filter(OS < 365)%>% nrow()->OS_death_num_1y
  dplyr::filter(OS_infor,Status == 1)%>%
    dplyr::filter(OS < 1096)%>% nrow()->OS_death_num_3y
  dplyr::filter(OS_infor,Status == 1)%>%
    dplyr::filter(OS < 1825)%>% nrow()->OS_death_num_5y
  
  dplyr::filter(PFS_infor,PFS == 1)%>%
    dplyr::filter(PFS.time < 365)%>% nrow()->PFS_death_num_1y
  dplyr::filter(PFS_infor,PFS == 1)%>%
    dplyr::filter(PFS.time < 1096)%>% nrow()->PFS_death_num_3y
  dplyr::filter(PFS_infor,PFS == 1)%>%
    dplyr::filter(PFS.time < 1825)%>% nrow()->PFS_death_num_5y
  if(nrow(TCGA_all[[4]][[i]]) == nrow(TCGA_all[[5]][[i]])){
    normal_num <- 0
  }else{
    normal_num <- nrow(TCGA_all[[6]][[i]])
  }
  
  Death_ratio[i,] <- c(TCGA_all$cancer_types[i],nrow(TCGA_all[[4]][[i]]),normal_num,
                       signif(OS_death_num_1y/nrow(OS_infor),3),signif(OS_death_num_3y/nrow(OS_infor),3),signif(OS_death_num_5y/nrow(OS_infor),3),
                       signif(PFS_death_num_1y/nrow(PFS_infor),3),signif(PFS_death_num_3y/nrow(PFS_infor),3),signif(PFS_death_num_5y/nrow(PFS_infor),3))
}

write_tsv(Death_ratio,"/home/liull/TCGA_nomograph/TCGA_death_ratio_summary.tsv")

