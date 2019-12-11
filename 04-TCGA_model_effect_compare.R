#######compare model function effect using the best features
##do not define high/low,use the number to do the survival analysis
##use the uni survival result
##not consider Age now

rm(list=ls())
library(magrittr)
options(stringsAsFactors = FALSE)

fn_select_samples <- function(data){
  
  
  dplyr::filter(data,Status==1)%>%dplyr::select(barcode)%>%as.matrix()%>%as.character->dead_samples
  setdiff(data$barcode,dead_samples)->alive_samples
  
  sample(dead_samples, floor(length(dead_samples)*0.7), replace = FALSE)->train_dead_samples
  setdiff(dead_samples,train_dead_samples)->test_dead_samples
  
  sample(alive_samples, floor(length(alive_samples)*0.7), replace = FALSE)->train_alive_samples
  setdiff(alive_samples,train_alive_samples)->test_alive_samples
  
  
  return(list(c(train_alive_samples,train_dead_samples),c(test_alive_samples,test_dead_samples)))
  
}

fn_nomogram_model <- function(data){
  
  if(max(data$OS,na.rm = TRUE) >= 1825){
    f_cph <- cph(Surv(OS,Status) ~., data=data,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825,singular.ok = TRUE)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  data)[[1]]%>%signif(3)->C_index
    
    
    # validate(f_cph, method="boot", B=10, dxy=T)->b
    # b[1,5]/2+0.5->C_index
    
    
    surv <- Survival(f_cph)
    nom_cph <- nomogram(f_cph, fun=list(function(x) surv(1095, x),function(x) surv(1825, x)),
                        fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                        funlabel=c("3-Year-Survival","5-Year-Survival"),
                        lp=F
    )
    
  }else{
    f_cph <- cph(Surv(OS,Status) ~., data=data,surv=TRUE,x=TRUE, y=TRUE,time.inc=1095,singular.ok = TRUE)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  data)[[1]] ->C_index
    
  }
  
  return(list(f_cph,C_index))
}

readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model_20191204/TNM_C_index.tsv")%>%
  dplyr::select(cancer_types,TNM_C_index)%>%
  dplyr::filter(!is.na(TNM_C_index))->TNM_index
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model_20191204/ssGSEA_C_index.tsv")%>%
  dplyr::select(cancer_types,new_gene_set_features)%>%
  dplyr::filter(!is.na(new_gene_set_features))->ssGSEA_feature
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model_20191204/pure_TIL_C_index.tsv")%>%
  dplyr::select(cancer_types,new_TIL_features)%>%
  dplyr::filter(!is.na(new_TIL_features))->pure_TIL_feature
merge(TNM_index,ssGSEA_feature)%>%merge(pure_TIL_feature)->all_features

TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_fixed_TIL_ssGSEA_data_20191204.rds.gz")%>%
  dplyr::filter(cancer_types %in% unique(all_features$cancer_types))

index_value <- data.frame(cancer_types=NA,
                          TNM_train_C=NA,TNM_test_C=NA,
                          gene_sets_train_C=NA,gene_sets_test_C=NA,
                          pure_TIL_train_C=NA,pure_TIL_test_C=NA,
                          TNM_gene_sets_train_C=NA,TNM_gene_sets_test_C=NA,
                          TNM_pure_TIL_train_C=NA,TNM_pure_TIL_test_C=NA,
                          gene_sets_pure_TIL_train_C=NA,gene_sets_pure_TIL_test_C=NA,
                          All_train_C=NA,All_test_C=NA)

for (i in 1:nrow(TCGA_all)) {
  
  TCGA_all$OS_stage[[i]]%>%as.data.frame()%>%
    dplyr::select(-Age)-> OS_stage
  OS_stage$OS=as.numeric(OS_stage$OS)
  OS_stage$Status=as.numeric(OS_stage$Status)
  
  gsub("^stage x$","NA",OS_stage$Stage) %>%
    gsub("^i/ii nos$","NA",.) %>% gsub("^is$","NA",.)%>% 
    gsub("^stage 0$","0",.) %>% 
    gsub("^stage ib$","I",.) %>% gsub("^stage i$","I",.) %>% gsub("^stage ia$","I",.)%>%
    gsub("^stage ii$","II",.) %>% gsub("^stage iic$","II",.) %>% gsub("^stage iib$","II",.) %>% gsub("^stage iia$","II",.) %>% 
    gsub("^stage iiic$","III",.) %>% gsub("^stage iii$","III",.) %>% gsub("^stage iiib$","III",.) %>% gsub("^stage iiia$","III",.) %>% 
    gsub("^stage iv$","IV",.) %>% gsub("^stage iva$","IV",.) %>% gsub("^stage ivb$","IV",.) %>% gsub("^stage ivc$","IV",.) ->OS_stage$Stage
  dplyr::filter(OS_stage,Stage != "NA")->OS_stage
  
  strsplit(all_features$new_gene_set_features[i]," + ",fixed = TRUE)%>%unlist()->gene_set_features_i
  TCGA_all$ssGSEA[[i]]%>%as.data.frame()%>%
    dplyr::select(barcode,gene_set_features_i)->ssGSEA
  
  strsplit(all_features$new_TIL_features[i]," + ",fixed = TRUE)%>%unlist()->TIL_features_i
  TCGA_all$pure_tumor_Infiltration[[i]]%>%as.data.frame()%>%
    dplyr::select(barcode,TIL_features_i)->TIL
  
  #TNM single-----------------------------------------------------------------------------------------
  dplyr::select(OS_stage,-barcode)->TNM_1
  ddist <- datadist(TNM_1)
  options(datadist ='ddist')
  
  TNM_train_C <- 0
  TNM_test_C <- 0
  for (j in 1:1000) {
    fn_select_samples(OS_stage)->.splite_samples
    dplyr::filter(OS_stage,barcode %in% .splite_samples[[1]])->TNM_train
    dplyr::filter(OS_stage,barcode %in% .splite_samples[[2]])%>%
      dplyr::select(-barcode)->TNM_test
    dplyr::select(TNM_train,-barcode)->TNM_train
    
    TNM_train_C + fn_nomogram_model(TNM_train)[[2]]->TNM_train_C
    TNM_test_C + fn_nomogram_model(TNM_test)[[2]]->TNM_test_C
    j=j+1
  }
  
  #gene set single---------------------------------------------------------------------------
  dplyr::select(OS_stage,-Stage)%>%merge(ssGSEA)->ssGSEA_1
  ddist <- datadist(ssGSEA_1)
  options(datadist ='ddist')
  
  ssGSEA_train_C <- 0
  ssGSEA_test_C <- 0
  for (j in 1:1000) {
    fn_select_samples(ssGSEA_1)->.splite_samples
    dplyr::filter(ssGSEA_1,barcode %in% .splite_samples[[1]])->ssGSEA_train
    dplyr::filter(ssGSEA_1,barcode %in% .splite_samples[[2]])%>%
      dplyr::select(-barcode)->ssGSEA_test
    dplyr::select(ssGSEA_train,-barcode)->ssGSEA_train
    
    ssGSEA_train_C + fn_nomogram_model(ssGSEA_train)[[2]]->ssGSEA_train_C
    ssGSEA_test_C + fn_nomogram_model(ssGSEA_test)[[2]]->ssGSEA_test_C
    j=j+1
  }
  
  #TIL single--------------------------------------------------------------------------------------------
  dplyr::select(OS_stage,-Stage)%>%merge(TIL)->TIL_1
  ddist <- datadist(TIL_1)
  options(datadist ='ddist')
  
  TIL_train_C <- 0
  TIL_test_C <- 0
  for (j in 1:1000) {
    fn_select_samples(TIL_1)->.splite_samples
    dplyr::filter(TIL_1,barcode %in% .splite_samples[[1]])->TIL_train
    dplyr::filter(TIL_1,barcode %in% .splite_samples[[2]])%>%
      dplyr::select(-barcode)->TIL_test
    dplyr::select(TIL_train,-barcode)->TIL_train
    
    TIL_train_C + fn_nomogram_model(TIL_train)[[2]]->TIL_train_C
    TIL_test_C + fn_nomogram_model(TIL_test)[[2]]->TIL_test_C
    j=j+1
  }
  
  #TNM + gene set------------------------------------------------------------------------
  merge(OS_stage,ssGSEA)->TNM_ssGSEA
  ddist <- datadist(TNM_ssGSEA)
  options(datadist ='ddist')
  
  TNM_ssGSEA_train_C <- 0
  TNM_ssGSEA_test_C <- 0
  for (j in 1:1000) {
    fn_select_samples(TNM_ssGSEA)->.splite_samples
    dplyr::filter(TNM_ssGSEA,barcode %in% .splite_samples[[1]])->TNM_ssGSEA_train
    dplyr::filter(TNM_ssGSEA,barcode %in% .splite_samples[[2]])%>%
      dplyr::select(-barcode)->TNM_ssGSEA_test
    dplyr::select(TNM_ssGSEA_train,-barcode)->TNM_ssGSEA_train
    
    TNM_ssGSEA_train_C + fn_nomogram_model(TNM_ssGSEA_train)[[2]]->TNM_ssGSEA_train_C
    TNM_ssGSEA_test_C + fn_nomogram_model(TNM_ssGSEA_test)[[2]]->TNM_ssGSEA_test_C
    j=j+1
  }
  
  #TNM + TIL------------------------------------------------------------------------
  merge(OS_stage,TIL)->TNM_TIL
  ddist <- datadist(TNM_TIL)
  options(datadist ='ddist')
  
  TNM_TIL_train_C <- 0
  TNM_TIL_test_C <- 0
  for (j in 1:1000) {
    fn_select_samples(TNM_TIL)->.splite_samples
    dplyr::filter(TNM_TIL,barcode %in% .splite_samples[[1]])->TNM_TIL_train
    dplyr::filter(TNM_TIL,barcode %in% .splite_samples[[2]])%>%
      dplyr::select(-barcode)->TNM_TIL_test
    dplyr::select(TNM_TIL_train,-barcode)->TNM_TIL_train
    
    TNM_TIL_train_C + fn_nomogram_model(TNM_TIL_train)[[2]]->TNM_TIL_train_C
    TNM_TIL_test_C + fn_nomogram_model(TNM_TIL_test)[[2]]->TNM_TIL_test_C
    j=j+1
  }
  
  #gene set + TIL------------------------------------------------------------------------
  dplyr::select(OS_stage,-Stage)%>%
    merge(ssGSEA)%>%merge(TIL)->ssGSEA_TIL
  ddist <- datadist(ssGSEA_TIL)
  options(datadist ='ddist')
  
  ssGSEA_TIL_train_C <- 0
  ssGSEA_TIL_test_C <- 0
  for (j in 1:1000) {
    fn_select_samples(ssGSEA_TIL)->.splite_samples
    dplyr::filter(ssGSEA_TIL,barcode %in% .splite_samples[[1]])->ssGSEA_TIL_train
    dplyr::filter(ssGSEA_TIL,barcode %in% .splite_samples[[2]])%>%
      dplyr::select(-barcode)->ssGSEA_TIL_test
    dplyr::select(ssGSEA_TIL_train,-barcode)->ssGSEA_TIL_train
    
    ssGSEA_TIL_train_C + fn_nomogram_model(ssGSEA_TIL_train)[[2]]->ssGSEA_TIL_train_C
    ssGSEA_TIL_test_C + fn_nomogram_model(ssGSEA_TIL_test)[[2]]->ssGSEA_TIL_test_C
    j=j+1
  }
  
  #TNM + gene set +TIL------------------------------------------------------------------------
  merge(OS_stage,ssGSEA)%>%merge(TIL)->all_3
  ddist <- datadist(all_3)
  options(datadist ='ddist')
  
  all_train_C <- 0
  all_test_C <- 0
  for (j in 1:1000) {
    fn_select_samples(all_3)->.splite_samples
    dplyr::filter(all_3,barcode %in% .splite_samples[[1]])->all_train
    dplyr::filter(all_3,barcode %in% .splite_samples[[2]])%>%
      dplyr::select(-barcode)->all_test
    dplyr::select(all_train,-barcode)->all_train
    
    all_train_C + fn_nomogram_model(all_train)[[2]]->all_train_C
    all_test_C + fn_nomogram_model(all_test)[[2]]->all_test_C
    j=j+1
  }
  
  index_value[i,] <- c(TCGA_all$cancer_types[i],TNM_train_C/1000,TNM_test_C/1000,
                       ssGSEA_train_C/1000,ssGSEA_test_C/1000,
                       TIL_train_C/1000,TIL_test_C/1000,
                       TNM_ssGSEA_train_C /1000,TNM_ssGSEA_test_C/1000,
                       TNM_TIL_train_C/1000,TNM_TIL_test_C/1000,
                       ssGSEA_TIL_train_C/1000,ssGSEA_TIL_test_C/1000,
                       all_train_C/1000,all_test_C/1000)
}

