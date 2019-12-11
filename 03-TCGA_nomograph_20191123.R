#####build nomogram model use ICP genes' ssGSEA score and TIL(pure tumor/fixed tumor)
##define high/low in all_data
##use the uni survival result
library(magrittr)
library(rms)
options(stringsAsFactors = FALSE)


split_values <- function(data1,n){
  out_data <- data1
  for (j in n:ncol(data1)) {
    
    median(data1[,j],na.rm = T) -> Values
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
  return(out_data)
}


fn_select_samples <- function(data){
  
  Num =0
  for (j in 1:1000) {
    sample(data$barcode, ceiling(nrow(data)*0.7), replace = FALSE)->train_samples
    dplyr::filter(data,barcode %in%train_samples)->train_data
    
    setdiff(data$barcode,train_samples)->test_samples
    dplyr::filter(data,barcode %in%test_samples)->test_data
    for (i in 4:ncol(test_data)) {
      if(length(unique(as.character(as.matrix(train_data[,i]))))==1){
        Num=Num+1
      }
      if(length(unique(as.character(as.matrix(test_data[,i]))))==1){
        Num=Num+1
      }
    }
    j++
    
    if(Num == 0){
      break;
    }
    
  }
  return(list(train_samples,test_samples))
  
}

fn_nomogram_model <- function(data){
  
  if(max(data$OS,na.rm = TRUE) >= 1825){
    f_cph <- cph(Surv(OS,Status) ~., data=data,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825,singular.ok = TRUE)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  data)[[1]]%>%signif(3)->C_index
    
    # surv <- Survival(f_cph)
    # nom_cph <- nomogram(f_cph, fun=list(function(x) surv(1095, x),function(x) surv(1825, x)),
    #                     fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
    #                     funlabel=c("3-Year-Survival","5-Year-Survival"),
    #                     lp=F
    # )
  }else{
    f_cph <- cph(Surv(OS,Status) ~., data=data,surv=TRUE,x=TRUE, y=TRUE,time.inc=1095,singular.ok = TRUE)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  data)[[1]] ->C_index
    
    # surv <- Survival(f_cph)
    # nom_cph <- nomogram(f_cph, fun=list(function(x) surv(365, x),function(x) surv(1095, x)),
    #                     fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
    #                     funlabel=c("3-Year-Survival","5-Year-Survival"),
    #                     lp=F
    # )
  }
  
  return(list(f_cph,C_index))
  
}

#pure TIL ---------------------------------------------------------------------------------------
#uni_res.sig
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/muli_uni_survival_20191123/OS_pure_uni_res.sig.tsv")->OS_pure_uni_sig

TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_fixed_TIL_ssGSEA_data_20191123.rds.gz")%>%
  dplyr::filter(cancer_types %in% unique(OS_pure_uni_sig$cancer_types))
index_value <- data.frame(cancer_types=NA,features=NA,
                          TNM_train_C=NA,TNM_test_C=NA,
                          gene_sets_train_C=NA,gene_sets_test_C=NA,
                          pure_TIL_train_C=NA,pure_TIL_test_C=NA,
                          TNM_gene_sets_train_C=NA,TNM_gene_sets_test_C=NA,
                          TNM_pure_TIL_train_C=NA,TNM_pure_TIL_test_C=NA,
                          gene_sets_pure_TIL_train_C=NA,gene_sets_pure_TIL_test_C=NA,
                          All_train_C=NA,All_test_C=NA)

#not consider age
for (i in 1:nrow(TCGA_all)) {
  
  dplyr::filter(OS_pure_uni_sig,cancer_types == TCGA_all$cancer_types[i])%>%
    dplyr::select(features)%>%as.matrix()%>%as.character()%>%unique()->i_all_features
  
  TCGA_all[[2]][[i]]%>%as.data.frame()%>%
    dplyr::select(-Age)%>%
    dplyr::filter(Stage != "NA")-> OS_stage
  OS_stage$OS=as.numeric(OS_stage$OS)
  OS_stage$Status=as.numeric(OS_stage$Status)
  
  #OS_stage$Age <- as.numeric(OS_stage$Age)
  gsub("^stage x$","NA",OS_stage$Stage) %>%
    gsub("^i/ii nos$","NA",.) %>% gsub("^is$","NA",.)%>% 
    gsub("^stage 0$","0",.) %>% 
    gsub("^stage ib$","I",.) %>% gsub("^stage i$","I",.) %>% gsub("^stage ia$","I",.)%>%
    gsub("^stage ii$","II",.) %>% gsub("^stage iic$","II",.) %>% gsub("^stage iib$","II",.) %>% gsub("^stage iia$","II",.) %>% 
    gsub("^stage iiic$","III",.) %>% gsub("^stage iii$","III",.) %>% gsub("^stage iiib$","III",.) %>% gsub("^stage iiia$","III",.) %>% 
    gsub("^stage iv$","IV",.) %>% gsub("^stage iva$","IV",.) %>% gsub("^stage ivb$","IV",.) %>% gsub("^stage ivc$","IV",.) ->OS_stage$Stage
  dplyr::filter(OS_stage,Stage != "NA")-> OS_stage
  # dplyr::select(OS_stage,barcode,OS,Status,Stage,intersect(i_all_features,colnames(OS_stage)))%>%
  #   dplyr::filter(Stage != "NA")->OS_stage
  
  
  TCGA_all[[4]][[i]]%>%as.data.frame()  -> ssGSEA
  dplyr::select(ssGSEA,barcode,intersect(i_all_features,colnames(ssGSEA))) ->ssGSEA
  if(ncol(ssGSEA) > 1){
    split_values(ssGSEA,2) -> ssGSEA
  }

  TCGA_all[[5]][[i]]%>%as.data.frame()  -> pure_TIL
  dplyr::select(pure_TIL,barcode,intersect(i_all_features,colnames(pure_TIL))) ->pure_TIL
  if(ncol(pure_TIL) > 1){
    split_values(pure_TIL,2)->pure_TIL
  }
  
  
  if(nrow(OS_stage) == 0){
    index_value[i,c(3:4,9:12,15:16)]<- c(rep(0,8))
    
    TCGA_all[[2]][[i]]%>%as.data.frame()%>%
      dplyr::select(-Age)-> OS_stage
    OS_stage$OS=as.numeric(OS_stage$OS)
    OS_stage$Status=as.numeric(OS_stage$Status)
    merge(OS_stage,ssGSEA)%>%
      merge(pure_TIL)%>%dplyr::select(-Stage)->all_sigs
    ddist <- datadist(all_sigs)
    options(datadist ='ddist')
    fn_select_samples(all_sigs)->.splite_res
    .splite_res[[1]]->train_samples
    .splite_res[[2]]->test_samples
    
    if(ncol(ssGSEA) > 1 && ncol(pure_TIL) > 1){
      dplyr::filter(all_sigs,barcode %in% train_samples)%>%dplyr::select(-barcode)->ssGSEA_pure_TIL_train
      dplyr::filter(all_sigs,barcode %in% test_samples)%>%dplyr::select(-barcode)->ssGSEA_pure_TIL_test
      index_value[i,13:14] <- c(fn_nomogram_model(ssGSEA_pure_TIL_train)[2],fn_nomogram_model(ssGSEA_pure_TIL_test)[2])
    }else{
      index_value[i,13:14] <- c(rep(0,2))
    }

    
  }else{
    
    merge(OS_stage,ssGSEA)%>%
      merge(pure_TIL)->all_sigs
    ddist <- datadist(all_sigs)
    options(datadist ='ddist')
    fn_select_samples(all_sigs)->.splite_res
    .splite_res[[1]]->train_samples
    .splite_res[[2]]->test_samples

    dplyr::filter(OS_stage,barcode %in% train_samples)%>%dplyr::select(-barcode)->TNM_train
    dplyr::filter(OS_stage,barcode %in% test_samples)%>%dplyr::select(-barcode)->TNM_test
    fn_nomogram_model(TNM_train)->TNM_train_res
    fn_nomogram_model(TNM_test)->TNM_test_res
    index_value[i,3:4]<- c(TNM_train_res[2],TNM_test_res[2])
    
    if(ncol(ssGSEA) > 1){
      merge(OS_stage,ssGSEA)->TNM_gene_sets
      dplyr::filter(TNM_gene_sets,barcode %in% train_samples)%>%dplyr::select(-barcode)->TNM_gene_sets_train
      dplyr::filter(TNM_gene_sets,barcode %in% test_samples)%>%dplyr::select(-barcode)->TNM_gene_sets_test
      index_value[i,9:10] <- c(fn_nomogram_model(TNM_gene_sets_train)[2],fn_nomogram_model(TNM_gene_sets_test)[2])
    }else{
      index_value[i,9:10] <- c(rep(0,2))
    }
    
    
    if(ncol(pure_TIL) > 1){
      merge(OS_stage,pure_TIL) ->TNM_pure_TIL
      dplyr::filter(TNM_pure_TIL,barcode %in% train_samples)%>%dplyr::select(-barcode)->TNM_pure_TIL_train
      dplyr::filter(TNM_pure_TIL,barcode %in% test_samples)%>%dplyr::select(-barcode)->TNM_pure_TIL_test
      index_value[i,11:12] <- c(fn_nomogram_model(TNM_pure_TIL_train)[2],fn_nomogram_model(TNM_pure_TIL_test)[2])
    }else{
      index_value[i,11:12] <- c(rep(0,2))
    }
    
    
    if(ncol(ssGSEA) > 1 && ncol(pure_TIL) > 1){
      dplyr::filter(all_sigs,barcode %in% train_samples)%>%dplyr::select(-barcode,-Stage)->gene_sets_pure_TIL_train
      dplyr::filter(all_sigs,barcode %in% test_samples)%>%dplyr::select(-barcode,-Stage)->gene_sets_pure_TIL_test
      index_value[i,13:14] <- c(fn_nomogram_model(gene_sets_pure_TIL_train)[2],fn_nomogram_model(gene_sets_pure_TIL_test)[2])
      
      dplyr::filter(all_sigs,barcode %in% train_samples)%>%dplyr::select(-barcode)->all_sigs_train
      dplyr::filter(all_sigs,barcode %in% test_samples)%>%dplyr::select(-barcode)->all_sigs_test
      index_value[i,15:16] <- c(fn_nomogram_model(all_sigs_train)[2],fn_nomogram_model(all_sigs_test)[2])
    }else{
      index_value[i,13:14] <- c(rep(0,2))
      index_value[i,15:16] <- c(rep(0,2))
    }
    
  }
  
  

  
  index_value[i,1:2] <- c(TCGA_all$cancer_types[i],paste(i_all_features,collapse = " + "))
  if(ncol(ssGSEA) > 1){
    merge(OS_stage,ssGSEA)%>%
      dplyr::select(-Stage) ->gene_sets
    dplyr::filter(gene_sets,barcode %in% train_samples)%>%dplyr::select(-barcode)->gene_sets_train
    dplyr::filter(gene_sets,barcode %in% test_samples)%>%dplyr::select(-barcode)->gene_sets_test
    
    index_value[i,5:6] <- c(fn_nomogram_model(gene_sets_train)[2],fn_nomogram_model(gene_sets_test)[2])
    
  }else{
    index_value[i,5:6] <- c(rep(0,2))
  }
  
  if(ncol(pure_TIL) > 1){
    merge(OS_stage,pure_TIL)%>%
      dplyr::select(-Stage) ->pure_TIL
    dplyr::filter(pure_TIL,barcode %in% train_samples)%>%dplyr::select(-barcode)->pure_TIL_train
    dplyr::filter(pure_TIL,barcode %in% test_samples)%>%dplyr::select(-barcode)->pure_TIL_test
    
    index_value[i,7:8] <- c(fn_nomogram_model(pure_TIL_train)[2],fn_nomogram_model(pure_TIL_test)[2])
    
  }else{
    index_value[i,7:8] <- c(rep(0,2))
  }
  
}
for (j in 3:ncol(index_value)) {
  signif(index_value[,j],3)->index_value[,j]
}

index_value %>% readr::write_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/7_model_C_index_pure.tsv")

#compare C-index plot 
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/7_model_C_index_pure.tsv")%>%
  dplyr::filter(TNM_train_C != 0)%>%
  dplyr::filter(gene_sets_train_C != 0)%>%
  dplyr::filter(pure_TIL_train_C != 0)%>%
  dplyr::select(-features)->index_value
single_feature_train_C_plot <- rbind(data.frame(type=rep("TNM_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_train_C),
                      data.frame(type=rep("gene_sets_train_C",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_train_C))%>%
  rbind(data.frame(type=rep("pure_TIL_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$pure_TIL_train_C))
single_feature_test_C_plot <- rbind(data.frame(type=rep("TNM_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_test_C),
                                     data.frame(type=rep("gene_sets_test_C",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_test_C))%>%
  rbind(data.frame(type=rep("pure_TIL_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$pure_TIL_test_C))

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/single_feature_train_C_pure.pdf",width = 12,height = 6)
ggplot(single_feature_train_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/single_feature_test_C_pure.pdf",width = 12,height = 6)
ggplot(single_feature_test_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()


multi_feature_train_C_plot <- rbind(data.frame(type=rep("TNM_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_train_C),
                                    data.frame(type=rep("TNM_gene_sets_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_gene_sets_train_C))%>%
  rbind(data.frame(type=rep("TNM_pure_TIL_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_pure_TIL_train_C))%>%
  rbind(data.frame(type=rep("gene_sets_pure_TIL_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_pure_TIL_train_C))%>%
  rbind(data.frame(type=rep("All_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$All_train_C))
multi_feature_test_C_plot <- rbind(data.frame(type=rep("TNM_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_test_C),
                                   data.frame(type=rep("TNM_gene_sets_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_gene_sets_test_C))%>%
  rbind(data.frame(type=rep("TNM_pure_TIL_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_pure_TIL_test_C))%>%
  rbind(data.frame(type=rep("gene_sets_pure_TIL_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_pure_TIL_test_C))%>%
  rbind(data.frame(type=rep("All_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$All_test_C))

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/multi_feature_train_C_pure.pdf",width = 12,height = 6)
ggplot(multi_feature_train_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/multi_feature_test_C_pure.pdf",width = 12,height = 6)
ggplot(multi_feature_test_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()



#fixed_TIL-----------------------------------------------------------
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/muli_uni_survival_20191123/OS_fixed_uni_res.sig.tsv")->OS_fixed_uni_sig

TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_fixed_TIL_ssGSEA_data_20191123.rds.gz")%>%
  dplyr::filter(cancer_types %in% unique(OS_fixed_uni_sig$cancer_types))
filtered_cancer <- character()
for (i in 1:nrow(TCGA_all)) {
  if(nrow(TCGA_all[[6]][[i]]) > 30){
    filtered_cancer <- c(filtered_cancer,TCGA_all$cancer_types[i])
  }
}
dplyr::filter(TCGA_all,cancer_types %in% filtered_cancer)->filtered_TCGA

index_value <- data.frame(cancer_types=NA,features=NA,
                          TNM_train_C=NA,TNM_test_C=NA,
                          gene_sets_train_C=NA,gene_sets_test_C=NA,
                          fixed_TIL_train_C=NA,fixed_TIL_test_C=NA,
                          TNM_gene_sets_train_C=NA,TNM_gene_sets_test_C=NA,
                          TNM_fixed_TIL_train_C=NA,TNM_fixed_TIL_test_C=NA,
                          gene_sets_fixed_TIL_train_C=NA,gene_sets_fixed_TIL_test_C=NA,
                          All_train_C=NA,All_test_C=NA)

for (i in 1:nrow(filtered_TCGA)) {
  
  dplyr::filter(OS_fixed_uni_sig,cancer_types == filtered_TCGA$cancer_types[i])%>%
    dplyr::select(features)%>%as.matrix()%>%as.character()%>%unique()->i_all_features
  
  filtered_TCGA[[2]][[i]]%>%as.data.frame()%>%
    dplyr::select(-Age)%>%
    dplyr::filter(Stage != "NA")-> OS_stage
  OS_stage$OS=as.numeric(OS_stage$OS)
  OS_stage$Status=as.numeric(OS_stage$Status)
  
  #OS_stage$Age <- as.numeric(OS_stage$Age)
  gsub("^stage x$","NA",OS_stage$Stage) %>%
    gsub("^i/ii nos$","NA",.) %>% gsub("^is$","NA",.)%>% 
    gsub("^stage 0$","0",.) %>% 
    gsub("^stage ib$","I",.) %>% gsub("^stage i$","I",.) %>% gsub("^stage ia$","I",.)%>%
    gsub("^stage ii$","II",.) %>% gsub("^stage iic$","II",.) %>% gsub("^stage iib$","II",.) %>% gsub("^stage iia$","II",.) %>% 
    gsub("^stage iiic$","III",.) %>% gsub("^stage iii$","III",.) %>% gsub("^stage iiib$","III",.) %>% gsub("^stage iiia$","III",.) %>% 
    gsub("^stage iv$","IV",.) %>% gsub("^stage iva$","IV",.) %>% gsub("^stage ivb$","IV",.) %>% gsub("^stage ivc$","IV",.) ->OS_stage$Stage
  dplyr::filter(OS_stage,Stage != "NA")-> OS_stage
  # dplyr::select(OS_stage,barcode,OS,Status,Stage,intersect(i_all_features,colnames(OS_stage)))%>%
  #   dplyr::filter(Stage != "NA")->OS_stage
  
  
  filtered_TCGA[[4]][[i]]%>%as.data.frame()  -> ssGSEA
  dplyr::select(ssGSEA,barcode,intersect(i_all_features,colnames(ssGSEA))) ->ssGSEA
  if(ncol(ssGSEA) > 1){
    split_values(ssGSEA,2) -> ssGSEA
  }
  
  filtered_TCGA[[6]][[i]]%>%as.data.frame()  -> fixed_TIL
  dplyr::select(fixed_TIL,barcode,intersect(i_all_features,colnames(fixed_TIL))) ->fixed_TIL
  
  if(ncol(fixed_TIL) > 1){
    split_values(fixed_TIL,2)->fixed_TIL
  }
  
  
  if(nrow(OS_stage) == 0){
    index_value[i,c(3:4,9:14)]<- c(rep(0,8))
    
    filtered_TCGA[[2]][[i]]%>%as.data.frame()%>%
      dplyr::select(-Age)-> OS_stage
    OS_stage$OS=as.numeric(OS_stage$OS)
    OS_stage$Status=as.numeric(OS_stage$Status)
    merge(OS_stage,ssGSEA)%>%
      merge(fixed_TIL)%>%dplyr::select(-Stage)->all_sigs
    ddist <- datadist(all_sigs)
    options(datadist ='ddist')
    fn_select_samples(all_sigs)->.splite_res
    .splite_res[[1]]->train_samples
    .splite_res[[2]]->test_samples
    
    if(ncol(ssGSEA) > 1 && ncol(fixed_TIL) > 1){
      dplyr::filter(all_sigs,barcode %in% train_samples)%>%dplyr::select(-barcode)->ssGSEA_fixed_TIL_train
      dplyr::filter(all_sigs,barcode %in% test_samples)%>%dplyr::select(-barcode)->ssGSEA_fixed_TIL_test
      index_value[i,13:14] <- c(fn_nomogram_model(ssGSEA_fixed_TIL_train)[2],fn_nomogram_model(ssGSEA_fixed_TIL_test)[2])
    }else{
      index_value[i,13:14] <- c(rep(0,2))
    }
    
  }else{
    
    merge(OS_stage,ssGSEA)%>%
      merge(fixed_TIL)->all_sigs
    ddist <- datadist(all_sigs)
    options(datadist ='ddist')
    fn_select_samples(all_sigs)->.splite_res
    .splite_res[[1]]->train_samples
    .splite_res[[2]]->test_samples
    
    dplyr::filter(OS_stage,barcode %in% train_samples)%>%dplyr::select(-barcode)->TNM_train
    dplyr::filter(OS_stage,barcode %in% test_samples)%>%dplyr::select(-barcode)->TNM_test
    index_value[i,3:4]<- c(fn_nomogram_model(TNM_train)[2],fn_nomogram_model(TNM_test)[2])
    
    if(ncol(ssGSEA) > 1){
      merge(OS_stage,ssGSEA)->TNM_gene_sets
      dplyr::filter(TNM_gene_sets,barcode %in% train_samples)%>%dplyr::select(-barcode)->TNM_gene_sets_train
      dplyr::filter(TNM_gene_sets,barcode %in% test_samples)%>%dplyr::select(-barcode)->TNM_gene_sets_test
      index_value[i,9:10] <- c(fn_nomogram_model(TNM_gene_sets_train)[2],fn_nomogram_model(TNM_gene_sets_test)[2])
    }else{
      index_value[i,9:10] <- c(rep(0,2))
    }
    
    if(ncol(fixed_TIL) > 1){
      merge(OS_stage,fixed_TIL) ->TNM_fixed_TIL
      dplyr::filter(TNM_fixed_TIL,barcode %in% train_samples)%>%dplyr::select(-barcode)->TNM_fixed_TIL_train
      dplyr::filter(TNM_fixed_TIL,barcode %in% test_samples)%>%dplyr::select(-barcode)->TNM_fixed_TIL_test
      index_value[i,11:12] <- c(fn_nomogram_model(TNM_fixed_TIL_train)[2],fn_nomogram_model(TNM_fixed_TIL_test)[2])
    }else{
      index_value[i,11:12] <- c(rep(0,2))
    }
    
    
    
    if(ncol(ssGSEA) > 1 && ncol(fixed_TIL) > 1){
      dplyr::filter(all_sigs,barcode %in% train_samples)%>%dplyr::select(-barcode,-Stage)->gene_sets_fixed_TIL_train
      dplyr::filter(all_sigs,barcode %in% test_samples)%>%dplyr::select(-barcode,-Stage)->gene_sets_fixed_TIL_test
      index_value[i,13:14] <- c(fn_nomogram_model(gene_sets_fixed_TIL_train)[2],fn_nomogram_model(gene_sets_fixed_TIL_test)[2])
      
      dplyr::filter(all_sigs,barcode %in% train_samples)%>%dplyr::select(-barcode)->all_sigs_train
      dplyr::filter(all_sigs,barcode %in% test_samples)%>%dplyr::select(-barcode)->all_sigs_test
      index_value[i,15:16] <- c(fn_nomogram_model(all_sigs_train)[2],fn_nomogram_model(all_sigs_test)[2])
    }else{
      index_value[i,13:14] <- c(rep(0,2))
      index_value[i,15:16] <- c(rep(0,2))
    }
    
  }
  
  
  index_value[i,1:2] <- c(filtered_TCGA$cancer_types[i],paste(i_all_features,collapse = " + "))
  if(ncol(ssGSEA) > 1){
    merge(OS_stage,ssGSEA)%>%
      dplyr::select(-Stage) ->gene_sets
    dplyr::filter(gene_sets,barcode %in% train_samples)%>%dplyr::select(-barcode)->gene_sets_train
    dplyr::filter(gene_sets,barcode %in% test_samples)%>%dplyr::select(-barcode)->gene_sets_test
    index_value[i,5:6] <- c(fn_nomogram_model(gene_sets_train)[2],fn_nomogram_model(gene_sets_test)[2])
    
  }else{
    index_value[i,5:6] <- c(rep(0,2))
  }
  
  if(ncol(fixed_TIL) > 1){
    merge(OS_stage,fixed_TIL)%>%
      dplyr::select(-Stage) ->fixed_TIL
    dplyr::filter(fixed_TIL,barcode %in% train_samples)%>%dplyr::select(-barcode)->fixed_TIL_train
    dplyr::filter(fixed_TIL,barcode %in% test_samples)%>%dplyr::select(-barcode)->fixed_TIL_test
    
    index_value[i,7:8] <- c(fn_nomogram_model(fixed_TIL_train)[2],fn_nomogram_model(fixed_TIL_test)[2])
    
  }else{
    index_value[i,7:8] <- c(rep(0,2))
  }
  
}
for (j in 3:ncol(index_value)) {
  signif(index_value[,j],3)->index_value[,j]
}
index_value %>% readr::write_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/7_model_C_index_fixed.tsv")


#compare C-index plot 
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/7_model_C_index_fixed.tsv")%>%
  dplyr::filter(TNM_train_C != 0)%>%
  dplyr::filter(gene_sets_train_C != 0)%>%
  dplyr::filter(fixed_TIL_train_C != 0)%>%
  dplyr::select(-features)->index_value
single_feature_train_C_plot <- rbind(data.frame(type=rep("TNM_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_train_C),
                                     data.frame(type=rep("gene_sets_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_train_C))%>%
  rbind(data.frame(type=rep("fixed_TIL_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$fixed_TIL_train_C))
single_feature_test_C_plot <- rbind(data.frame(type=rep("TNM_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_test_C),
                                    data.frame(type=rep("gene_sets_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_test_C))%>%
  rbind(data.frame(type=rep("fixed_TIL_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$fixed_TIL_test_C))

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/single_feature_train_C_fixed.pdf",width = 12,height = 6)
ggplot(single_feature_train_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/single_feature_test_C_fixed.pdf",width = 12,height = 6)
ggplot(single_feature_test_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()


multi_feature_train_C_plot <- rbind(data.frame(type=rep("TNM_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_train_C),
                                    data.frame(type=rep("TNM_gene_sets_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_gene_sets_train_C))%>%
  rbind(data.frame(type=rep("TNM_fixed_TIL_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_fixed_TIL_train_C))%>%
  rbind(data.frame(type=rep("gene_sets_fixed_TIL_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_fixed_TIL_train_C))%>%
  rbind(data.frame(type=rep("All_train",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$All_train_C))
multi_feature_test_C_plot <- rbind(data.frame(type=rep("TNM_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_test_C),
                                   data.frame(type=rep("TNM_gene_sets_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_gene_sets_test_C))%>%
  rbind(data.frame(type=rep("TNM_fixed_TIL_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$TNM_fixed_TIL_test_C))%>%
  rbind(data.frame(type=rep("gene_sets_fixed_TIL_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$gene_sets_fixed_TIL_test_C))%>%
  rbind(data.frame(type=rep("All_test",nrow(index_value)),cancer_types=index_value$cancer_types,C_index=index_value$All_test_C))

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/multi_feature_train_C_fixed.pdf",width = 12,height = 6)
ggplot(multi_feature_train_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()

pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_ssGSEA/multi_feature_test_C_fixed.pdf",width = 12,height = 6)
ggplot(multi_feature_test_C_plot, aes(x=cancer_types,y=C_index,colour=type,group=type))+ geom_line() + theme_bw()
dev.off()

