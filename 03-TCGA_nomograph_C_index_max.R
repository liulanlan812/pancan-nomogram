
library(magrittr)
library(rms)
options(stringsAsFactors = FALSE)


split_values <- function(data1){
  
  out_data <- data1
  median_values <- data.frame(cancer_types=NA,features=NA,medians=NA)
  
  for (j in 2:ncol(data1)) {
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


readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/muli_uni_survival/OS_multi_res.tsv")->OS_multi_res

TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_TIL_exp_data.rds.gz")%>%
  dplyr::filter(cancer_types %in% unique(OS_multi_res$cancer_types))


#find the combination to max the C-index(1-5 features)----------------------------------------------
OS_C_index_value <- data.frame(cancer_types=NA,features=NA,Max_OS_C_index=NA)

for (i in 1:(nrow(TCGA_all))) {
  
  #one cancer's significant features
  dplyr::filter(OS_multi_res,cancer_types == TCGA_all$cancer_types[i])%>%
    dplyr::select(Features)%>%as.matrix()%>%as.character()->i_all_features
  gsub("_high","",i_all_features)-> i_all_features
  
  
  #stage,survival,expr,TIL =>all_sigs
  TCGA_all[[4]][[i]]%>%as.data.frame() -> OS_stage
  OS_stage$Age <- as.numeric(OS_stage$Age)
  gsub("^stage x$","NA",OS_stage$Stage) %>%
    gsub("^i/ii nos$","NA",.) %>% gsub("^is$","NA",.)%>% 
    gsub("^stage 0$","0",.) %>% 
    gsub("^stage ib$","I",.) %>% gsub("^stage i$","I",.) %>% gsub("^stage ia$","I",.)%>%
    gsub("^stage ii$","II",.) %>% gsub("^stage iic$","II",.) %>% gsub("^stage iib$","II",.) %>% gsub("^stage iia$","II",.) %>% 
    gsub("^stage iiic$","III",.) %>% gsub("^stage iii$","III",.) %>% gsub("^stage iiib$","III",.) %>% gsub("^stage iiia$","III",.) %>% 
    gsub("^stage iv$","IV",.) %>% gsub("^stage iva$","IV",.) %>% gsub("^stage ivb$","IV",.) %>% gsub("^stage ivc$","IV",.) ->OS_stage$Stage
  dplyr::select(OS_stage,barcode,OS,Status,Stage,intersect(i_all_features,colnames(OS_stage)))%>%
    dplyr::filter(Stage != "NA")->OS_stage
  
  
  TCGA_all[[2]][[i]]%>%as.data.frame()  -> TIL
  dplyr::select(TIL,barcode,intersect(i_all_features,colnames(TIL))) ->TIL
  if(ncol(TIL) > 1){
    split_values(TIL) ->split_res
    split_res[[1]]->TIL
    split_res[[2]]->i_median_values
  }
  
  
  TCGA_all[[3]][[i]]%>%as.data.frame() -> Expr
  rownames(Expr) <- Expr$symbol
  Expr[,-1] ->Expr
  gsub("-", ".", rownames(Expr),fixed = TRUE) ->rownames(Expr)
  t(Expr)%>%as.data.frame()%>%
    dplyr::select(intersect(i_all_features,rownames(Expr)))%>%
    tibble::rownames_to_column()->Expr
  colnames(Expr)[1] <- "rowname"
  
  if(ncol(Expr) > 1){
    Samples_ID <- Expr$rowname
    Features <- colnames(Expr)
    Expr[,-1]%>%as.data.frame() -> Expr
    log2(Expr + 0.01)->Expr
    cbind(rowname = Samples_ID,Expr)->Expr
    colnames(Expr) <- Features
    
    split_values(Expr)[[1]] -> Expr
  }
  
  merge(OS_stage,TIL) %>%
    merge(Expr,by.x = "barcode",by.y = "rowname")->all_sigs
  dplyr::select(all_sigs,-barcode)->all_sigs
  all_sigs$OS=as.numeric(all_sigs$OS)
  all_sigs$Status=as.numeric(all_sigs$Status)
  
  
  # small_sigs <- character()
  
  
  if(ncol(all_sigs) >=8){
    Features_limit <- 5
  }else{
    Features_limit <- ncol(all_sigs)-3
  }
  All_C_index <-list()
  for (n in 1:Features_limit) {
    print(n)
    C_index <- character()
    
    combn(4:ncol(all_sigs),n)->all_3_zuhe
    
    for (m in 1:ncol(all_3_zuhe)) {
      selected_features <- colnames(all_sigs)[all_3_zuhe[,m]]
      all_sigs %>% dplyr::select(OS,Status,Stage,selected_features)->Filtered_sigs
      
      ddist <- datadist(Filtered_sigs)
      options(datadist='ddist')
      
      f_cph <- cph(Surv(OS,Status) ~., data=Filtered_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825)
      1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  Filtered_sigs)[[1]]->m_C_index
      
      c(C_index,m_C_index)->C_index
      
    }
    All_C_index[[n]] <- C_index
    
  }
  Max_C_index <- character()
  for (p in 1:length(All_C_index)) {
    c(Max_C_index,max(All_C_index[[p]]))->Max_C_index
  }
  n <- which(Max_C_index[] == max(Max_C_index))
  m <- which(All_C_index[[n]][] == max(All_C_index[[n]][]))
  selected_features <- colnames(all_sigs)[combn(4:ncol(all_sigs),n)[,m]]
  OS_C_index_value[i,] <- c(TCGA_all$cancer_types[i],paste(colnames(all_sigs)[combn(4:ncol(all_sigs),n)[,m]],collapse = " + "),
                            max(Max_C_index))
  
  all_sigs %>% dplyr::select(OS,Status,Stage,selected_features)->Filtered_sigs
  
  
  ddist <- datadist(Filtered_sigs)
  options(datadist='ddist')
  
  f_cph <- cph(Surv(OS,Status) ~., data=Filtered_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825)
  1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  Filtered_sigs)[[1]]->i_C_index
  
  surv <- Survival(f_cph)
  nom_cph <- nomogram(f_cph, fun=list(function(x) surv(1095, x),function(x) surv(1825, x)),
                      fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                      funlabel=c("3-Year-Survival","5-Year-Survival"),
                      lp=F
  )
  pdf(paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_C_index/","OS_",TCGA_all$cancer_types[i],"_nomograph.pdf",sep = ""),
      width = 18,height = 10)
  plot(nom_cph)
  dev.off()
  
}

write.table(OS_C_index_value,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_C_index/OS_C_index.txt",
            sep = "\t",row.names = F,col.names = T,quote=FALSE)


#delete the features has small efficients(small range)----------------------------------------------
read.table("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_C_index/OS_C_index.txt",
           sep = "\t",header = T)->All_max_features

OS_C_index_value_2  <- data.frame(cancer_types=NA,features=NA,Max_OS_C_index=NA,features_2=NA,Max_OS_C_index_2=NA)

for (i in 1:(nrow(TCGA_all))) {
  
  #one cancer's significant features
  dplyr::filter(OS_multi_res,cancer_types == TCGA_all$cancer_types[i])%>%
    dplyr::select(Features)%>%as.matrix()%>%as.character()->i_all_features
  gsub("_high","",i_all_features)-> i_all_features
  
  
  #stage,survival,expr,TIL =>all_sigs
  TCGA_all[[4]][[i]]%>%as.data.frame() -> OS_stage
  OS_stage$Age <- as.numeric(OS_stage$Age)
  gsub("^stage x$","NA",OS_stage$Stage) %>%
    gsub("^i/ii nos$","NA",.) %>% gsub("^is$","NA",.)%>% 
    gsub("^stage 0$","0",.) %>% 
    gsub("^stage ib$","I",.) %>% gsub("^stage i$","I",.) %>% gsub("^stage ia$","I",.)%>%
    gsub("^stage ii$","II",.) %>% gsub("^stage iic$","II",.) %>% gsub("^stage iib$","II",.) %>% gsub("^stage iia$","II",.) %>% 
    gsub("^stage iiic$","III",.) %>% gsub("^stage iii$","III",.) %>% gsub("^stage iiib$","III",.) %>% gsub("^stage iiia$","III",.) %>% 
    gsub("^stage iv$","IV",.) %>% gsub("^stage iva$","IV",.) %>% gsub("^stage ivb$","IV",.) %>% gsub("^stage ivc$","IV",.) ->OS_stage$Stage
  dplyr::select(OS_stage,barcode,OS,Status,Stage,intersect(i_all_features,colnames(OS_stage)))%>%
    dplyr::filter(Stage != "NA")->OS_stage
  
  
  TCGA_all[[2]][[i]]%>%as.data.frame()  -> TIL
  dplyr::select(TIL,barcode,intersect(i_all_features,colnames(TIL))) ->TIL
  if(ncol(TIL) > 1){
    split_values(TIL) ->split_res
    split_res[[1]]->TIL
    split_res[[2]]->i_median_values
  }
  
  
  TCGA_all[[3]][[i]]%>%as.data.frame() -> Expr
  rownames(Expr) <- Expr$symbol
  Expr[,-1] ->Expr
  gsub("-", ".", rownames(Expr),fixed = TRUE) ->rownames(Expr)
  t(Expr)%>%as.data.frame()%>%
    dplyr::select(intersect(i_all_features,rownames(Expr)))%>%
    tibble::rownames_to_column()->Expr
  colnames(Expr)[1] <- "rowname"
  
  if(ncol(Expr) > 1){
    Samples_ID <- Expr$rowname
    Features <- colnames(Expr)
    Expr[,-1]%>%as.data.frame() -> Expr
    log2(Expr + 0.01)->Expr
    cbind(rowname = Samples_ID,Expr)->Expr
    colnames(Expr) <- Features
    
    split_values(Expr)[[1]] -> Expr
  }
  
  merge(OS_stage,TIL) %>%
    merge(Expr,by.x = "barcode",by.y = "rowname")->all_sigs
  dplyr::select(all_sigs,-barcode)->all_sigs
  all_sigs$OS=as.numeric(all_sigs$OS)
  all_sigs$Status=as.numeric(all_sigs$Status)
  
  selected_features <- strsplit(All_max_features$features[i],"[ + ]")%>%unlist()%>%unique()
  all_sigs %>% dplyr::select(OS,Status,Stage,intersect(colnames(all_sigs),selected_features))->Filtered_sigs
  
  ddist <- datadist(Filtered_sigs)
  options(datadist='ddist')

  
  if(max(Filtered_sigs$OS,na.rm = TRUE) >= 1825){
    
    f_cph <- cph(Surv(OS,Status) ~., data=Filtered_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  Filtered_sigs)[[1]]->i_C_index
    surv <- Survival(f_cph)
    nom_cph <- nomogram(f_cph, fun=list(function(x) surv(1095, x),function(x) surv(1825, x)),
                        fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                        funlabel=c("3-Year-Survival","5-Year-Survival"),
                        lp=F
    )
  }else{
    
    f_cph <- cph(Surv(OS,Status) ~., data=Filtered_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1095)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  Filtered_sigs)[[1]]->i_C_index
    surv <- Survival(f_cph)
    nom_cph <- nomogram(f_cph, fun=list(function(x) surv(365, x),function(x) surv(1095, x)),
                        fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                        funlabel=c("1-Year-Survival","3-Year-Survival"),
                        lp=F
    )
  }
  
  #select features that has (coef >= 0.2*Stage: can not control age's effect)
  # max(f_cph$coefficients[grep("Stage",names(f_cph$coefficients))])->max_stage_coef
  # abs(f_cph$coefficients[setdiff(1:length(f_cph$coefficients),grep("Stage",names(f_cph$coefficients)))]) ->other_coef
  # sapply(names(other_coef[which(other_coef[] >= 0.2*max_stage_coef)]), 
  #        function(x) unlist(strsplit(x,"="))[1])%>%as.character()->selected_features_2

  #select features that has range >= 0.2*Stage range
  max(nom_cph$Stage$points)-min(nom_cph$Stage$points) ->Stage_range
  small_sigs <- character()
  for (j in c(2:(length(nom_cph)-3))) {
    max(nom_cph[[j]]$points)-min(nom_cph[[j]]$points) -> j_Range
    if(j_Range < Stage_range*0.2){
      small_sigs <- c(small_sigs,names(nom_cph)[[j]])
    }
  }
  
  Filtered_sigs %>% dplyr::select(-small_sigs)->Filtered_sigs_2
  
  ddist <- datadist(Filtered_sigs_2)
  options(datadist='ddist')

  if(max(Filtered_sigs$OS,na.rm = TRUE) >= 1825){
    
    f_cph_2 <- cph(Surv(OS,Status) ~., data=Filtered_sigs_2,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph_2), data =  Filtered_sigs_2)[[1]]->i_C_index_2
    surv <- Survival(f_cph_2)
    nom_cph_2 <- nomogram(f_cph_2, fun=list(function(x) surv(1095, x),function(x) surv(1825, x)),
                        fun.at = c(0.1,0.3,0.5,0.7,0.9),
                        funlabel=c("3-Year-Survival","5-Year-Survival"),
                        lp=F
    )
    
  }else{
    
    f_cph_2 <- cph(Surv(OS,Status) ~., data=Filtered_sigs_2,surv=TRUE,x=TRUE, y=TRUE,time.inc=1095)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph_2), data =  Filtered_sigs_2)[[1]]->i_C_index_2
    surv <- Survival(f_cph_2)
    nom_cph_2 <- nomogram(f_cph_2, fun=list(function(x) surv(365, x),function(x) surv(1095, x)),
                        fun.at = c(0.1,0.3,0.5,0.7,0.9),
                        funlabel=c("1-Year-Survival","3-Year-Survival"),
                        lp=F
    )
  }
  
  write.table(Filtered_sigs_2,paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_C_index/","OS_",TCGA_all$cancer_types[i],"_features_2.txt",sep = ""),
              sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  
  pdf(paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_C_index/","OS_",TCGA_all$cancer_types[i],"_nomograph_2.pdf",sep = ""))
  plot(nom_cph_2,xfrac=.45, cex.axis=1, tcl=0.25)
  #xfrac= length from features to the right plot; cex.axis= size of number up on the axis; tcl=height of axis
  dev.off()
  
  
  OS_C_index_value_2[i,] <- c(TCGA_all$cancer_types[i],paste(intersect(colnames(all_sigs),selected_features),collapse = " + "),i_C_index,
                            paste(colnames(Filtered_sigs_2)[-c(1,2,3)],collapse = " + "),i_C_index_2)
  
}

write.table(OS_C_index_value_2,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_C_index/OS_C_index_2.txt",
            sep = "\t",row.names = F,col.names = T,quote=FALSE)
