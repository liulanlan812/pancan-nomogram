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

# #bootstrap adjust c-index
# v <- validate(f_cph_0, dxy=TRUE, B=1000)
# Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
# orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]
# bias_corrected_c_index  <- abs(Dxy)/2+0.5
# orig_c_index <- abs(orig_Dxy)/2+0.5



# #compare two model
# rcorrp.senc()
# r<-rcorrp.senc(x1,x1,y)
# pValue<-1-pnorm((r[11]-r[12])/(r[2]/r[5])*1.96) 

#OS----------------------------------------------------------------------

readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/muli_uni_survival/OS_multi_res.tsv")->OS_multi_res

TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_TIL_exp_data.rds.gz")%>%
  dplyr::filter(cancer_types %in% unique(OS_multi_res$cancer_types))

OS_C_index_value <- data.frame(cancer_types=NA,features=NA,OS_C_index=NA)

for (i in 1:nrow(TCGA_all)) {
  
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
    
    split_values(Expr) -> split_res
    split_res[[1]]->Expr
    rbind(i_median_values,split_res[[2]])->i_median_values
  }
  i_median_values$cancer_types <- rep(TCGA_all$cancer_types[i],nrow(i_median_values))
  
  if(i==1){
    Median_values <- i_median_values
  }else{
    Median_values <- rbind(Median_values,i_median_values)
  }
  
  merge(OS_stage,TIL) %>%
    merge(Expr,by.x = "barcode",by.y = "rowname")->all_sigs

  # write.table(all_sigs,
  #             paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/Cancer_sig_features/","OS_",TCGA_all$cancer_types[i],"_all_sigs.txt",sep = ""),
  #             sep = "\t",row.names = F,col.names = T,quote=FALSE)

  
  # model
  all_sigs$barcode -> i_samples
  dplyr::select(all_sigs,-barcode)->all_sigs
  
  ddist <- datadist(all_sigs)
  options(datadist='ddist')
  all_sigs$OS=as.numeric(all_sigs$OS)
  all_sigs$Status=as.numeric(all_sigs$Status)
  
  # if(max(all_sigs$OS,na.rm = T) < 1825){
  #   
  #   f_cph_0 <- cph(Surv(OS,Status) ~., data=all_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1095)
  #   1-rcorrcens(Surv(OS,Status) ~ predict(f_cph_0), data =  all_sigs)[[1]] ->C_index_0
  #   print(C_index_0)
  #   
  #   surv_0 <- Survival(f_cph_0)
  #   nom_cph_0 <- nomogram(f_cph_0, fun=list(function(x) surv_0(365, x),function(x) surv_0(1095, x)),
  #                         fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
  #                         funlabel=c("1-Year-Survival","3-Year-Survival"),
  #                         lp=F
  #   )
  #   
  # }else{
  #   f_cph_0 <- cph(Surv(OS,Status) ~., data=all_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825)
  #   1-rcorrcens(Surv(OS,Status) ~ predict(f_cph_0), data =  all_sigs)[[1]] ->C_index_0
  #   print(C_index_0)
  #   
  #   surv_0 <- Survival(f_cph_0)
  #   nom_cph_0 <- nomogram(f_cph_0, fun=list(function(x) surv_0(1095, x),function(x) surv_0(1825, x)),
  #                         fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
  #                         funlabel=c("3-Year-Survival","5-Year-Survival"),
  #                         lp=F
  #   )
  # }
  
  
  #filter too small contribution features

  small_sigs <- character()
  for (m in 4:ncol(all_sigs)) {
    
    selected_features <- colnames(all_sigs)[m]
    all_sigs %>% dplyr::select(OS,Status,Stage,selected_features)->Filtered_sigs
    
    ddist <- datadist(Filtered_sigs)
    options(datadist='ddist')
    
    f_cph <- cph(Surv(OS,Status) ~., data=Filtered_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825)
    1-rcorrcens(Surv(OS,Status) ~ predict(f_cph), data =  Filtered_sigs)[[1]]
    
    surv <- Survival(f_cph)
    nom_cph <- nomogram(f_cph, fun=list(function(x) surv(1095, x),function(x) surv(1825, x)),
                        fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                        funlabel=c("3-Year-Survival","5-Year-Survival"),
                        lp=F
    )
    
    max(nom_cph$Stage$points)-min(nom_cph$Stage$points) ->Stage_range
    
    for (j in c(2:(length(nom_cph)-3))) {
      max(nom_cph[[j]][[3]])-min(nom_cph[[j]][[3]]) -> j_Range
      if(j_Range < Stage_range*0.3){
        small_sigs <- c(small_sigs,names(nom_cph)[[j]])
      }
    }
    
  }

  
  all_sigs %>% dplyr::select(-small_sigs)->Filtered_sigs

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
  
  OS_C_index_value[i,] <- c(TCGA_all$cancer_types[i],paste(colnames(Filtered_sigs)[-c(1,2)],collapse = " + "),i_C_index)

  pdf(paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_point_range/","OS_",TCGA_all$cancer_types[i],"_nomograph.pdf",sep = ""),
      width = 18,height = 10)
  plot(nom_cph)
  dev.off()

  
}
write.table(OS_C_index_value,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_point_range/OS_C_index.txt",
            sep = "\t",row.names = F,col.names = T,quote=FALSE)

write.table(Median_values,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/OS_Median_values.txt",
            sep = "\t",row.names = F,col.names = T,quote=FALSE)
write.table(OS_C_index_value,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/OS_C_index.txt",
            sep = "\t",row.names = F,col.names = T,quote=FALSE)

# readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/OS_C_index.txt")-> C_index
# colnames(C_index)[3] <- "OS_C_index"
# 
# pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/OS_C_index.pdf",width = 8,height = 6)
# ggplot(C_index,aes(cancer_types,OS_C_index))+ ylim(0.5,1) +
#   geom_point() + 
#   theme_bw()
# #theme(panel.grid =element_blank())
# dev.off()


# cal2 <- calibrate(f_cph_0, cmethod='KM', method="boot", u=3650, m=round(nrow(all_sigs)/4), b= nrow(all_sigs))
# 
# par(mar=c(8,5,3,2),cex = 1.0)
# plot(cal2,lwd=2,lty=1,
#      #errbar.col=c(rgb(0,118,192,maxColorValue=255)),
#      xlim=c(0,1.0),ylim=c(0,1.0),
#      xlab="Nomogram-Predicted Probability of 10-Year OS",
#      ylab="Actual 10-Year OS",
#      col=c(rgb(192,98,83,maxColorValue=255)),subtitles=FALSE,riskdist=FALSE)
# abline(0,1,lty =3)


#PFS----------------------------------------------------------
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/muli_uni_survival/PFS_multi_res.tsv")%>%
  dplyr::filter(coxp < 0.1)->PFS_multi_res

TCGA_all <- readr::read_rds("/home/liull/TCGA_nomograph/TCGA_result/TCGA_combined_clinical_TIL_exp_data.rds.gz")%>%
  dplyr::filter(cancer_types %in% unique(PFS_multi_res$cancer_types))

PFS_C_index_value <- data.frame(cancer_types=NA,features=NA,PFS_C_index=NA)

#i=17 SKCM,only one nTreg,can not plot
#i=19 THCA,only one KIR2DL1,can not plot
for (i in 1:nrow(TCGA_all)) {
  
  #one cancer's significant features
  dplyr::filter(PFS_multi_res,cancer_types == TCGA_all$cancer_types[i])%>%
    dplyr::select(Features)%>%as.matrix()%>%as.character()->i_all_features
  gsub(".high","",i_all_features)-> i_all_features
  
  
  #stage,survival,expr,TIL =>all_sigs
  TCGA_all[[4]][[i]]%>%as.data.frame()%>% dplyr::select(barcode,Age,Stage) -> stage_status
  TCGA_all[[5]][[i]]%>%as.data.frame() -> PFS_status
  dplyr::left_join(PFS_status,stage_status)->PFS_stage
  
  PFS_stage$Age <- as.numeric(PFS_stage$Age)
  gsub("^stage x$","NA",PFS_stage$Stage) %>%
    gsub("^i/ii nos$","NA",.) %>% gsub("^is$","NA",.)%>% 
    gsub("^stage 0$","0",.) %>% 
    gsub("^stage ib$","I",.) %>% gsub("^stage i$","I",.) %>% gsub("^stage ia$","I",.)%>%
    gsub("^stage ii$","II",.) %>% gsub("^stage iic$","II",.) %>% gsub("^stage iib$","II",.) %>% gsub("^stage iia$","II",.) %>% 
    gsub("^stage iiic$","III",.) %>% gsub("^stage iii$","III",.) %>% gsub("^stage iiib$","III",.) %>% gsub("^stage iiia$","III",.) %>% 
    gsub("^stage iv$","IV",.) %>% gsub("^stage iva$","IV",.) %>% gsub("^stage ivb$","IV",.) %>% gsub("^stage ivc$","IV",.) ->PFS_stage$Stage
  
  
  if(length(grep(pattern="Stage", x=i_all_features, value=TRUE)) != 0){
    dplyr::select(PFS_stage,barcode,PFS,PFS.time,Stage,intersect(i_all_features,colnames(PFS_stage))) ->PFS_stage
  }else{
    dplyr::select(PFS_stage,barcode,PFS,PFS.time,intersect(i_all_features,colnames(PFS_stage))) ->PFS_stage
  }
  
  
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
    
    split_values(Expr) -> split_res
    split_res[[1]]->Expr
    rbind(i_median_values,split_res[[2]])->i_median_values
  }
  i_median_values$cancer_types <- rep(TCGA_all$cancer_types[i],nrow(i_median_values))
  
  if(i==1){
    Median_values <- i_median_values
  }else{
    Median_values <- rbind(Median_values,i_median_values)
  }
  
  merge(PFS_stage,TIL) %>%
    merge(Expr,by.x = "barcode",by.y = "rowname")->all_sigs
  
  write.table(all_sigs,
              paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/Cancer_sig_features/","PFS_",TCGA_all$cancer_types[i],"_all_sigs.txt",sep = ""),
              sep = "\t",row.names = F,col.names = T,quote=FALSE)
  
  
  # model
  all_sigs$barcode -> i_samples
  dplyr::select(all_sigs,-barcode)->all_sigs
  
  ddist <- datadist(all_sigs)
  options(datadist='ddist')
  all_sigs$PFS.time=as.numeric(all_sigs$PFS.time)
  all_sigs$PFS=as.numeric(all_sigs$PFS)
  
  if(max(all_sigs$PFS.time,na.rm = T) < 1825){
    
    f_cph_0 <- cph(Surv(PFS.time,PFS) ~., data=all_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1095)
    1-rcorrcens(Surv(PFS.time,PFS) ~ predict(f_cph_0), data =  all_sigs)[[1]] ->C_index_0
    #print(C_index_0)
    
    surv_0 <- Survival(f_cph_0)
    nom_cph_0 <- nomogram(f_cph_0, fun=list(function(x) surv_0(365, x),function(x) surv_0(1095, x)),
                          fun.at = c(0.1,0.3,0.5,0.7,0.9),
                          funlabel=c("1-Year-Survival","3-Year-Survival"),
                          lp=F
    )
    
  }else if(max(all_sigs$PFS.time,na.rm = T) >= 3650){
    f_cph_0 <- cph(Surv(PFS.time,PFS) ~ . , data=all_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=3650)
    1-rcorrcens(Surv(PFS.time,PFS) ~ predict(f_cph_0), data =  all_sigs)[[1]] ->C_index_0
    #print(C_index_0)
    
    surv_0 <- Survival(f_cph_0)
    nom_cph_0 <- nomogram(f_cph_0, fun=list(function(x) surv_0(1095, x),function(x) surv_0(1825, x),function(x) surv_0(3650, x)),
                          fun.at = c(0.1,0.3,0.5,0.7,0.9),
                          funlabel=c("3-Year-Survival","5-Year-Survival","10-Year-Survival"),
                          lp=F
    )
  }else{
    f_cph_0 <- cph(Surv(PFS.time,PFS) ~., data=all_sigs,surv=TRUE,x=TRUE, y=TRUE,time.inc=1825)
    1-rcorrcens(Surv(PFS.time,PFS) ~ predict(f_cph_0), data =  all_sigs)[[1]] ->C_index_0
    #print(C_index_0)
    
    surv_0 <- Survival(f_cph_0)
    nom_cph_0 <- nomogram(f_cph_0, fun=list(function(x) surv_0(365, x),function(x) surv_0(1095, x),function(x) surv_0(1825, x)),
                          fun.at = c(0.1,0.3,0.5,0.7,0.9),
                          funlabel=c("1-Year-Survival","3-Year-Survival","5-Year-Survival"),
                          lp=F
    )
  }
  # readr::write_rds(f_cph_0,
  #                  paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/Cancer_sig_features/",
  #                        "PFS_",TCGA_all$cancer_types[i],"_model.rds.gz",sep = ""))
  
  PFS_C_index_value[i,] <- c(TCGA_all$cancer_types[i],paste(i_all_features,collapse = " + "),C_index_0)
  
  pdf(paste("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot/","PFS_",TCGA_all$cancer_types[i],"_nomograph.pdf",sep = ""),
      width = 18,height = 10)
  plot(nom_cph_0)
  dev.off()
  
  
}

write.table(Median_values,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/PFS_Median_values.txt",
           sep = "\t",row.names = F,col.names = T,quote=FALSE)
write.table(PFS_C_index_value,"/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/PFS_C_index.txt",
            sep = "\t",row.names = F,col.names = T,quote=FALSE)

# readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/PFS_C_index.txt")-> C_index
# colnames(C_index)[3] <- "PFS_C_index"
# 
# pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/PFS_C_index.pdf",width = 8,height = 6)
# ggplot(C_index,aes(cancer_types,PFS_C_index))+ ylim(0.5,1) +
#   geom_point() + 
#   theme_bw()
# #theme(panel.grid =element_blank())
# dev.off()

#OS PFS c-index compare-----------------------------------------------------------
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/OS_C_index.txt")-> OS_C_index
readr::read_tsv("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/PFS_C_index.txt")-> PFS_C_index
colnames(OS_C_index)[3] <- "C_index"
colnames(PFS_C_index)[3] <- "C_index"

intersect(OS_C_index$cancer_types,PFS_C_index$cancer_types)->all_cancers
dplyr::filter(OS_C_index,cancer_types %in% all_cancers)->OS_C_index
dplyr::filter(PFS_C_index,cancer_types %in% all_cancers)->PFS_C_index

rbind(cbind(Type=rep("OS",nrow(OS_C_index)),OS_C_index),cbind(Type=rep("PFS",nrow(PFS_C_index)),PFS_C_index)) %>%dplyr::select(-features)->all_index
pdf("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/All_C_index.pdf",width = 12,height = 6)
ggplot(all_index, aes(x=cancer_types,y=C_index,colour=Type,group=Type))+ geom_line() + theme_bw()
dev.off()
