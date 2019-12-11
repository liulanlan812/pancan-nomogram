#survival plots
#move to .3
library(magrittr)
library(survminer)
library(survival)
options(stringsAsFactors = FALSE)

list.files("/home/liull/TCGA_nomograph/TCGA_result/nomograph_model/nomograph_plot_C_index",
           pattern = "features_2.txt",full.names = TRUE)->all_features_file

for (i in 1:length(all_features_file)) {
  
  read.table(all_features_file[i],header = T,sep = "\t")->i_data
  i_data$Status[which(i_data$Status[] == 1)] <- 2
  i_data$Status[which(i_data$Status[] == 0)] <- 1
  
  
  survfit(Surv(OS, Status) ~Stage, data=i_data) ->Stage_fit
  ggsurvplot(Stage_fit, data = i_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
}