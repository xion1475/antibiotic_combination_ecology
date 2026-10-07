#CFU-based species ratio statistical tests. Fig4
setwd("/Users/apple/Desktop/")
library(ggplot2)
library(pracma)

fold_change_log_ratio_data=read.delim("Fig4A_log_fold_change_comp.txt",header=TRUE)
#row.names(fold_change_log_ratio_data)=fold_change_log_ratio_data$Antibiotic_condition

for (i in 2:ncol(fold_change_log_ratio_data)) {
  current_drugs=colnames(fold_change_log_ratio_data)[i]
  print(current_drugs)
  #a=TukeyHSD(aov(fold_change_log_ratio_data[,i]~fold_change_log_ratio_data[,1]))
  a=t.test(fold_change_log_ratio_data[1:3,i],fold_change_log_ratio_data[7:9,i]) #galK only
  print(a)
}
p.adjust(p_value)
