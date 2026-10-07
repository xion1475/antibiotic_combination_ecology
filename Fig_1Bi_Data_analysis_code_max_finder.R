library(ggplot2)
setwd(("/Users/apple/Desktop/"))
#Fig2B predicting combo growth rate with single growth rate
total_gr=read.delim("FigS2B.txt")
mono_data=total_gr[total_gr$Ecology=="Monoculture",]
summary(lm(Combo_XY~Smaller_Single*Bigger_Single,data=mono_data))
summary(lm(Combo_XY~Smaller_Single,data=mono_data))
summary(lm(Combo_XY~Bigger_Single,data=mono_data))
summary(lm(Combo_XY~Predicted,data=mono_data))
mono_spear=cor.test(mono_data$Combo_XY,mono_data$Predicted,method="spearman",exact=FALSE)
ggplot(mono_data,aes(x=Smaller_Single,y=Combo_XY))+geom_point()+ggtitle("Monoculture")
ggplot(mono_data,aes(x=Predicted,y=Combo_XY))+geom_point()+ggtitle("Monoculture")

co_data=total_gr[total_gr$Ecology=="Mutualism",]
summary(lm(Combo_XY~Smaller_Single*Bigger_Single,data=co_data))     
summary(lm(Combo_XY~Smaller_Single,data=co_data))
summary(lm(Combo_XY~Bigger_Single,data=co_data))
summary(lm(Combo_XY~Predicted,data=co_data))
co_spear=cor.test(co_data$Combo_XY,co_data$Predicted,method="spearman",exact=FALSE)
ggplot(co_data,aes(x=Predicted,y=Combo_XY))+geom_point()+ggtitle("Mutualism")
ggplot(co_data,aes(x=Smaller_Single,y=Combo_XY))+geom_point()+ggtitle("Mutualism")


#This block of code converts all time point after the max fluorescent/OD value to be identical to that.
for (j in 1:1) {
  #source("/Users/apple/Desktop/Drug_Interaction_Cross_Feeding/Nikon_death_functions.R")
  setwd(("/Users/xiongxia/Desktop/"))
  raw_data=read.delim("XX20250822_non_grower_raw.txt", header = FALSE,row.names = 1)
  
  finalized_data=raw_data
  total_time_point=ncol(raw_data)
  for (i in 4:nrow(raw_data)) {
    max_FP=max(raw_data[i,])
    max_colm=which(raw_data[i,]==max_FP)
    finalized_data[i,max_colm:total_time_point]=max_FP
    finalized_data[i,] = finalized_data[i,] - finalized_data[i,1]
    finalized_data[i,] = ifelse(finalized_data[i,] < 0, 0.09, finalized_data[i,] + 0.09)
  }
  
  write.table(finalized_data,"pr.txt",sep='\t')
}

#This block of code calculates the interaction index in the Yeh et al. (2006) manner.
for (j in 1:1) {
  setwd(("/Users/apple/Desktop/"))
  source("Drug_Interaction_Cross_Feeding/Yeh_interaction_function.R")
  raw_data=read.delim("total_interaction_data.txt",header=TRUE)
  yeh_metric=c()
  
  for (i in 1:nrow(raw_data)) {
    g_X=min(1,raw_data[i,4])
    g_Y=min(1,raw_data[i,5])
    g_XY=raw_data[i,7]
    epsilon_1=Yeh_paper_interaction_method(g_X,g_Y,g_XY)
    #epsilon_1=Yeh_paper_interaction_method_simpler(g_X,g_Y,g_XY)
    yeh_metric=c(yeh_metric,epsilon_1)
  }
  write.table(yeh_metric,"yeh_epsilon2.txt",sep='\t')
  
}
#measured interaction value differences.
for (j in 1:1) {
  setwd(("/Users/apple/Desktop/"))
  divided_interaction_data=read.delim("Interaction_divided_all.txt",header=TRUE)
  for (i in 1:nrow(divided_interaction_data)) {
    a=data.frame(matrix(ncol = 2, nrow = 6))
    colnames(a) = c('Condition', 'Value')
    a[,1]=c("Mono","Mono","Mono","Co","Co","Co")
    a[1:3,2]=t(divided_interaction_data[i,2:4])
    a[4:6,2]=t(divided_interaction_data[i,5:7])
    testresult=t.test(Value~Condition, paired=TRUE,alternative="two.sided",data=a)
    print(testresult$p.value)
  }
}

#measured interaction type lists.
for (j in 1:1) {
  setwd(("/Users/apple/Desktop/"))
  divided_interaction_data=read.delim("Interaction_divided_all.txt",header=TRUE)
  for (i in 1:nrow(divided_interaction_data)) {
    #monoculture
    mono_mean=sum(divided_interaction_data[i,2:4])/3
    if (mono_mean>0) {
      mono=t.test(divided_interaction_data[i,2:4], mu=0.05,alternative="greater")
    } else {
      mono=t.test(divided_interaction_data[i,2:4], mu=-0.05,alternative="less")
    }
    #monoculture
    co_mean=sum(divided_interaction_data[i,5:7])/3
    if (co_mean>0) {
      co=t.test(divided_interaction_data[i,5:7], mu=0.05,alternative="greater")
    } else {
      co=t.test(divided_interaction_data[i,5:7], mu=-0.05,alternative="less")
    }
    #co=t.test(divided_interaction_data[i,5:7], mu=0,alternative="two.sided")
    print(c(mono$p.value,co$p.value))
  }
}



#FigS1C drug degradation
setwd("/Users/apple/Desktop/")
db=read.delim("total_S_drug_degrade.txt")

a=db[db$Antibiotic=="RIF",]
summary(aov(MIC~Antibiotic_condition,data=a))

#Fig4A ratio vs S_mutual
FLP_ratio_data=read.delim("FigSX_A.txt",header=TRUE)
FLP_ratio_data=read.delim("FigS7A.txt",header=TRUE)

no_drug=c(0.626195732, 0.545325779, 0.443737441)  
no_drug=c(-0.203289897,-0.263343972,-0.352873926)
for (i in 1:nrow(FLP_ratio_data)) {
  ratios=FLP_ratio_data[i,3:5]
  t_test_result=t.test(ratios,no_drug,alternative = "greater")
  result=t_test_result$p.value
  print(result)
}
summary(aov(p_val_rank~Interaction_different_on_ecology,data=FLP_ratio_data))
summary(aov(p_value_against_no_drug_one_sided~Interaction_different_on_ecology,data=FLP_ratio_data))
mean(FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="YES",]$p_value_against_no_drug_one_sided)
#[1] 0.1783798
mean(FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="NO",]$p_value_against_no_drug_one_sided)
#[1] 0.4514202

no_diff=FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="NO",]$Mean_FLP_ratio_rank
diff=FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="YES",]$Mean_FLP_ratio_rank
wilcox.test(diff,no_diff)

no_diff=FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="NO",]$Mean_YCFP_ratio_rank
diff=FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="YES",]$Mean_YCFP_ratio_rank
t.test(no_diff,diff)

diff=FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="YES",]$P_value_rank_log_YCFP
no_diff=FLP_ratio_data[FLP_ratio_data$Interaction_different_on_ecology=="NO",]$P_value_rank_log_YCFP

additive=FLP_ratio_data[FLP_ratio_data$Interaction_Class=="Additivity",]$Mean_YCFP_ratio_rank
antagonism=FLP_ratio_data[FLP_ratio_data$Interaction_Class=="Antagonism",]$Mean_YCFP_ratio_rank
wilcox.test(additive,antagonism)

#Fig1E Fisher's exact test on percentage distribution of drug interaction types.
db=data.frame(
  Antagonism=c(6,25),
  Synergy=c(14,6),
  Additivity=c(80,69),
  row.names=c("Monoculture","Mutualism"))
fisher.test(db)

diff_intxn_more_S=data.frame(
  moreS=c(8,4),
  not_more_S=c(6,16),
  row.names=c("diff_intxn","not_diff_intxn")
)

diff_intxn_more_S=data.frame(
  diff_int=c(8,9),
  not_diff_int=c(6,13),
  row.names=c("more_S","not_more_S")
)