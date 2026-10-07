#Open library of functions.
library(ggplot2)
library(pracma)

#Open folder with data. Read file names for all samples.
folder_path="/Users/xiongxia//Desktop/SMX+TMP_acetate_TZ060624"
setwd(folder_path)
file_names=list.files(folder_path)

#Make a whole file
total_mAU=c()
total_time=c()
total_ID=c()
for (file_name in file_names) {
  HPLC_sample_data=read.csv(file_name)
  total_mAU=c(total_mAU,HPLC_sample_data[,2])
  total_time=c(total_time,HPLC_sample_data[,1])
  total_ID=c(total_ID,rep(file_name,nrow(HPLC_sample_data)))
}

#Make a whole file including data from all sample.
whole_file=data.frame(Time=total_time,mAU=total_mAU,Sample_ID=total_ID)
ggplot(whole_file,aes(x=Time,y=mAU))+geom_point()
write.table(whole_file,"TZ06062024_HPLC_data_file.txt",sep='\t')
#whole_file=read.delim("TZ05212024_HPLC_data_file.txt",header=TRUE)

#Process standard curve samples.
standard_curve=file_names[15:23]
standard_curve_data=whole_file[whole_file$Sample_ID %in% standard_curve,]
ggplot(standard_curve_data,aes(x=Time,y=mAU,color=Sample_ID))+geom_line()

#Based on the most concentrated acetate standard sample, we 
#determined the start and end of the retention time for acetate.
start_time=5.6195
end_time=8.2
peak_duration=end_time-start_time
total_AUC_data=c();debackground_peak_AUC_list=c()
for (file_name in file_names) {
  
  #This block of code calculates the area under the curve (AUC) for the acetate peak
  integrated_area=0
  HPLC_sample_data=whole_file[whole_file$Sample_ID==file_name,]
  acetate_peak=subset(HPLC_sample_data, Time >= start_time & Time <= end_time)
  integrated_area=trapz(acetate_peak$Time, acetate_peak$mAU)
  total_AUC_data=c(total_AUC_data,integrated_area)
  
  #This next block of code removes the background below the peak.
  starting_peak_value=acetate_peak[1,2]
  ending_peak_value=acetate_peak[nrow(acetate_peak),2]
  background_area=(starting_peak_value+ending_peak_value)*peak_duration/2
  debackground_peak_AUC=integrated_area-max(0,background_area)
  debackground_peak_AUC_list=c(debackground_peak_AUC_list,debackground_peak_AUC)
}
#Save the de-backgrounded peak AUC data.
write.table(debackground_peak_AUC_list,"AUC_result_acetate.txt",sep='\t') 


#From this point on, we plot the whole peak data for everyone.
single_file=whole_file[whole_file$Sample_ID=="TZ01.CSV" | whole_file$Sample_ID=="TZ04.CSV" | whole_file$Sample_ID=="TZ07.CSV" | whole_file$Sample_ID=="TZ10.CSV", ]
single_file=whole_file[whole_file$Sample_ID=="TZ02.CSV" | whole_file$Sample_ID=="TZ05.CSV" | whole_file$Sample_ID=="TZ08.CSV" | whole_file$Sample_ID=="TZ11.CSV", ]
single_file=whole_file[whole_file$Sample_ID=="TZ03.CSV" | whole_file$Sample_ID=="TZ06.CSV" | whole_file$Sample_ID=="TZ09.CSV" | whole_file$Sample_ID=="TZ12.CSV", ]

ggplot(single_file,aes(x=Time,y=mAU,color=Sample_ID))+
  geom_line()+
  ylim(0,17)+
  geom_vline(xintercept = start_time, linetype="dotted", color = "blue", size=1.5)+
  geom_vline(xintercept = end_time, linetype="dotted", color = "blue", size=1.5)

#plot average
AUC=read.delim("acetate_run_data1.txt",header=TRUE)
AUC=read.delim("after_10_min_peak_AUC.txt",header=TRUE)
boxplot(AUC.OD~SAMPLE,data=AUC)
title("After_10_min_peak")
TukeyHSD(aov(AUC.OD~SAMPLE,data=AUC))
subset_AUC=AUC[AUC$SAMPLE=="NO" | AUC$SAMPLE=="SMX+TMP",]
summary(aov(AUC.OD~SAMPLE,data=subset_AUC))
