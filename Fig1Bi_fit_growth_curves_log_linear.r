library(tidyverse)
rm(list = ls())
source("/Users/xiongxia/Desktop/Colony counting/log_linear_functions.R")
setwd("/Users/xiongxia/Desktop")
OD_path = "XX20260423_gr_soil_Gore.txt"

OD = data_normalization_by_first_data_point(OD_path, 1)
OD=read.delim(OD_path, header = FALSE, sep = "\t", row.names = 1)
OD = data_normalization_by_max_average_rate(OD, "negative")  

write.csv(OD,"XX20260423_gr_soil_Gore_normalized.csv")
# transpose it so wells are columns, then convert to dataframe
OD = t(OD)
OD = data.frame(OD)
names(OD)[1:3] = c("cycle","seconds","temp")
# convert to long format by "gathering"
OD = OD %>%
  gather(well, OD, -cycle, -seconds, -temp) %>%
  mutate(hour = seconds / 60 / 60)
#OD$OD = ifelse(OD$OD < 0.125, 0.125, OD$OD)


#Plot raw growth curves
ggplot(OD,aes(x=hour,y=log(OD),col=well))+
  geom_point()+
  guides(fill="none")

first=OD[OD$well=="Pseudomonas.ASV19",]
ggplot(first,aes(x=hour,y=log(OD),col=well))+
  geom_point()+
  guides(fill="none")


# fit piece-wise log-linear
OD = OD[OD$hour <15,]
growth_rate = OD %>%
  group_by(well) %>%
  summarize(model_fit = fit_loglinear(hour, log(OD), 200),
            fit_variable = c("growth_rate", "y0", "start", "end")) %>%
  ungroup() %>%
  pivot_wider(names_from = fit_variable, values_from = model_fit)

all_data = left_join(OD, growth_rate) %>%
  group_by(well) %>%
  mutate(OD_pred = log_linear(hour, growth_rate[1], y0[1],start[1], end[1])) %>%
  ungroup

all_data %>%
  #filter(hour < 10) %>%
  ggplot(aes(x = hour, y = exp(OD_pred), color = well))+
  #geom_line(aes(y = OD), linetype = "dashed", size = 1)+
  geom_line()+
  scale_y_log10()
write.table(all_data,"E_SgalK_drug_gr_TZ10042024.txt",sep="\t")
write.table(growth_rate,"XX20260108_gr.txt",sep="\t")

