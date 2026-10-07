library(ggplot2)
setwd("/Users/xiongxia/Desktop/")
heatmap_data=read.delim("FigS12D.txt")

ggplot(heatmap_data, aes(as.factor(folP), as.factor(folA), fill=E_co_gal)) + 
  geom_tile()+
  xlab("*folP* repression")+
  ylab("*folA* repression")+
  labs(fill="")+
  scale_fill_gradientn(colors = c("black", "red", "white")) +
  theme(
    panel.background = element_rect(fill = "grey", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = ggtext::element_markdown(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20),
    legend.title = element_text(size=20)
  )

ggplot(heatmap_data, aes(as.factor(folP), as.factor(folA), fill=S_co_gal)) + 
  geom_tile()+
  xlab("*folP* repression")+
  ylab("*folA* repression")+
  labs(fill="")+
  scale_fill_gradientn(colors = c("black", "red", "white")) +
  theme(
    panel.background = element_rect(fill = "grey", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = ggtext::element_markdown(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20),
    legend.title = element_text(size=20)
  )

ggplot(heatmap_data, aes(as.factor(folP), as.factor(folA), fill=E_co_ace)) + 
  geom_tile()+
  xlab("*folP* repression")+
  ylab("*folA* repression")+
  labs(fill="")+
  scale_fill_gradientn(colors = c("black", "red", "white")) +
  theme(
    panel.background = element_rect(fill = "grey", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = ggtext::element_markdown(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20),
    legend.title = element_text(size=20)
  )
ggplot(heatmap_data, aes(as.factor(folP), as.factor(folA), fill=S_co_ace)) + 
  geom_tile()+
  xlab("*folP* repression")+
  ylab("*folA* repression")+
  labs(fill="")+
  scale_fill_gradientn(colors = c("black", "red", "white")) +
  theme(
    panel.background = element_rect(fill = "grey", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = ggtext::element_markdown(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20),
    legend.title = element_text(size=20)
  )

ggplot(heatmap_data, aes(as.factor(folP), as.factor(folA), fill=S_co_total_C)) + 
  geom_tile()+
  xlab("*folP* repression")+
  ylab("*folA* repression")+
  labs(fill="")+
  scale_fill_gradientn(colors = c("black", "red", "white")) +
  theme(
    panel.background = element_rect(fill = "grey", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = ggtext::element_markdown(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20),
    legend.title = element_text(size=20)
  )
