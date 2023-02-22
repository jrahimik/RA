## Get plot
library(ggplot2)
library(RColorBrewer)

h2_data <- read.table("Data/h2Estimates.csv",header=F,sep=",")
h2_data <- h2_data[,c(1:3)]
colnames(h2_data)<- c("category","estimation","std_dev")

z <- h2_data$estimation/h2_data$std_dev
p_values_data <- 2 * (1 - pnorm(abs(z)))
h2_data$pvalue <- p_values_data

#p_values_data <- c(0.01,1,1,1,1,0.01,0.01,1,1,1,1)

# p_values_data <- data.frame(category = c("A","B","C"), p_value = c(0.1, 0.05, 0.01))
# data <- data.frame(category = c("A","B","C"), estimation = c(0.1, 0.05, 0.01),std_dev=c(0.01,0.02,0.03))

map_data <- data.frame(number=c(1,2,3,4,5,6,7,8,9,10,11,12,13),name=c("CYB561","HLA","ACVR1B","CCL","IL1","AKR1B10","CLEC18A","CDC42BPB","IL15","CCRL2","COL6A1","FIGNL1","CHRNA3"))



ggplot(h2_data, aes(x =category , y = estimation,fill=category)) +
  geom_bar(stat="identity", width=0.6) +
  geom_errorbar(aes(ymin=estimation, ymax=estimation+std_dev), width=0.2, size=0.5) +
  annotate("text", x = 1:13, y = 0.12, label = ifelse(h2_data$pvalue < 0.13, "*", ""), size = 6)+
  annotate("text", x = 1:13, y = 0.15, label = ifelse(h2_data$pvalue < 0.05, "**", ""), size = 6,color = "red")+
  scale_fill_gradientn(colors = rainbow(length(h2_data$category)))+guides(fill=FALSE)+
  geom_text(aes(label = map_data$name), size = 4, hjust = 0.5, vjust = -0.5, angle = 90)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(x = "Module number", y = "h2")+scale_x_continuous(breaks = 1:13)


