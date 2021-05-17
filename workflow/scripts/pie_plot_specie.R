library(ggplot2)
data <- read.table("species.txt", header = F, col.names= c('Specie','Hits','pHits'), dec = '.', sep = "\t", stringsAsFactors = FALSE)
data<-as.data.frame(data)
head(data)
data$pHits <- as.numeric(data$Hits/sum(data$Hits))
head(data)

data<-data[order(data$pHits, decreasing = TRUE),]
head(data)
data2 <- head(data[order(data$pHits, decreasing = TRUE),], n=6)
t <- c("other",sum(data[c(7:dim(data)[1]),2]),sum(data[c(7:dim(data)[1]),2])/sum(data$Hits))
t
data2 <-rbind(data2,t)
head(data2, n=7)
data2$pHits<-as.numeric(data2$pHits)
data2<-as.data.frame(data2)
head(data2, n=7)

# Create a basic bar
pie = ggplot(data2, aes(x="", y=pHits, fill=Specie)) + geom_bar(stat="identity", width=1)
pie
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round(pHits*100), "%")), position = position_stack(vjust = 0.5), size=8)

# Add color scale (hex colors)
#pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999")) 

# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "Species distribution")

# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    legend.text = element_text(size = 15),
                                    plot.title = element_text(hjust = 0.5, color = "#666666"))
svg(filename = "species_dist.svg", width = 21, height = 9.5)
pie
dev.off()