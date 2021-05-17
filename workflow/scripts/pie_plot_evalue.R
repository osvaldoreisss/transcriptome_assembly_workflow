library(ggplot2)
data <- read.table("evalue.txt", header = F, col.names= c('Evalue','value'), dec = '.', sep = "\t", stringsAsFactors = FALSE)
data<-as.data.frame(data)
head(data)
data$pEvalue <- as.numeric(data$value/sum(data$value))
head(data)

# Create a basic bar
pie = ggplot(data, aes(x="", y=pEvalue, fill=Evalue)) + geom_bar(stat="identity", width=1)

# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round(pEvalue*100), "%")), position = position_stack(vjust = 0.5), size=8)

# Add color scale (hex colors)
#pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999")) 

# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "E-value distribution")

# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    legend.text = element_text(size = 15),
                                    plot.title = element_text(hjust = 0.5, color = "#666666"))
svg(filename = "evalue_dist.svg", width = 21, height = 9.5)
pie
dev.off()