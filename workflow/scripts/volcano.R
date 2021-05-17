library(EnhancedVolcano)
data <- read.table("DESeq2.6194.dir/kallisto.gene.counts.matrix.control_vs_treated.DESeq2.DE_results", header = T, dec = '.', sep = "\t", stringsAsFactors = FALSE)

up<-dim(data[data$log2FoldChange>=-1 & data$padj<0.05,])[1]
down<-dim(data[data$log2FoldChange<=-1 & data$padj<0.05,])[1]
unchange<-dim(data)[1]-up-down

keyvals <- rep('black', nrow(data))

# set the base name/label as 'Mid'
names(keyvals) <- rep(paste('unchanged:',unchange), nrow(data))

# modify keyvals for transcripts with fold change > 2.5
keyvals[which(data$log2FoldChange >= 1 & data$padj<0.05)] <- 'red2'
names(keyvals)[which(data$log2FoldChange >= 1 & data$padj<0.05)] <- paste('up:',up)

# modify keyvals for transcripts with fold change < -2.5
keyvals[which(data$log2FoldChange <= -1 & data$padj<0.05)] <- 'green'
names(keyvals)[which(data$log2FoldChange <= -1 & data$padj<0.05)] <- paste('down:',down)

unique(names(keyvals))
keyvals

p <- EnhancedVolcano(data,
                     lab = rownames(data),
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-15,15),
                     ylim = c(0,20),
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote('-'~Log[10]~'(FDR)'),
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     transcriptPointSize = 1.0,
                     transcriptLabSize = 0,
                     colOverride = keyvals,
                     colAlpha = 1,
                     legendLabSize = 15,
                     legendIconSize = 5.0,
                     border = 'partial',
                     borderWidth = 1.5,
                     borderColour = 'black',
                     #legend = c(paste('unchanged:',unchange),paste('up:',up),paste('down:',down)),
                     legendPosition = 'right')

svg(filename = "volcano.svg", width = 10, height = 10)
p
dev.off()