args = commandArgs(trailingOnly=TRUE)
library(ggplot2)

data <- read.table(args[1], header = F, dec = '.', sep = "\t", stringsAsFactors = FALSE)

dat <- data.frame(
  FunctionClass = factor(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z"), levels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z")),
  Legend = c("A: RNA processing and modification", "B: Chromatin structure and dynamics", "C: Energy production and conversion", "D: Cell cycle control, cell division, chromosome partitioning", "E: Amino acid transport and metabolism", "F: Nucleotide transport and metabolism", "G: Carbohydrate transport and metabolism", "H: Coenzyme transport and metabolism", "I: Lipid transport and metabolism", "J: Translation, ribosomal structure and biogenesis", "K: Transcription", "L: Replication, recombination and repair", "M: Cell wall/membrane/envelope biogenesis", "N: Cell motility", "O: Posttranslational modification, protein turnover, chaperones", "P: Inorganic ion transport and metabolism", "Q: Secondary metabolites biosynthesis, transport and catabolism", "R: General function prediction only", "S: Function unknown", "T: Signal transduction mechanisms", "U: Intracellular trafficking, secretion, and vesicular transport", "V: Defense mechanisms", "W: Extracellular structures", "Y: Nuclear structure", "Z: Cytoskeleton"),
  Frequency=data[,2]
)

p <- ggplot(data=dat, aes(x=FunctionClass, y=Frequency, fill=Legend))+
  geom_bar(stat="identity", position=position_dodge(), colour="seashell")+
  guides(fill = guide_legend(ncol = 1))+
  theme(axis.line = element_blank(), axis.text = element_text(size = 15),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15), axis.title = element_text(size = 15))+
  xlab("Factor Class")


svg(filename = "kog_out.svg", width = 21, height = 21)
p
dev.off()