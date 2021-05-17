args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(ggplot2)
library(forcats)

WEGO <- read.table(args[1], header = FALSE, sep = "\t")
head(WEGO)


colnames(WEGO) <- c("ID", "Domain", "Number of genes", "Percentage of genes", "Description")
head(WEGO)

WEGO$`Percentage of genes` <- WEGO$`Percentage of genes`*100
head(WEGO)

WEGO2 <- WEGO %>%
## Group the entries by "Domain"
group_by(Domain) %>%
## Take the top 5 entries per "Domain" according to "Percentage of genes"
top_n(7, `Percentage of genes`) %>% 
## Ungroup the entries
ungroup() %>% 
## Arrange the entries by "Domain", then by "Percentage of genes"
arrange(Domain, `Percentage of genes`) %>% 
## Take note of the arrangement by creating a "Position" column
mutate(Position = n():1)

head(WEGO2)

normalizer <- max(WEGO2$`Number of genes`)/max(WEGO2$`Percentage of genes`)

p <- ggplot(data = WEGO2, aes(x = fct_reorder(Description, desc(Position)), y = `Percentage of genes`, fill = Domain)) +
## Plot "Description" in the x-axis following the order stated in the "Position" column
## vs normalized "Number of genes" in the second y-axis
geom_col(data = WEGO2, aes(x = fct_reorder(Description, desc(Position)), y = `Number of genes`/normalizer)) +
## Add a second y-axis based on the transformation of "Percentage of genes" to "Number of genes".
## Notice that the transformation undoes the normalization for the earlier geom_col.
scale_y_continuous(sec.axis = sec_axis(trans = ~.*normalizer, name = "Number of genes")) +
## Modify the aesthetic of the theme
theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15),
legend.text = element_text(size = 15), legend.title = element_text(size = 15),
legend.key.size =  unit(0.2, "in"), plot.title = element_text(size = 15, hjust = 0.5)) +
## Add a title to the plot
labs(x = NULL, title = "Gene Ontology (GO) Annotation")

svg(filename = "gene_per_go.svg", width = 21, height = 9.5)
p
dev.off()