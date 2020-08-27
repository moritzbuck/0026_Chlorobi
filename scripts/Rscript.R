library(plyr)
library(pheatmap)
library(ggplot)
library(ggplot2)

options(width = 230)

data = read.csv("data/motu_table_cutoff.csv", row.names=1, check.names=FALSE)
data = data[rowSums(data) > 0, colSums(data) > 0]

md = read.csv("/home/moritz/kadath/data/data_submit/metadata/samples_contextual_working_copy.csv", row.names=1, comment.char="#")
row.names(md) = gsub("_","-",row.names(md))

tdat = as.data.frame(t(data) >0)
tdat$region = region[row.names(tdat)]
tdat = tdat[!grepl("^[A-I]-[0-9]$", row.names(tdat)),]
tdat = tdat[!grepl("^SB-21", row.names(tdat)),]
tdat[grep("^IH", row.names(tdat)),"region"] = "Wisconsin"
tdat[grep("^SRR", row.names(tdat)),"region"] = "Kenora"
tdat[grep("^SA2", row.names(tdat)),"region"] = "Äspö"
per_region = ddply(tdat, .(region), colwise(mean))

row.names(per_region) = per_region$region

pheatmap(per_region[,-1])
pheatmap((per_region[,-1] > 0) + 0)
