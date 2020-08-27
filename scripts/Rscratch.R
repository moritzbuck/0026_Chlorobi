library(plyr)
library(pheatmap)
library(ggplot)
library(ggplot2)
read.csv("data/motu_table_cutoff.csv", row.names=1)
data = read.csv("data/motu_table_cutoff.csv", row.names=1)
md = read.csv("/home/moritz/kadath/data/data_submit/metadata/samples_contextual_working_copy.csv", row.names=1, comment.char="#")
md$region
region = md$region
names(regions) = row.names(md)
names(region) = row.names(md)
region
data
data[rowSums(data) > 0, colSums(data) > 0]
data = data[rowSums(data) > 0, colSums(data) > 0]
data
t(data)
data = t(data)
dim(data)
data = t(data)
tdat = t(data)
tdat$region = region[row.names(tdat)]
tdat
tdat = as.data.frame(t(data))
tdat$region = region[row.names(tdat)]
tdat
ddply(tdat, .(region), colwise(sum))
ddply(tdat > 0, .(region), colwise(sum))
tdat
tdat = as.data.frame(t(data) >0)
tdat =
tdat$region = region[row.names(tdat)]
ddply(tdat > 0, .(region), colwise(sum))
ddply(tdat, .(region), colwise(sum))
ddply(tdat, .(region), colwise(mean))
per_region = ddply(tdat, .(region), colwise(mean))
row.names(per_region) = per_region$region
per_region$region
md
md[is.na(md$region),]
row.names(md)
row.names(data)
names(data)
colnames(data)
tdat
colnames(tdat)
head(tdat)
row.names(tdat)
replace(".", "_", row.names(tdat))
gsub(".", "_", row.names(tdat))
gsub(".", "_", row.names(tdat), fixed=TRUE)
row.names(md)
data = read.csv("data/motu_table_cutoff.csv", row.names=1, as.is=TRUE)
colnames(data)
?read.csv
data = read.csv("data/motu_table_cutoff.csv", row.names=1, check.names=FALSE)
?read.csv
colnames(data)
data = data[rowSums(data) > 0, colSums(data) > 0]
tdat = as.data.frame(t(data) >0)
tdat$region = region[row.names(tdat)]
tdat$region
tdat[is.na(tdat$region),]
row.names(md)
row.names(tdat)
row.names(md)
gsub("_","-",row.names(md))
row.names(md) = gsub("_","-",row.names(md))
region = md$region
names(region) = row.names(md)
tdat$region = region[row.names(tdat)]
tdat[is.na(tdat$region),]
tdat[is.na(tdat$region),1:4]
options(width = 230)
tdat[is.na(tdat$region),1:4]
tdat[grep("^[A-I]-[0-9]$", row.names(tdat)),"region"]
tdat[grep("^[A-I]-[0-9]$", row.names(tdat)),]
tdat[!grep("^[A-I]-[0-9]$", row.names(tdat)),]
tdat[!grepl("^[A-I]-[0-9]$", row.names(tdat)),]
tdat = tdat[!grepl("^[A-I]-[0-9]$", row.names(tdat)),]
tdat
tdat[grep("^IH", row.names(tdat)),1:3]
tdat[grep("^IH", row.names(tdat)),"region"]

tdat = tdat[!grepl("^SB-21$", row.names(tdat)),]
tdat
per_region = ddply(tdat, .(region), colwise(mean))
row.names(per_region) = per_region$region
tdat[is.na(tdat$region),1:4]
tdat = tdat[!grepl("^SB-21", row.names(tdat)),]
tdat[is.na(tdat$region),1:4]
per_region = ddply(tdat, .(region), colwise(mean))
row.names(per_region) = per_region$region
per_region$region
head(per_region)
per_region[,-1]
pheatmap(per_region[,-1])
pheatmap((per_region[,-1] > 0) + 0)
savehistory("Rscratch.R")
