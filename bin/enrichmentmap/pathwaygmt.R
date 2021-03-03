#!/usr/bin/env Rscript

# input parameters
args = commandArgs(trailingOnly=TRUE)
genes_file = args[1]
out = args[2]

# load libraries
library(graflex)
library(org.Hs.eg.db)
library(GO.db)

# read genes
genes_df = read.csv(genes_file, sep="\t")
mygenes = genes_df$ENSEMBL
mygenes2 =  mapIds(org.Hs.eg.db, keys=mygenes, keytype = "ENSEMBL", column="SYMBOL")

# pathway-gene matrix
K = getK(mygenes)
rownames(K) = mygenes2

# pathway short description
go = select(GO.db, keys=colnames(K), columns=c("TERM","ONTOLOGY"), keytype="GOID")
terms = go$TERM
go = go[which(go$ONTOLOGY == "BP"), ]

# write output
lines = c()
for (i in 1:ncol(K)){
    goid = colnames(K)[i]
    if (goid %in% go$GOID){
        term = terms[i]
        genes = names(K[which(K[,i]==1),i])
        genes = paste(genes, collapse="\t")
        line = paste(goid, term, genes, sep="\t")
        lines = c(lines, line)
    }
}
o<-file(out)
lines = paste(lines, collapse="\n")
writeLines(lines, o)
close(o)
