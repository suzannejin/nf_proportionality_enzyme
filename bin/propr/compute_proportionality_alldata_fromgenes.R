# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option(c("--data"), 
              type="character", 
              default=NULL, 
              help="Input GTEx data file.", 
              metavar="character"),
  make_option(c("--genes"), 
              type="character", 
              default=NULL, 
              help="Input file with the genes of interest.", 
              metavar="character"),
  make_option(c("-o", "--outdir"),
              type="character",
              default=getwd(),
              help="Output directory. Default = current working directory [%default].",
              metavar="character"),
  make_option(c("--genecol"),
              type="character",
              default="ENSEMBL",
              help="Gene ID column name. Default=ENSEMBL",
              metavar="character"),
  make_option(c("--metric"), 
              type="character", 
              default="rho", 
              help="Method used to compute the gene pair associations. Default=rho. Options=[rho,phi,phs,cor,vlr]", 
              metavar="character"),
  make_option(c("--norm"), 
              type="character", 
              default="NA", 
              help="Data normalization. Default=NA. Options=[rpkm,tmm,NA]", 
              metavar="character"),
  make_option(c("--ivar"), 
              type="character", 
              default="clr", 
              help="Data transformation for propr. Default=clr(centered log-ratio trasformation). Options=[clr,log2,NA]", 
              metavar="character"),
  make_option(c("--tissue"),
              type="character",
              default=NULL,
              help="File with tissue name",
              metavar="character"),
  make_option("--cutoff_interval",
              type="double",
              default=0.005,
              help="Interval for updateCutoff step"),
  make_option("--interval_min",
              type="double",
              default=0.1,
              help="Minimum cutoff for FDR computation"),
  make_option("--interval_max",
              type="double",
              default=0.995,
              help="Maximum cutoff for FDR computation"),
  make_option("--permutation",
              type="integer",
              default=20,
              help="Permutation number for proportionality analysis. Default = %default",
              metavar="number"),
  make_option("--sampseed",
              type="character",
              default=NULL,
              help="Table with the seed for each nsample",
              metavar="character"),
  make_option("--nsamp",
              type="integer",
              default=NULL,
              help="Number of random samples used",
              metavar="number"),
  make_option("--local",
              action="store_true",
              default=FALSE,
              help="Load local libraries.")
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Compute proportionality analysis on a set of genes of interest.")
opt = parse_args(opt_parser)


# ========= #
# LIBRARIES #
# ========= #

if (opt$local){
  # library paths
  .libPaths(c("/nfs/software/R/packages","/nfs/users2/cn/ierb/R/", "/nfs/users2/cn/sjin/R/x86_64-redhat-linux-gnu-library/3.3/"))
}

library(recount)
library(zCompositions)
library(propr)
library(edgeR)


# ===== #
# INPUT #
# ===== #

Sys.time()
print("----reading input files")

# read genes of interest
mygenes_df = read.csv(opt$genes, sep="\t")
mygenes = mygenes_df[,opt$genecol]

# read gtex
load(file.path(opt$data))

# count matrix
m = assay(rse_gene)

# filter low-expression genes
av=apply(m,1,mean)
mygenes_av=rownames(m)[which(av>boxplot(av)$stats[2,1])] 
rse_gene = rse_gene[mygenes_av,] 
m = assay(rse_gene)


# normalize the counts if required
# https://www.reneshbedre.com/blog/expression_units.html
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
Sys.time()
print("----reading and processing (normalize if required) gene expression matrix")
if (opt$norm == "rpkm"){
  m = getRPKM(rse_gene)
}else if (opt$norm == "tmm"){
  # Note that normalization in edgeR is model-based, and the original read counts are not themselves transformed. This means that users should not transform the read counts in any way
  # before inputing them to edgeR. For example, users should not enter RPKM or FPKM values to edgeR in place of read counts. Such quantities will prevent edgeR from correctly
  # estimating the mean-variance relationship in the data, which is a crucial to the statistical
  # strategies underlying edgeR. Similarly, users should not add artificial values to the counts
  # before inputing them to edgeR

  # For further information check: https://www.biostars.org/p/84087/ and https://www.biostars.org/p/317701/
  y = DGEList(counts=m)
  y = calcNormFactors(y)
  m = cpm(y)
}

# data transformation
# https://www.biostars.org/p/100926/
# https://www.biostars.org/p/242573/
# https://blog.qbaseplus.com/seven-tips-for-bio-statistical-analysis-of-gene-expression-data
Sys.time()
print(paste("----transforming data (ivar=", opt$ivar, ")", sep=""))
if (opt$ivar == "clr"){
  M = t(m)
  M = cmultRepl(M,method="CZM",label=0,output="p-counts")  # replace zeros
  M = propr:::proprCLR(M)  # get clr
  m = t(M)
}else if (opt$ivar == 'log2'){
  m = log2(m+1)
}
ivar = NA

# genes and samples
rownames(m)=substring(rownames(m),1,15)
mysamples = seq(1:ncol(m))


# filter samples per tissue
if (!is.null(opt$tissue)){

  # get tissue name
  fil = opt$tissue 
  tissue = readLines(fil) 

  # get tissue samples
  mysamples=which(colData(rse_gene)[,"smtsd"]==tissue)
}

# random sample
if (is.integer(opt$nsamp)){
  print(paste("----getting ", opt$nsamp, " random samples from ", length(mysamples), sep=""))
  samp_df = read.csv(opt$sampseed, sep="\t")
  sampseed = samp_df[which(samp_df[,"nsamp"]==opt$nsamp),"seed"]
  set.seed(sampseed); mysamples=sample(mysamples, opt$nsamp)
  print(mysamples)
}

# filter matrix
Sys.time()
print(paste("----subsetting gene expression matrix with [samples:", length(mysamples), "][genes:", length(mygenes),"]", sep=""))
M=t(m[mygenes,mysamples])

# remove variables
rm(rse_gene)
rm(av)
rm(m)
rm(mygenes)
rm(mysamples)


# =============== #
# PROPORTIONALITY #
# =============== #

Sys.time()
print(paste("----running proportionality analysis on matrix size of [", nrow(M), "][", ncol(M), "]", sep=""))

# calculate proportionality 
pro=propr(M, metric=opt$metric, ivar=ivar, p=opt$permutation)  # use p=20 for faster computation. Moreover, p is only useful for the calculation of FDR, and now we don't care that much on the FDR.
pro=updateCutoffs(pro,seq(opt$interval_min, opt$interval_max, opt$cutoff_interval))
# pro = updateCutoffs(pro, c(seq(0.1, 0.7, 0.01), seq(0.705, 0.995, 0.005)))
# pro = updateCutoffs(pro, seq(0.1, 0.995, 0.005))

rm(M)

# ====== #
# OUTPUT #
# ====== #

Sys.time()
print("----writing output files")

# create output directory
if (!dir.exists(opt$outdir)) {dir.create(opt$outdir, recursive=TRUE)}

# write r data
rdata = file.path(opt$outdir, "enzyme_results.Rdata")
save(pro, file=rdata)

# write matrix
matout = file.path(opt$outdir, "enzyme_results.mat")
write.table(pro@matrix, matout)

# select proportional pairs above smallest cutoff
mypairs=which(pro@results[,"propr"]>=min(pro@fdr[, "cutoff"]))
results = pro@results[mypairs,]

# write pairs
out1 = file.path(opt$outdir, "enzyme_results.csv")
write.table(results, out1, row.names = FALSE, quote=FALSE, sep=",", dec=".")

# write cutoff
cutoff = min(pro@fdr[which(pro@fdr[,"FDR"]<0.05), "cutoff"])
out2 = file.path(opt$outdir, "enzyme_cutoff.txt")
writeLines(as.character(cutoff), out2)

# write FDR
out3 = file.path(opt$outdir, "enzyme_fdr.csv")
write.table(pro@fdr, out3, row.names = FALSE, quote=FALSE, sep=",", dec=".")
