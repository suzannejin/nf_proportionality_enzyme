# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option(c("-i", "--input"), 
              type="character", 
              default=NULL, 
              help="Input proportionality results.", 
              metavar="character"),
  make_option(c("-o", "--outdir"),
             type="character",
             default=getwd(),
             help="Output dir",
             metavar="character"),
  make_option(c("-g","--mygenes"),
              type="character",
              default=NULL,
              help="Genes names list",
              metavar="character"),
  make_option(c("-c", "--cutoff"),
              type="double",
              default=NULL,
              help="Proportionality cutoff"),
  make_option(c("-p","--permutation"),
              type="integer",
              default=100,
              help="Number permutation",
              metavar="number"),
  make_option(c("--filter"),
              action="store_true",
              default=FALSE,
              help="Filter according to GO ontology and evidence"),
  make_option("--local",
              action="store_true",
              default=FALSE,
              help="Load local libraries.")
); 
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

if (opt$local){
  # library paths
  .libPaths(c("/nfs/software/R/packages","/nfs/users2/cn/ierb/R/", "/nfs/users2/cn/sjin/R/x86_64-redhat-linux-gnu-library/3.3/"))
}

library(data.table)
library(igraph)
library(graflex)
library(GO.db)
library(org.Hs.eg.db)


# ========= #
# FUNCTIONS #
# ========= #

get_adj <- function(ppairs,mygenes,idx,ngenes_1){
  
  # initialize matrix
  # note that the matrix size is equal to the original number of genes used for propr step (2223) x (2223)
  mat = matrix(0,nrow=ngenes_1,ncol=ngenes_1)  
  
  # update proportional pairs
  for (i in 1:nrow(ppairs)){
    a=as.integer(ppairs[i, "Partner"]) 
    b=as.integer(ppairs[i, "Pair"])
    mat[a,b] = 1
    mat[b,a] = 1
  }

  # matrix with the filtered (if required) genes
  mat = mat[idx,idx] 
  
  # rename
  rownames(mat) = mygenes
  colnames(mat) = mygenes
  
  return(mat)
}



annotate_go_term <- function(df, go="Concept"){
    terms = c()
    for (i in 1:nrow(df)){

        # GO concept
        concept=df[i, go]

        # GO term
        term = as.character(select(GO.db, keys=concept, columns="TERM", keytype="GOID")["TERM"][1])
        terms = c(terms, term)
    }
    df["Term"] = terms
    
    return(df)
}

getK2 <- function (keys, columns = "GO", keytype = "ENSEMBL", minK = 10, ont=NA, evi=NA) 
{
    if (length(columns) > 1) 
        stop("Please provide a single column.")
    if (!is.na(ont)){
        columns = c(columns, "ONTOLOGY")
    }
    if (!is.na(evi)){
        columns = c(columns, "EVIDENCE")
    }
    # check packages
    library("AnnotationDbi")
    library("org.Hs.eg.db")
    db <- org.Hs.eg.db::org.Hs.eg.db
    
    # get pathways-gene data frame
    godf <- AnnotationDbi::select(db, keys = keys, columns = columns, 
        keytype = keytype)
    
    # filter if required
    if (!is.na(ont)){
        godf = godf[which(godf$ONTOLOGY %in% ont),]
    }
    if (!is.na(evi)){
        godf = godf[which(godf$EVIDENCE %in% evi),]
    }
    # reorganize data frame
    gotab <- table(godf[, 1:2])
    gotab[gotab > 1] <- 1
#     gotab <- gotab[keys, ]
    if (length(rownames(gotab)) > length(which(rownames(gotab) %in% keys))) {
        stop("Uh oh! Unexpected mapping.")
    }
    gotab[, colSums(gotab) >= minK]
}


# ===== #
# INPUT #
# ===== #

Sys.time()
print("----reading input files")

# read proportional pairs
pairs = fread(opt$input)
ppairs = pairs[which(pairs[,"propr"]>=opt$cutoff),]

# read gene names
mygenes_df = read.csv(opt$mygenes, sep="\t")
mygenes = mygenes_df$ENSEMBL
mygenes_df$idx = 1:nrow(mygenes_df)
ngenes_1 = nrow(mygenes_df)  # number of genes used for the propr step


# ======= #
# GRAFLEX #
# ======= #

Sys.time()
print("----get GO reference graph")

# gene-pathway matrix
if (opt$filter){  ## filtered by GO ontology and evidence
  K = getK2(mygenes, ont=c("BP"), evi=c("EXP", "IDA", "IPI", "IMP", "IGI", "TAS", "IC"))
  mygenes_df = mygenes_df[which(mygenes_df$ENSEMBL %in% rownames(K)),]
  mygenes = mygenes_df$ENSEMBL
}else{
  K = getK(mygenes)
}

# adj matrix
A = get_adj(ppairs, mygenes, mygenes_df$idx, ngenes_1)

Sys.time()
print(paste("----GREA using graflex  [K:", nrow(K), "x", ncol(K), "] [A:", nrow(A), "x", ncol(A), "]", sep=""))

# graflex
or = graflex(A, K, p=opt$permutation)

# add GO term description
or = annotate_go_term(or)

# write output
Sys.time()
print("----writing output files")
out = paste(opt$outdir, "/", "graflex_enzyme_c", opt$cutoff*100, "_p", opt$permutation, ".tsv", sep="")
write.table(or, out, row.names = FALSE, quote=FALSE, sep="\t", dec=".")
