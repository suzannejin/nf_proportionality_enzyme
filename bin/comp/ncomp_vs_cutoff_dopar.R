# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option("--input1", 
              type="character", 
              default=NULL, 
              help="Input file 1: propr pairs list", 
              metavar="character"),
  make_option("--input2",
              type="character",
              default=NULL,
              help="Input file 2: cutoff-FDR list",
              metavar="character"),
  make_option("--outdir",
              type="character",
              default=NULL,
              help="Where the output files will be stored",
              metavar="character"),
  make_option("--prefix",
              type="character",
              default=NULL, 
              help="File prefix",
              metavar="character"),
  make_option("--title",
              type="character",
              default=NULL,
              help="Plot title",
              metavar="character"),
  make_option("--ncores",
              type="integer",
              default=1,
              help="Number of cpus",
              metavar="number")
); 
 
opt_parser = OptionParser(usage="Plot number of component vs proportionality cutoff.", option_list=option_list)
opt = parse_args(opt_parser)

# some variabless
prefix = opt$prefix
outdir = opt$outdir
if (!dir.exists(outdir)) {dir.create(outdir, recursive=TRUE)}
if (is.null(opt$title)){
    title = prefix
}else{
    title = opt$title
}


# ======= #
# LIBRARY #
# ======= #

library(igraph)
library(ggplot2)
library(doMC)


# ======= #
# DATASET #
# ======= #

# Read input files
print("----loading files")
file1 = read.csv(opt$input1)
file2 = read.csv(opt$input2)

# get cutoff list
l = file2[,c("cutoff","FDR")]
l[,"FDR"] = round( l[,"FDR"]*100, 2)
rownames(l) = 1:nrow(l)
l["ncom"] = 0; l["between20-100"]=0; l["ngenes"]=0; l["npairs"]=0
for (i in 1:2){
    l[paste("comsize", i, sep="")] = 0
}

# list 2 with the size of all components
l2 = list()
l2$cutoff = c()
l2$FDR = c()
l2$comsize = c()

# parallelize computation for i cutoff
doMC::registerDoMC(cores = opt$ncores)
`%dopar%` <- foreach::`%dopar%`
x <- foreach::foreach(i = 1:nrow(l)) %dopar% {

    cutoff = l[i, "cutoff"]
    print(paste("compute component for cutoff ", cutoff, sep=""))

    # get proportional pairs
    mypairs = which(file1[,"propr"]>=cutoff)
    ppairs = file1[mypairs,1:2]
    npairs = nrow(ppairs)

    # compute graph
    graph=graph_from_data_frame(ppairs,directed=FALSE)

    # compute components
    com = components(graph)

    # return components
    list(com, npairs)
}

# update values
print("----updating values")
for (i in 1:nrow(l)){

    # components
    com = x[[i]][[1]]

    # update list1
    l[i, "ncom"] = com$no
    l[i, "between20-100"] = length(com$csize[com$csize>=20 & com$csize<=100])
    for (j in 1:2){
        if (com$no>=j){
            l[i, paste("comsize", j, sep="")] = sort(com$csize,decreasing=TRUE)[j]
        }
    }
    l[i, "ngenes"] = length(com$membership)
    l[i, "npairs"] = x[[i]][[2]]

    # update list2
    l2$cutoff = c(l2$cutoff, rep(l[i, "cutoff"], com$no))
    l2$FDR = c(l2$FDR, rep(l[i, "FDR"], com$no))
    l2$comsize = c(l2$comsize, com$csize)
}
l2 = data.frame(l2)

# write list1
out = paste(outdir, "/ncomp_list.txt", sep="")
write.table(l, out, row.names = FALSE, quote=FALSE, sep=",", dec=".")

# write list2
out2 = paste(outdir, "/ncomp_list2.txt", sep="")
write.table(l2, out2, row.names = FALSE, quote=FALSE, sep=",", dec=".")


# ===== #
# PLOTS #
# ===== #

# plot ncom, comsize1, comsize2, between20-100
gl = data.frame(
        "cutoff" = rep( l[,"cutoff"], 4),
        "value" = c(l[,"ncom"], l[,"comsize1"], l[,"comsize2"], l[,"between20-100"]),
        "type" = c(rep("ncom", nrow(l)), rep("comsize1", nrow(l)), rep("comsize2", nrow(l)), rep("between20-100", nrow(l)))
    )
g = ggplot(gl, aes(x=cutoff, y=value)) + geom_bar(stat="identity") +
                facet_wrap(~type, scales="free_y", ncol=2) + 
                geom_text(aes(label=value), color="black", size=1.5, hjust=0, angle=90) +
                ggtitle(title)
ggsave(paste(outdir, "/ncomp_plot.png", sep=""), width=18, height=10)

# plot ncom, comsize1, ngenes, npairs
gl = data.frame(
        "cutoff" = rep( l[,"cutoff"], 4),
        "value" = c(l[,"ncom"], l[,"comsize1"], l[,"ngenes"], l[,"npairs"]),
        "type" = c(rep("ncom", nrow(l)), rep("comsize1", nrow(l)), rep("ngenes", nrow(l)), rep("npairs", nrow(l)))
    )
g = ggplot(gl, aes(x=cutoff, y=value)) + geom_bar(stat="identity") +
                facet_wrap(~type, scales="free_y", ncol=2) + 
                geom_text(aes(label=value), color="black", size=1.5, hjust=0, angle=90) +
                ggtitle(title)
ggsave(paste(outdir, "/ncomp_plot2.png", sep=""), width=18, height=10)


# plot comsize vs cutoff
g = ggplot(l2, aes(x=cutoff, y=comsize, group=cutoff)) + geom_boxplot() + stat_summary(fun.y="mean", shape=23, fill="white") + ggtitle(title)
ggsave(paste(outdir, "/ncomp_plot_comsize.png", sep=""), width=18, height=10)