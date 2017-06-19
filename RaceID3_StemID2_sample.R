## install required packages (only at first time)
## install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit","vegan","Rtsne","scran","randomForest","rgl"))
## source("https://bioconductor.org/biocLite.R")
## biocLite("scran")
## biocLite("DESeq2")
## biocLite("biomaRt")

## !!!! novel:
## normalization to minimum instead of median as alternative to downsampling, sum factor normalization (sfn), house keeping gene normalization (sfn), improved faster clustering (k-medoids), new possible metrics (e. g. logpearson), feature selection (FSelect), novel plotexptsne (highest expression on top), comptsne (fast,perplexity), cdiff with negative binomial, fast version of StemID: projback, lineagetree, comppvalue (fast), plot StemID score and lineage trees with linkscore threshold (scthr), filter out single genes and correlating groups of genes for clustering, random forest based (OOB) correction, PCA-based cell cycle correction, correcting batch effect by regression, FateID + Pseudotemporal ordering


## !!!

## load class definition and functions

source("RaceID3_StemID2_class.R")

## input data
x <- read.csv("transcript_counts_intestine_5days_YFP.xls",sep="\t",header=TRUE)
rownames(x) <- x$GENEID

# prdata: data.frame with transcript counts for all genes (rows) in all cells (columns); with rownames == gene ids; remove ERCC spike-ins 
prdata <- x[grep("ERCC",rownames(x),invert=TRUE),-1]

## RaceID3
# initialize SCseq object with transcript counts
sc <- SCseq(prdata)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=FALSE, sfn=FALSE, hkn=FALSE, dsn=1, rseed=17000, CGenes=NULL, FGenes=NULL, ccor=.4)

# regress out the batch effect
# optional:
#vars <- data.frame(row.names=names(sc@fdata),batch=sub("(_|)\\d.+","",names(sc@fdata)))
#sc@fdata <- varRegression(sc@fdata,vars)

# correct for cell cycle, proliferation, and expression of degradation markers by PCA
# optional:
require(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
g   <- sub("__chr.+","",rownames(sc@fdata));
k   <- getBM(attributes = c("external_gene_name", "go_id","name_1006"),filters="external_gene_name",values=g,mart=mart)
gCC <- name2id( k$external_gene_name[k$name_1006 == "cell cycle"],rownames(sc@fdata)) 
gCP <- name2id( k$external_gene_name[k$name_1006 == "cell proliferation"],rownames(sc@fdata))
vset <- list(gCC,gCP)
x <- CCcorrect(sc@fdata,vset=vset,CGenes=NULL,ccor=.4,nComp=NULL,pvalue=.05,quant=.01,mode="pca")
# number of principal components that have been removed
x$n
# loadings of the first principal component that has been removed
y <- x$pca$rotation[,x$n[1]]
# genes from vset are either enriched in the head or the tail of this list
tail(y[order(y,decreasing=TRUE)],10)
# reassign the corrected expression matrix to sc@fdata
sc@fdata <- x$xcor

# k-medoids clustering
sc <- clustexp(sc,clustnr=30,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids",FSelect=TRUE)
# compute t-SNE map
sc <- comptsne(sc,rseed=15555,sammonmap=FALSE,initial_cmd=TRUE,fast=TRUE,perplexity=30)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.95)
# reassign clusters based on random forest
sc <- rfcorrect(sc,rfseed=12345,final=TRUE,nbfactor=5)

## diagnostic plots
# gap statistics: only if do.gap == TRUE
##plotgap(sc)
# plot within-cluster dispersion as a function of the cluster number: only if sat == TRUE
plotsaturation(sc,disp=TRUE)
# plot change of the within-cluster dispersion as a function of the cluster number: only if sat == TRUE
plotsaturation(sc)
# silhouette of k-medoids clusters
plotsilhouette(sc)
# Jaccard's similarity of k-medoids clusters
plotjaccard(sc)
# barchart of outlier probabilities
plotoutlierprobs(sc)
# regression of background model
plotbackground(sc)
# dependence of outlier number on probability threshold (probthr)
plotsensitivity(sc)
# heatmap of k-medoids cluster
clustheatmap(sc,final=FALSE,hmethod="single")
# heatmap of final cluster
clustheatmap(sc,final=TRUE,hmethod="single")
# highlight k-medoids clusters in t-SNE map
plottsne(sc,final=FALSE)
# highlight final clusters in t-SNE map
plottsne(sc,final=TRUE)
# highlight cell labels in t-SNE map
plotlabelstsne(sc,labels=sub("(\\_\\d+)","",names(sc@ndata)))
# highlight groups of cells by symbols in t-SNE map
plotsymbolstsne(sc,types=sub("(\\_\\d+)$","", names(sc@ndata)))
# highlight transcirpt counts of a set of genes in t-SNE map, e. g. all Apoa genes
g <- c("Apoa1__chr9", "Apoa1bp__chr3", "Apoa2__chr1", "Apoa4__chr9", "Apoa5__chr9")
plotexptsne(sc,g,n="Apoa genes",logsc=TRUE)

## identification of marker genes
# differentially regulated genes in each cluster compared to the full ensemble
cdiff <- clustdiffgenes(sc,pvalue=.01)

## write results to text files
# final clusters 
x <- data.frame(CELLID=names(sc@cpart),cluster=sc@cpart)
write.table(x[order(x$cluster,decreasing=FALSE),],"cell_clust.xls",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
# differentially expressed genes in cluster
for ( n in names(cdiff) ) write.table(data.frame(GENEID=rownames(cdiff[[n]]),cdiff[[n]]),paste(paste("cell_clust_diff_genes",sub("\\.","\\_",n),sep="_"),".xls",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

# differentially expressed genes between two sets of clusters, e. g. cluster 3 and clusters 2,11
d <- diffgenes(sc,cl1=3,cl2=c(2,11),mincount=1)
plotdiffgenes(d,gene=names(d$z)[1])

# differentially expressed genes between two sets of cells, e. g. cluster 3 and clusters 2,11 by a different approach based on the overlap of transctipt count distributions
A <- names(sc@cpart)[sc@cpart == 3]
B <- names(sc@cpart)[sc@cpart %in% c(2,11)]
x <- diffexpnb(sc@ndata, A=A, B=B, method="pooled",norm=FALSE, DESeq=FALSE, vfit=sc@background$vfit, locreg=FALSE)
plotdiffgenesnb(x,pthr=.05,lthr=1,mthr=1,Aname="Cl.3",Bname="Cl.2,11",show_names=TRUE,padj=TRUE)

# plot 3D tSNE
plot3dtsne(sc,perplexity=30,fast=TRUE,x=NULL,g=NULL,logsc=FALSE,final=TRUE,ret=FALSE,tp=1)



## StemID2

# initialization
ltr <- Ltree(sc)
# computation of the entropy
ltr <- compentropy(ltr)
# computation of the projections for all cells
ltr <- projcells(ltr,cthr=2,nmode=FALSE)
# computation of the projections for all cells after randomization
ltr <- projback(ltr,pdishuf=2000,nmode=FALSE,fast=FALSE,rseed=17000)
# assembly of the lineage tree
ltr <- lineagetree(ltr,pthr=0.05,nmode=FALSE,fast=FALSE)
# determination of significant differentiation trajectories
ltr <- comppvalue(ltr,pethr=0.05,nmode=FALSE,fast=FALSE)

## diagnostic plots
# histogram of ratio between cell-to-cell distances in the embedded and the input space
plotdistanceratio(ltr)
# t-SNE map of the clusters with more than cthr cells including a minimum spanning tree for the cluster medoids
plotmap(ltr)
# visualization of the projections in t-SNE space overlayed with a minimum spanning tree connecting the cluster medoids
plotmapprojections(ltr)
# lineage tree showing the projections of all cells in t-SNE space
plottree(ltr,showCells=TRUE,nmode=FALSE,scthr=.3)
# lineage tree without showing the projections of all cells
plottree(ltr,showCells=FALSE,nmode=FALSE,scthr=.3)
# heatmap of the enrichment p-values for all inter-cluster links
plotlinkpv(ltr)
# heatmap of the link score for all inter-cluster links
plotlinkscore(ltr)
# heatmap showing the fold enrichment (or depletion) for significantly enriched or depleted links
projenrichment(ltr)

## extract projections onto all links for all cells in a given cluster i
x <- getproj(ltr,i=3)
# heatmap of all projections for cluster i
pheatmap(x$pr)
# heatmap of z-score for all projections for cluster i
pheatmap(x$prz)

## extracting all cells on two branches sharing the same cluster and computing differentially expressed genes between these two branches
x <- branchcells(ltr,list("3.7","2.3"))
# z-scores for differentially expressed genes
head(x$diffgenes$z)
# plotting the cells on the two branches as additional clusters in the t-SNE map
plottsne(x$scl)


## computing the StemID2 score
x <- compscore(ltr,nn=1,scthr=0)
#plotting the StemID2 score
plotscore(ltr,nn=1,scthr=0)

# retrieve cells from branch in pseudo-temporal order as inferred by the projection coordinates
n <- cellsfromtree(ltr,c(3,2,11,4))

# filter out lowly expressed genes
fs  <- filterset(ltr@sc@ndata,n=n$f,minexpr=2,minnumber=1)
# compute self organizing map (SOM) of co-expressed genes
s1d <- getsom(fs,nb=1000,k=5,locreg=TRUE,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)

# coloring scheme for clusters ( vector with colours)
fcol <- ltr@sc@fcol
y    <- ltr@sc@cpart[n$f]
# plot average z-score for all modules derived from the SOM
plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
# plot z-score profile of each gene
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
# plot normalized expression profile of each gene
plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
# plot binarized expression profile of each gene (z-score < -1, -1 < z-score < 1, z-score > 1)
plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

# extract all genes from node 1 of the SOM
g <- names(ps$nodes)[ps$nodes == 1]
# plot average expression profile of these genes along the trajectory
plotexpression(fs,y,g,n$f,k=25,col=fcol,name="Node 1",cluster=FALSE,locreg=TRUE,alpha=.5,types=NULL)
# plot expression of a single gene
plotexpression(fs,y,"Clca4__chr3",n$f,k=25,col=fcol,cluster=FALSE,locreg=TRUE,alpha=.5,types=NULL)
# plot average expression profile of these genes along the trajectory, highlighting batch origin
plotexpression(fs,y,g,n$f,k=25,col=fcol,name="Node 1",cluster=FALSE,locreg=TRUE,alpha=.5,types=sub("\\_\\d+","",n$f))








