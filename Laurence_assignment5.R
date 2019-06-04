# install standard R packages:
install.packages("ggplot2") # for plotting
install.packages("reshape") # for Venn diagrams
install.packages("Vennerable", repos="http://R-Forge.R-project.org") # for
# Venn diagrams
# vector of Bioconductor package names:
pkgs <- c(
  "rtracklayer",
  "TFBSTools",
  "JASPAR2014",
  "BSgenome.Hsapiens.UCSC.hg18",
  "TxDb.Hsapiens.UCSC.hg18.knownGene"
)
# install BioConductor packages:
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs, suppressUpdates=TRUE)

# set two lists with the names of the Biocodnuctor and the CRAN packages
bioc.libs <- c("Biobase", "affy", "annotate", "hgu95av2.db", "hgu95av2cdf",
               "estrogen", "limma", "genefilter", "vsn", "GEOquery", "annaffy",
               "GO.db","KEGG.db")
cran.libs <- c("RColorBrewer", "R.utils", "gplots")

source("http://bioconductor.org/biocLite.R") # get BiocInstaller
biocLite() # if they are outdated packages, choose the option to update all
# the packages

library(BiocInstaller)
biocLite(bioc.libs)
biocLite(cran.libs)

library(BiocInstaller)
biocValid()

source("http://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")

FOXA1.df <- read.table("http://www.carroll-lab.org.uk/FreshFiles/Data/Data_Sheet_3/MCF7_FOXA1%20binding.bed",header=TRUE)
# show the first six lines of the data.frame
head(FOXA1.df)

require(GenomicRanges) # load the GenomicRanges pakcage
FOXA1 <- GRanges(
  FOXA1.df$chr,
  IRanges(FOXA1.df$star, FOXA1.df$end),
  strand="*"
)
# we can add more data to each peak subsequently
names(FOXA1) <- paste("FOXA1_peak", 1:nrow(FOXA1.df), sep="_")
score(FOXA1) <- FOXA1.df[,"X.10.log10.pvalue."]

# show the first and last lines of the GRanges object
FOXA1

length(FOXA1)
mean(width(FOXA1))
median(width(FOXA1))
max(width(FOXA1))
FOXA1 <- FOXA1[width(FOXA1) <= 2000]
hist(width(FOXA1))

# get the -log_10 transformed p-values from the score column
mlPvals <- score(FOXA1)
hist(mlPvals, xlab="-log_10(p-value)", col="gray")

require(rtracklayer) # load the rtracklayer package
ER <- import.bed("http://www.carroll-lab.org.uk/FreshFiles/Data/Data_Sheet_3/MCF7_ER_binding.bed")
# assign scores by converting the 'name' feald to type numeric
score(ER) <- as.numeric(ER$name)
# overwrite the name column
ER$name <- paste("ER_peaks", 1:length(ER), sep="_")
# use the names() function
names(ER) <- ER$name
ER
bp <- barplot(c(length(ER), length(FOXA1)), names=c("ER", "FOXA1"))
# add actual values as text lables to the plot
text(bp, c(length(ER), length(FOXA1)), labels=c(length(ER),
                                                length(FOXA1)), pos=1)
findOverlaps(FOXA1, ER)

# get subsets of binding sites
ER.uniq <- setdiff(ER, FOXA1)
FOXA1.uniq <- setdiff(FOXA1, ER)
# plot overlap of binding sites as Venn diagram
install.packages("Vennerable")
require(Vennerable)
# build objects with the numbers of sites in the subsets
venn <- Venn(SetNames=c("ER", "FOXA1"),
             Weight=c(
               '10'=length(ER.uniq),
               '11'=length(ovl),
               '01'=length(FOXA1.uniq)
             )
)
# plot Venn Diagram
plot(venn)

# load gene annotation from UCSC for human genome build hg18
require(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene # just a shortcut

require(ChIPpeakAnno) # laod package to annotate peaks with genes
# calculate the overlap with features
ER.features <- assignChromosomeRegion(ER, TxDb=txdb, nucleotideLevel=FALSE)
# show the results
ER.features

# make pie-chart
pie(ER.features$percentage)# make pie-chart
# make barplot
bp <- barplot(ER.features$percentage, ylab="%")
text(bp, ER.features$percentage, signif(ER.features$percentage, 4),
     pos=1)
# get all genes from the data base as GRanges object
genes <- genes(txdb)
# take the region around the gene start as promoter
prom <- promoters(genes, upstream=2000, downstream=200)
prom
# get only those ER peaks that overlap a promoter region
ERatProm <- subsetByOverlaps(ER, prom)
# subset size
length(ERatProm)
## [1] 684
# percent of all ER peaks
length(ERatProm) / length(ER) * 100
# We search for overlap between ER peaks and promoters
ERatProm.Hits = findOverlaps(ER, prom)
ERprom = genes[subjectHits(ERatProm.Hits)]
# take only unique ids
gene.ids <- unique(names(ERprom))
# write names to an output file
write.table(gene.ids, file="ER_regulated_genes.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
# get TSS as only the start coordinate of genes.
tss <- resize(genes, width=1, fix="start")
# calculate the distances from peaks to tss
d <- distanceToNearest(ER, tss)
# show object d
d
# show the metacolumns as DataFrame object
mcols(d)
# extract the distance column as vector
dist <- mcols(d)[,1]
# get the average distance in kb
mean(dist) * 10^-3
# subset hits object by distance
close.Hits <- d[dist <= 10000,]
# show the subset
close.Hits

# get the indeces of genes
ER.genes <- genes[subjectHits(close.Hits)]
# extract the vector of names from the GRanges object
gene.ids <- names(ER.genes)
# take only unique ids
gene.ids <- unique(gene.ids)
# write names to an output file
write.table(gene.ids, file="ER_regulated_genes.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
# extract the vector of names from the GRanges object
gene.ids <- names(ER.genes)
# take only unique ids
gene.ids <- unique(gene.ids)
# write names to an output file
write.table(gene.ids, file="ER_regulated_genes.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
export.bed(ER, "ER_peaks.bed")
export.bed(FOXA1, "FOXA1_peaks.bed")

"7157" %in% gene.ids

require(rtracklayer)
ER <- import.bed("ER_peaks.bed")
FOXA1 <- import.bed("FOXA1_peaks.bed")
# take only the subset of 100 peaks with the highest scores
ER.top <- ER[order(ER$score, decreasing = TRUE)][1:100]
FOXA1.top <- FOXA1[order(FOXA1$score, decreasing = TRUE)][1:100]
require(BSgenome.Hsapiens.UCSC.hg18) # human reference genome
# get DNA sequence of peaks (this might take some time)
ER.seq <- getSeq(Hsapiens, ER.top)
ER.seq
FOXA1.seq <- getSeq(Hsapiens, FOXA1.top)
# write sequences to fasta file
writeXStringSet(ER.seq, "./ER_seq.fasta")
writeXStringSet(FOXA1.seq, "./FOXA1_seq.fasta")

# get ER motif from JASPAR
require(JASPAR2014) # load the JASPAR2014 package
require(TFBSTools)
# search for "ESR1" by using species human (9606) and take the first result
pfm <- getMatrixSet(JASPAR2014, list(species=9606, name="ESR1"))[[1]]
pfm
# we convert the matrix to the information content using the toICM
# function.
icm <- toICM(pfm)
# plot as sequence logo plot
seqLogo(icm)
matchPWM(as.matrix(pfm), ER.seq[[3]])
hitsER <- lapply(ER.seq, function(s) matchPWM(as.matrix(pfm), s, min.score="75%") )
hitsFOXA1 <- lapply(FOXA1.seq, function(s) matchPWM(as.matrix(pfm), s,
min.score="75%") )
countsER <- sapply(hitsER, length)
countsFOXA1 <- sapply(hitsFOXA1, length)
sum(countsER >= 1)
# create a data.frame with the number of motif hits and the corresponding peak
# source group as columns.
plotDF <- data.frame(
  motif_hits = c(countsER, countsFOXA1),
  group = rep(c("ER peaks", "FOXA1 peaks"), c(length(countsER),
                                              length(countsFOXA1)))
)
# load the ggplot2 package
require(ggplot2)
qplot(x=motif_hits, fill=group, data=plotDF, geom="bar", facets=group~.)
