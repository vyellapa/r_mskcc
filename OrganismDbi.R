### R code from vignette source 'OrganismDbi.Rnw'

###################################################
### code chunk number 1: columns
###################################################
library(Homo.sapiens)
columns(Homo.sapiens)


###################################################
### code chunk number 2: keys
###################################################
keytypes(Homo.sapiens)


###################################################
### code chunk number 3: keys
###################################################
head(keys(Homo.sapiens, keytype="ENTREZID"))


###################################################
### code chunk number 4: select
###################################################
k <- head(keys(Homo.sapiens, keytype="ENTREZID"),n=3)
select(Homo.sapiens, keys=k, columns=c("TXNAME","SYMBOL"), keytype="ENTREZID")


###################################################
### code chunk number 5: transcripts
###################################################
transcripts(Homo.sapiens, columns=c("TXNAME","SYMBOL"))


###################################################
### code chunk number 6: transcriptsBy
###################################################
transcriptsBy(Homo.sapiens, by="gene", columns=c("TXNAME","SYMBOL"))


###################################################
### code chunk number 7: setupColData
###################################################
gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
           join2 = c(org.Hs.eg.db="ENTREZID",
                     TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID"))


###################################################
### code chunk number 8: makeOrganismPackage (eval = FALSE)
###################################################
## destination <- tempfile()
## dir.create(destination)
## makeOrganismPackage(pkgname = "Homo.sapiens",
##   graphData = gd,
##   organism = "Homo sapiens",
##   version = "1.0.0",
##   maintainer = "Package Maintainer<maintainer@somewhere.org>",
##   author = "Some Body",
##   destDir = destination,
##   license = "Artistic-2.0")


