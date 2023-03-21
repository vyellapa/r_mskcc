gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
           join2 = c(org.Hs.eg.db="ENTREZID",
                     TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID"))

require("Homo.sapiens")

test_extractPkgsAndCols <- function(){
  gdm <- OrganismDbi:::.mungeGraphData(gd)
  res <- OrganismDbi:::.extractPkgsAndCols(gdm)
  checkEquals(names(res), c("GO.db","org.Hs.eg.db","org.Hs.eg.db",
                          "TxDb.Hsapiens.UCSC.hg19.knownGene"))
  checkTrue(all(res %in%  c("GOID","ENTREZID","GO","GENEID")))
}

test_testKeys <- function(){
  gdm <- OrganismDbi:::.mungeGraphData(gd)
  res <- OrganismDbi:::.extractPkgsAndCols(gdm)
  ## the following should not blow up
  OrganismDbi:::.testKeys(res)

  ## the following should blow up
  res2 <- c("GO.db","org.Hs.eg.db","org.Hs.eg.db",
            "TxDb.Hsapiens.UCSC.hg19.knownGene")
  names(res2) <- c("WRONGVALUE","ENTREZID","GO","GENEID")
  checkException(OrganismDbi:::.testKeys(res2))

  ## note that the case of a bad DB would have been caught elsewhere
  ## (and earlier)
}

## test_TxDb <- function(){
##     require("Homo.sapiens")
##     txdbMan <- loadDb(system.file('extdata','myTxDb.sqlite',
##                                     package='OrganismDbi'))
##     ## I have to be explicit here.
##     txdbFull <- loadDb(system.file('extdata',
##                                    'TxDb.Hsapiens.UCSC.hg19.knownGene.sqlite',
##                                    package='TxDb.Hsapiens.UCSC.hg19.knownGene'))
##     odb <- makeOrganismDbFromTxDb(txdb=txdbMan)
##     res <- resources(odb)
##     checkTrue(res[['TxDb.Hsapiens.UCSC.hg19.knownGene']] == dbfile(txdbMan))
##     m1 <- metadata(TxDb(odb))
##     checkTrue(as.integer(m1[m1$name=='transcript_nrow','value']) < 10)
    
##     ## NOW change it:
##     TxDb(odb) <- txdbFull
##     res <- resources(odb)
##     checkTrue(res[['TxDb.Hsapiens.UCSC.hg19.knownGene']] == dbfile(txdbFull))
##     m2 <- metadata(TxDb(odb))
##     checkTrue(as.integer(m2[m2$name=='transcript_nrow','value']) > 82000)

##     ## check that a TxDb made on the fly (stored in local RAM) also works OK
##     transcript_ids <- c(
##                         "uc001aaa.3",
##                         "uc010nxq.1",
##                         "uc010nxr.1",
##                         "uc001aal.1",
##                         "uc001aaq.2"
##                         )
##     txdbInstant <- makeTxDbFromUCSC(genome="hg19", tablename="knownGene",
##                                     transcript_ids=transcript_ids)
##     TxDb(odb) <- txdbInstant
##     res <- resources(odb)
##     checkTrue(res[['TxDb.Hsapiens.UCSC.hg19.knownGene']] == dbfile(txdbInstant))
##     checkTrue(res[['TxDb.Hsapiens.UCSC.hg19.knownGene']] == "")
##     m3 <- metadata(TxDb(odb))
##     checkTrue(as.integer(m3[m3$name=='transcript_nrow','value']) < 10)

##     ## NOW (CRITICALLY) Set things back to the other TxDb (for both)
##     TxDb(odb) <- txdbFull
##     TxDb.Hsapiens.UCSC.hg19.knownGene <- TxDb(odb)
## }

## Fast testing (AWESOME):
## BiocGenerics:::testPackage(pattern="^test_construction.*\\.R$")


## Slower testing that requires installation (Not so awesome):
## BiocGenerics:::testPackage("OrganismDbi", pattern="^test_construction.*\\.R$")

## TODO: fix the bug that is changing the value of the name (assigning to the original symbols name inside the function)  And then uncomment this test...
