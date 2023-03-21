## unit tests for the "meta-select"
## These tests are testing the software and not the indiv. packages
## Will base testing on humans for now
require(RSQLite)
require("Homo.sapiens")
x <- Homo.sapiens
## debug(OrganismDbi:::.select)
## debug(OrganismDbi:::.mergeSelectResults)
require("Rattus.norvegicus") 
r <- Rattus.norvegicus 

test_keytypes <- function(){
  res <- OrganismDbi:::.keytypes(x)
  checkTrue("GOID" %in% res)
  checkTrue("TXID" %in% res)
  checkTrue("ENTREZID" %in% res)
}

test_makekeytypeMapping <- function(){
  res <- OrganismDbi:::.makekeytypeMapping(x)
  checkTrue("GO.db" == names(res)[res=='GOID'] )
  checkTrue("TxDb.Hsapiens.UCSC.hg19.knownGene" == names(res)[res=='TXID'] )
  checkTrue("org.Hs.eg.db" == names(res)[res=='ENTREZID'] )
}

test_lookupDbFromKeytype <- function(){
  res <- OrganismDbi:::.lookupDbFromKeytype(x, "GOID")
  checkTrue(class(res)=="GODb")  
  res <- OrganismDbi:::.lookupDbFromKeytype(x, "TXID")
  checkTrue(class(res)=="TxDb")
  res <- OrganismDbi:::.lookupDbFromKeytype(x, "ENTREZID")
  checkTrue(class(res)=="OrgDb")
}

test_lookupDbFromKeytype2 <- function(){
  res <- OrganismDbi:::.keys(x, "GOID")
  checkTrue(is.character(head(res)))
  checkTrue(length(res) > 33000)

  res <- OrganismDbi:::.keys(x, "TXID")
  checkTrue(is.character(head(res)))
  checkTrue(length(res) > 70000)

  res <- OrganismDbi:::.keys(x, "ENTREZID")
  checkTrue(is.character(head(res)))
  checkTrue(length(res) > 42000)
}

test_mkeys <- function(){
  tbl1 <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
  tbl2 <- "org.Hs.eg.db"
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl1")
  checkTrue("GENEID"==res)
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl2")
  checkTrue("ENTREZID"==res)

  tbl1 <- "GO.db"
  tbl2 <- "org.Hs.eg.db"
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl1")
  checkTrue("GOID"==res)
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl2")
  checkTrue("GO"==res)

  tbl1 <- "GO.db"
  tbl2 <- "org.Hs.eg.db"
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="both")
  res2 <- c("GOID","GO")
  names(res2) <- c("GO.db","org.Hs.eg.db")
  checkEquals(res, res2)
}


test_getSelects <- function(){
  allCols <- OrganismDbi:::.colsByNodes(x)

  ## start at one end case
  cols <- c("TERM", "ALIAS")
  keytype <- "GENEID"
  keys <- head(keys(x, keytype), n=2)
  
  subgr <- OrganismDbi:::.getRelevantSubgraph(x, cols=cols, keys,
                                              keytype=keytype)
  root <- OrganismDbi:::.lookupDbNameFromKeytype(x, keytype)
  fKeys <- OrganismDbi:::.getForeignKeys(x, subgr)
  selectCols <- unique(c(keytype, fKeys, cols))
  needCols <- OrganismDbi:::.getColsByNodes(subgr, selectCols, allCols)
  visitNodes <- OrganismDbi:::.bfs(subgr, root)
  
  res <- OrganismDbi:::.getSelects(x, keytype, keys, needCols, visitNodes)
 
  checkTrue(length(res)==3)
  checkTrue(class(res)=="list")
  checkTrue("GENEID" %in% colnames(res[[1]]))
  checkTrue("GO" %in% colnames(res[[2]]))
  checkTrue("TERM" %in% colnames(res[[3]]))

  ## The very simple case
  cols <- c("SYMBOL")
  keytype <- "OMIM"
  keys <- head(keys(x, keytype), n=2)
  
  subgr <- OrganismDbi:::.getRelevantSubgraph(x, cols=cols, keys,
                                              keytype=keytype)
  root <- OrganismDbi:::.lookupDbNameFromKeytype(x, keytype)
  fKeys <- OrganismDbi:::.getForeignKeys(x, subgr)
  selectCols <- unique(c(keytype, fKeys, cols))
  needCols <- OrganismDbi:::.getColsByNodes(subgr, selectCols, allCols)
  visitNodes <- OrganismDbi:::.bfs(subgr, root)

  res <- OrganismDbi:::.getSelects(x, keytype, keys, needCols, visitNodes)

  checkTrue(length(res)==1)
  checkTrue(class(res)=="list")
  checkTrue("OMIM" %in% colnames(res[[1]]))
  checkTrue("SYMBOL" %in% colnames(res[[1]]))


  
  ## Then there is this case (start in the middle case)
  cols <- c("GOID" ,  "SYMBOL", "TXNAME")
  keytype <- "ENTREZID"
  keys <- head(keys(x, "ENTREZID"))
  
  subgr <- OrganismDbi:::.getRelevantSubgraph(x, cols=cols, keys,
                                              keytype=keytype)
  root <- OrganismDbi:::.lookupDbNameFromKeytype(x, keytype)
  fKeys <- OrganismDbi:::.getForeignKeys(x, subgr)
  selectCols <- unique(c(keytype, fKeys, cols))
  needCols <- OrganismDbi:::.getColsByNodes(subgr, selectCols, allCols)
  visitNodes <- OrganismDbi:::.bfs(subgr, root)

  res <- OrganismDbi:::.getSelects(x, keytype, keys, needCols, visitNodes)

  checkTrue(length(res)==3)
  checkTrue(class(res)=="list")
  checkTrue("SYMBOL" %in% colnames(res[[1]]))
  checkTrue("TXNAME" %in% colnames(res[[2]]))
  checkTrue("GOID" %in% colnames(res[[3]]))
  
}

test_mergeSelectResults <- function(){
  allCols <- OrganismDbi:::.colsByNodes(x)

  cols <- c("GOID" ,  "SYMBOL", "TXNAME")
  keytype <- "ENTREZID"
  keys <- head(keys(x, "ENTREZID"))
  subgr <- OrganismDbi:::.getRelevantSubgraph(x, cols=cols, keys,
                                              keytype=keytype)
  root <- OrganismDbi:::.lookupDbNameFromKeytype(x, keytype)
  fKeys <- OrganismDbi:::.getForeignKeys(x, subgr)
  selectCols <- unique(c(keytype, fKeys, cols))
  needCols <- OrganismDbi:::.getColsByNodes(subgr, selectCols, allCols)
  visitNodes <- OrganismDbi:::.bfs(subgr, root)
  selected <- OrganismDbi:::.getSelects(x, keytype, keys, needCols, visitNodes)
  res <- OrganismDbi:::.mergeSelectResults(x, selected, visitNodes)
  
  checkTrue(dim(res)[2]==8)
  checkTrue(class(res)=="data.frame")
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("TXNAME" %in% colnames(res))  
}


## MANY more tests
test_select <- function(){
  cls <- c("GO","ALIAS")
  keys <- head(keys(x, "ENTREZID"))
  keytype <- "ENTREZID"
  res <- OrganismDbi:::.select(x, keys, cls, keytype)  
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==5)
  checkTrue(class(res)=="data.frame")
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("EVIDENCE" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 

  cls <- c("IPI", "ALIAS", "CDSSTART") 
  res <- OrganismDbi:::.select(x, keys, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==4) 
  checkTrue("IPI" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("CDSSTART" %in% colnames(res)) 

  cls <- c("GOID","ENTREZID")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==4)
  checkTrue("GOID" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
 
  cls <- c("ALIAS","TXCHROM","EXONNAME")
  res <- OrganismDbi:::.select(x, keys, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==4) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("EXONNAME" %in% colnames(res)) 
  
  cls <- c("ACCNUM","CDSSTART") 
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))

  cls <- c("ACCNUM", "ALIAS")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))
  checkTrue("ALIAS" %in% colnames(res))

  cls <- c("CDSSTART","CDSEND")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))
  checkTrue("CDSEND" %in% colnames(res))

  cls <- c("CDSSTART")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==2)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))

  cls <- c("ENTREZID")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==1)
  checkTrue("ENTREZID" %in% colnames(res))

  keys <- head(keys(x, "ENTREZID"))
  keytype <- "ENTREZID"
  cls <- c("GOID" ,  "SYMBOL", "TXNAME")
  res <- select(Homo.sapiens, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==6)
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("GOID" %in% colnames(res)) 
  checkTrue("SYMBOL" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res)) 

##   ## This tests for fields that are not in the final Homo.sapiens pkg
##   ## I am keeping it because it may be of use if I decide to add hom pkgs
##   cls <- c("ALIAS", "ORYZA_SATIVA")
##   res <- select(Homo.sapiens, keys, cls, keytype)
##   checkTrue(dim(res)[1] >0)
##   checkTrue(dim(res)[2]==3)
##   checkTrue("ENTREZID" %in% colnames(res)) 
##   checkTrue("ALIAS" %in% colnames(res)) 
##   checkTrue("ORYZA_SATIVA" %in% colnames(res))

  ## Getting an error here 
  keys <- head(keys(x, "TXNAME"))
  keytype <- "TXNAME"
  cls <- c("ENTREZID" , "TXNAME")
  res <- select(Homo.sapiens, keys, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==2)
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res))   
}


## Also need to test a species with other keys to join DBs

test_rattus <- function(){ 
  cls <- c("GO","ALIAS") 
  k <- head(keys(r, "ENTREZID")) 
  keytype <- "ENTREZID" 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==5) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ONTOLOGY" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 

  cls <- c("GO","ALIAS","TXCHROM","TXNAME") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==7) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("TXCHROM" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res)) 

  cls <- c("TXCHROM","TXNAME") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==3) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("TXCHROM" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res)) 

  ## now test different keytype
  k <- head(keys(r, keytype="ENSEMBL"))
  keytype <- "ENSEMBL"
  cls <- c("GO","ALIAS","TXCHROM","TXNAME") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==7) 
  checkTrue("ENSEMBL" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ONTOLOGY" %in% colnames(res)) 
  checkTrue("TXCHROM" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res)) 

  ## now test key that starts us from TxDb
  k <- head(keys(r, keytype="TXNAME"))
  keytype <- "TXNAME"
  cls <- c("GO","ALIAS","TXCHROM") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==6) 
  checkTrue("TXNAME" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("EVIDENCE" %in% colnames(res)) 
  checkTrue("TXCHROM" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 

  ## now test key that starts us from Go
  ## TODO: A cleanup bug??? <- Row of NAs in 1st line
  k <- head(keys(r, keytype="GOID"))
  keytype <- "GOID"
  cls <- c("GOID","ALIAS","TXCHROM") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==5) 
  checkTrue("GOID" %in% colnames(res)) 
  checkTrue("ONTOLOGY" %in% colnames(res)) 
  checkTrue("TXCHROM" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res))


  ## what happens when we use a key from the middle?
  k <- keys <- head(keys(r,keytype="ENTREZID"))
  cls <- c("GOID","SYMBOL","TXNAME")
  keytype <- "ENTREZID"
  res <- OrganismDbi:::.select(r, k, cls, keytype)
  checkTrue(dim(res)[1] >0)
  checkTrue(dim(res)[2]==6) 
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("GOID" %in% colnames(res)) 
  checkTrue("SYMBOL" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res))   
} 






test_dbconn_and_dbfile <- function(){
    resc <- dbconn(x)[[1]]
    m <- dbGetQuery(resc, "SELECT * FROM metadata")
    checkTrue(dim(m)[2] ==2)
    checkTrue(dim(m)[1] > 6)
              
    resf <- dbfile(x)[[1]]
    mf <- dbGetQuery(dbConnect(SQLite(), resf), "SELECT * FROM metadata")
    checkTrue(all(mf == m))
}


test_selectByRanges <- function(){
    ranges <-  GRanges(seqnames=Rle(c('chr11'), c(2)),IRanges(start=c(107899550, 108025550), end=c(108291889, 108050000)), strand='+', seqinfo=seqinfo(Homo.sapiens))
    
    res <- selectByRanges(x, ranges, 'SYMBOL', 'tx')
    checkTrue(length(res) > 2)
    checkTrue('ACAT1' %in% res$SYMBOL[[1]])
    checkTrue(class(res)=="GRanges")

    res2 <- selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), '5utr')
    checkTrue(length(res2) > 1)
    checkTrue('ACAT1' %in% res2$SYMBOL[[1]])
    checkTrue('04110' %in% res2$PATH[[2]])
    checkTrue(length(res2$PATH[[2]]) > 1)
    checkTrue(class(res2)=="GRanges")
}


test_selectRangesById <- function(){
    ## notice one of my keys is 'bad' (but things still work)
    res <- selectRangesById(x, c('bob','A1BG'), columns='SYMBOL',
                            keytype='SYMBOL')
    checkTrue(length(res[[1]]) > 1)
    checkTrue('A1BG' %in% res[[1]]$SYMBOL[[1]])
    checkTrue('A1BG' %in% names(res))
    checkTrue(class(res[[1]])=="GRanges")
    
    ##In the case of bad keys it should still work, but some will be filtered
    symbols <- keys(x, 'SYMBOL')[1:10]
    res2 <- selectRangesById(x,keys=symbols,keytype='SYMBOL')
    checkTrue(class(res2[[1]])=="GRanges")
    checkTrue('A1BG' %in% names(res2))
}











## Fast testing: BiocGenerics:::testPackage(pattern="^test_select.*\\.R$")


## TODO: add something to fix the the cosmetic bug where then the GOID is to the right of the columns that come with it (like EVIDENCE and/or ONTOLOGY.  This should really be handled in a general way (even though it ONLY happens with GOID)


