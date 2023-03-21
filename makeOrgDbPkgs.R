## script to test my package code generator
require(OrganismDbi)
version = "1.3.0"

## for human
gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
           join2 = c(org.Hs.eg.db="ENTREZID",
                     TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID"))

makeOrganismPackage(pkgname = "Homo.sapiens",
                    graphData = gd,
                    organism = "Homo sapiens",
                    version = version,
                    maintainer =
              "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                    author = "Bioconductor Core Team",
                    destDir = ".",
                    license = "Artistic-2.0")



## for mouse
gd <- list(join1 = c(GO.db="GOID", org.Mm.eg.db="GO"),
           join2 = c(org.Mm.eg.db="ENTREZID",
                     TxDb.Mmusculus.UCSC.mm10.knownGene="GENEID"))

makeOrganismPackage(pkgname = "Mus.musculus",
                    graphData = gd,
                    organism = "Mus musculus",
                    version = version,
                    maintainer =
              "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                    author = "Bioconductor Core Team",
                    destDir = ".",
                    license = "Artistic-2.0")


## for rat
gd <- list(join1 = c(GO.db="GOID", org.Rn.eg.db="GO"),
           join2 = c(org.Rn.eg.db="ENTREZID",
                     TxDb.Rnorvegicus.UCSC.rn5.refGene="GENEID"))

makeOrganismPackage(pkgname = "Rattus.norvegicus",
                    graphData = gd,
                    organism = "Rattus norvegicus",
                    version = version,
                    maintainer =
              "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                    author = "Bioconductor Core Team",
                    destDir = ".",
                    license = "Artistic-2.0")
