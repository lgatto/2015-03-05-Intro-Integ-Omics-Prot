# Integrating proteomics data: caveats


**Mapping** of *peptides along protein sequences (although not
  explicitly considered a mapping task)* and *short reads along genome
  coordinates*.

But...

## Protein and gene identifers

The protein database and the genome are _independent_, i.e. the
proteins do not make explicitly reference to the genome they originate
from.


```r
## The UniProt human proteome (release 2015_02)
library("Pbase")
```

```
## Loading required package: methods
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist, unsplit
## 
## Loading required package: Rcpp
## Loading required package: Gviz
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: IRanges
## Loading required package: GenomicRanges
## Loading required package: GenomeInfoDb
## Loading required package: grid
## 
## This is Pbase version 0.6.5
```

```r
up <- Proteins("data/HUMAN_2015_02.fasta.gz")

## Using the accession number to Ensembl Biomart query
## for transcript identifiers
library("biomaRt")
ens <- useMart("ensembl", "hsapiens_gene_ensembl")
ens
```

```
## Object of class 'Mart':
##  Using the ensembl BioMart database
##  Using the hsapiens_gene_ensembl dataset
```

```r
upbm <- select(ens, keys = seqnames(up),
               keytype = "uniprot_swissprot_accession",
               columns = c(
                   "uniprot_swissprot_accession",
                   "ensembl_transcript_id"))
```

```r
## How many UniProt accession with Ensembl transcripts were found?
table(seqnames(up) %in% unique(upbm$uniprot_swissprot_accession))
```

```
## 
## FALSE  TRUE 
## 70781 19015
```


```r
## How many transcripts per accession do we find?
table(table(upbm$uniprot_swissprot_accession))
```

```
## 
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
## 8237 5056 2602 1404  713  383  219  112   78   56   35   17   16   15   12 
##   16   17   18   19   20   21   22   23   24   27   29   30   31   32   33 
##    8    7    4    6    5    5    2    2    3    1    1    1    2    3    1 
##   35   42   44   50   59   90 
##    3    2    1    1    1    1
```

## Coverage

**coverage**: % vs X


```r
cvg <- data.table::fread("./data/Ensembl_76.csv", skip = 17,
                         stringsAsFactors = FALSE)
summary(cvg$coverage)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.000   6.813  38.890  39.870  67.620 100.000
```

This has an impact on **protein inference* and **missing values**.
