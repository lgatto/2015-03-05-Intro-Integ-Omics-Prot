# Integrating transcriptomics and proteomics data

Different approaches to data integration exist: (i) model- or
network-based approaches that identify common patterns in different
data sources as well as (ii) reference-based approaches, that map
different sources of data against a common reference.

The former are very versatile and rely on experiment-wide
clustering/modelling and crucially depend on reliably linking features
(explicitly via common identifiers or through functional
contextualisation). Transcript and protein measurement have previously
been combined and compared by linking the respective features by a
**common (gene) identifier** 
([1](http://www.ncbi.nlm.nih.gov/pubmed/21179022),
[2](http://www.ncbi.nlm.nih.gov/pubmed/22068331)). Such approaches are
often difficult to track and are susceptible to inconsistencies in the
relation between different data sources when, for example, multiple
transcripts are compared to ambiguous protein groups. The latter
approaches are a natural choice for data stemming from genomics,
transcriptomics, epigenomics, etc that directly rely on **mapping** their
data features along a **genome reference**.

# Mapping peptides to genomic coordinates



The **goal** is to map peptides from protein coordinates (1 to *L_p*)
to genomic coordinates.

![plot of chunk mapplot](figure/mapplot-1.png) 

Illustration with the
[`Pbase`](http://bioconductor.org/packages/devel/bioc/html/Pbase.html)
Bioconductor package (devel version).

## Data

We have an example data, named `p`, composed of 9
proteins, with UniProt accession numbers and Ensembl transcipt
identifiers and each protein has a set experimentally observed peptides
(see table below). This `p` object is generated from the protein
database (fasta file) and the MS identification results (`mzIdentML`
file) against this very same protein database.


|Acc      |ENST            |npep |
|:--------|:---------------|:----|
|A4UGR9   |ENST00000409195 |36   |
|A6H8Y1   |ENST00000358731 |23   |
|O43707   |ENST00000252699 |6    |
|O75369   |ENST00000295956 |13   |
|P00558   |ENST00000373316 |5    |
|P02545   |ENST00000368300 |12   |
|P04075   |ENST00000338110 |21   |
|P04075-2 |ENST00000395248 |20   |
|P60709   |ENST00000331789 |1    |

For example, P00558:

![plot of chunk pplot1](figure/pplot1-1.png) 

## Genomic coordinates of the transcripts/exons


```r
grl <- etrid2grl(acols(p)$ENST)
pcgrl <- proteinCoding(grl)
```

![plot of chunk gviz1](figure/gviz1-1.png) 

## Mapping peptides to the genome


```r
res <- pmapToGenome(p, pcgrl)
```

![plot of chunk gviz2](figure/gviz2-1.png) 

## Detailed annotation tracks

Maintaining access to the raw MS data (used as input with the fasta
file to generate the identification results).

![plot of chunk gviz3](figure/gviz3-1.png) 


# Multiple transcipts per protein

## Data

We use our example data, named `p`, composed of 9
proteins, with UniProt accession numbers and Ensembl transcipt
identifers and each protein has a set experimentally observed
peptides:


|Acc      |ENST            |npep |
|:--------|:---------------|:----|
|A4UGR9   |ENST00000409195 |36   |
|A6H8Y1   |ENST00000358731 |23   |
|O43707   |ENST00000252699 |6    |
|O75369   |ENST00000295956 |13   |
|P00558   |ENST00000373316 |5    |
|P02545   |ENST00000368300 |12   |
|P04075   |ENST00000338110 |21   |
|P04075-2 |ENST00000395248 |20   |
|P60709   |ENST00000331789 |1    |

If we hadn't the curated UniProt accession/Ensembl transcript
identifer maps, we would, for example, query an online repositary such
as the Ensembl Biomart instance.


```r
library("biomaRt")
ens <- useMart("ensembl", "hsapiens_gene_ensembl")
ens
```

```
## Object of class 'Mart':
##  Using the ensembl BioMart database
##  Using the hsapiens_gene_ensembl dataset
```

<!-- ## Updated on 16-06-2015: -->
<!-- ## replaced uniprot_swissprot_accession with uniprot_swissprot -->


```r
bm <- select(ens, keys = seqnames(p),
             keytype = "uniprot_swissprot",
             columns = c(
                 "uniprot_swissprot",
                 "ensembl_transcript_id"))

bm
```

```
##    uniprot_swissprot ensembl_transcript_id
## 1             A4UGR9       ENST00000409043
## 2             A4UGR9       ENST00000409728
## 3             A4UGR9       ENST00000409195
## 4             A4UGR9       ENST00000409273
## 5             A4UGR9       ENST00000409605
## 6             A6H8Y1       ENST00000617085
## 7             A6H8Y1       ENST00000358731
## 8             O43707       ENST00000252699
## 9             O43707       ENST00000390009
## 10            O75369       ENST00000490882
## 11            O75369       ENST00000295956
## 12            O75369       ENST00000358537
## 13            O75369       ENST00000429972
## 14            P00558       ENST00000373316
## 15            P02545       ENST00000368301
## 16            P02545       ENST00000368300
## 17            P02545       ENST00000368299
## 18            P02545       ENST00000448611
## 19            P02545       ENST00000473598
## 20            P02545       ENST00000347559
## 21            P04075       ENST00000338110
## 22            P04075       ENST00000395248
## 23            P04075       ENST00000566897
## 24            P04075       ENST00000569545
## 25            P04075       ENST00000563060
## 26            P04075       ENST00000412304
## 27            P04075       ENST00000564546
## 28            P04075       ENST00000564595
## 29            P60709       ENST00000331789
```



If we focus on P02545 for example, we see that we
retrieve 6 possible transcript identifers, including our
annotated ENST00000368300.


|uniprot_swissprot |ensembl_transcript_id |
|:-----------------|:---------------------|
|P02545            |ENST00000368301       |
|P02545            |ENST00000368300       |
|P02545            |ENST00000368299       |
|P02545            |ENST00000448611       |
|P02545            |ENST00000473598       |
|P02545            |ENST00000347559       |

## Genomic coordinates

Let's fetch the coordinates of all possible transcipts, making sure
that the names of the Ensembl identifiers are used to name the grl
ranges (using `use.names = TRUE`). We obtain 30 sets of ranges for 9
proteins.


```r
eid
```

```
##            P02545            P02545            P02545            P02545 
## "ENST00000368301" "ENST00000368300" "ENST00000368299" "ENST00000448611" 
##            P02545            P02545 
## "ENST00000473598" "ENST00000347559"
```

```r
grl <- etrid2grl(eid, ens, use.names = TRUE)
pcgrl <- proteinCoding(grl)
```


```r
grTr <- lapply(pcgrl, function(i)
    GeneRegionTrack(i, name = mcols(i)$transcript[1]))

plotTracks(grTr)
```

![plot of chunk givtr](figure/givtr-1.png) 

## Discriminating transcripts



We extract the transcript sequences, translate them into protein
sequences and align each to our original protein sequence.


```r
library("BSgenome.Hsapiens.NCBI.GRCh38")
lseq <- lapply(getSeq(BSgenome.Hsapiens.NCBI.GRCh38, pcgrl),
               function(s) translate(unlist(s)))

laln <- sapply(lseq, pairwiseAlignment, aa(p[k]))
sapply(laln, nmatch)/width(aa(p[k]))
```

```
##    P02545    P02545    P02545    P02545    P02545    P02545 
## 0.8614458 1.0000000 0.9246988 0.8298193 0.8358434 0.9548193
```



We see that transcript number 2, ENST00000368300, perfectly aligns
with our protein sequence. This is also the transcipt that corresponds
to the curated Ensembl transcript in `acols(p)$ENST`.




```r
res <- pmapToGenome(p[k], pcgrl[ki])
```

![plot of chunk pepcoords2](figure/pepcoords2-1.png) 

# Mapping MS peptides and RNA-Seq short reads

The last step of the mapping process is the combine the newly mapped
peptides and reads from RNA-Seq experiments. The figures below
illustrate this with data from Sheynkman et
al. ([2013](http://www.ncbi.nlm.nih.gov/pubmed/23629695),
[2014](http://www.ncbi.nlm.nih.gov/pubmed/25149441)) from the Jurkat
cell line (TIB-152). The
[mass spectrometry](https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/PASS_View?identifier=PASS00215)
and
[RNA-Seq](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45428)
(SRR791580) where processed with standard pipelines. 


![plot of chunk peptrack](figure/peptrack-1.png) ![plot of chunk peptrack](figure/peptrack-2.png) 

![plot of chunk alntrack](figure/alntrack-1.png) 

### References

Laurent Gatto and Sebastian Gibb
(2014). [`Pbase`](http://bioconductor.org/packages/devel/bioc/html/Pbase.html):
Manipulating and exploring protein and proteomics data. R package
version 0.6.9. https://github.com/ComputationalProteomicsUnit/Pbase

Pang et al. Tools to covisualize and coanalyze proteomic data with
genomes and transcriptomes: validation of genes and alternative mRNA
splicing. J Proteome Res. 2014 Jan 3;13(1):84-98. doi:
10.1021/pr400820p. Epub 2013 Nov 12. PubMed
[PMID: 24152167](http://www.ncbi.nlm.nih.gov/pubmed/24152167).

## Session information


```r
sessionInfo()
```

```
## R version 3.2.0 Patched (2015-04-22 r68234)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.2 LTS
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] grid      stats4    parallel  methods   stats     graphics  grDevices
##  [8] utils     datasets  base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.NCBI.GRCh38_1.3.1000
##  [2] BSgenome_1.37.1                       
##  [3] Biostrings_2.37.2                     
##  [4] biomaRt_2.25.1                        
##  [5] ggplot2_1.0.1                         
##  [6] XVector_0.9.1                         
##  [7] knitr_1.10.5                          
##  [8] rtracklayer_1.29.10                   
##  [9] Pbase_0.9.0                           
## [10] Gviz_1.13.2                           
## [11] GenomicRanges_1.21.15                 
## [12] GenomeInfoDb_1.5.7                    
## [13] IRanges_2.3.11                        
## [14] S4Vectors_0.7.5                       
## [15] Rcpp_0.11.6                           
## [16] BiocGenerics_0.15.2                   
## 
## loaded via a namespace (and not attached):
##  [1] Biobase_2.29.1             vsn_3.37.1                
##  [3] splines_3.2.0              foreach_1.4.2             
##  [5] Formula_1.2-1              highr_0.5                 
##  [7] affy_1.47.1                latticeExtra_0.6-26       
##  [9] Rsamtools_1.21.8           impute_1.43.0             
## [11] RSQLite_1.0.0              lattice_0.20-31           
## [13] biovizBase_1.17.1          limma_3.25.9              
## [15] chron_2.3-45               digest_0.6.8              
## [17] RColorBrewer_1.1-2         colorspace_1.2-6          
## [19] preprocessCore_1.31.0      plyr_1.8.2                
## [21] MALDIquant_1.12            XML_3.98-1.2              
## [23] zlibbioc_1.15.0            scales_0.2.4              
## [25] affyio_1.37.0              cleaver_1.7.0             
## [27] BiocParallel_1.3.25        SummarizedExperiment_0.1.5
## [29] GenomicFeatures_1.21.13    nnet_7.3-9                
## [31] proto_0.3-10               survival_2.38-1           
## [33] magrittr_1.5               evaluate_0.7              
## [35] doParallel_1.0.8           MASS_7.3-40               
## [37] foreign_0.8-63             mzR_2.3.1                 
## [39] BiocInstaller_1.19.6       Pviz_1.3.0                
## [41] tools_3.2.0                data.table_1.9.4          
## [43] formatR_1.2                matrixStats_0.14.0        
## [45] stringr_1.0.0              MSnbase_1.17.5            
## [47] munsell_0.4.2              cluster_2.0.1             
## [49] AnnotationDbi_1.31.16      lambda.r_1.1.7            
## [51] pcaMethods_1.59.0          mzID_1.7.0                
## [53] futile.logger_1.4.1        RCurl_1.95-4.6            
## [55] dichromat_2.0-0            iterators_1.0.7           
## [57] VariantAnnotation_1.15.13  labeling_0.3              
## [59] bitops_1.0-6               gtable_0.1.2              
## [61] codetools_0.2-11           DBI_0.3.1                 
## [63] reshape2_1.4.1             GenomicAlignments_1.5.9   
## [65] gridExtra_0.9.1            Hmisc_3.16-0              
## [67] ProtGenerics_1.1.0         futile.options_1.0.0      
## [69] stringi_0.4-1              rpart_4.1-9               
## [71] acepack_1.3-3.3
```

| [Home](./README.md) | [Caveats](./caveats.md) | [Mapping](./mapping.md) | [Transfer learning](./transfer-learning.md) |
