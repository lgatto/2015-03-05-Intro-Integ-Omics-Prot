# Mapping peptides to genomic coordinates



## Goal


```r
Pbase:::mapplot()
```

![plot of chunk mapplot](figure/mapplot-1.png) 

Illustration with [`Pbase`](http://bioconductor.org/packages/devel/bioc/html/Pbase.html).

## Data

We have an example data composed of 9 proteins, with
UniProt accession numbers and Ensembl transcipt identifers and each
protein has a set experimentally observed peptides:


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

Maintaining access to the raw MS data

![plot of chunk gviz3](figure/gviz3-1.png) 

See also
[Pang et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/24152167),
Tools to covisualize and coanalyze proteomic data with genomes and
transcriptomes: validation of genes and alternative mRNA splicing.
