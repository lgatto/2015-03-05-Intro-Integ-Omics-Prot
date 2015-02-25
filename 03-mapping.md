# Mapping peptides to genomic coordinates

Illustration with [`Pbase`](http://bioconductor.org/packages/devel/bioc/html/Pbase.html).



We have an example data composed of 9 proteins, with
UniProt accession numbers and Ensembl transcipt identifers:


```r
kable(cbind(Acc = seqnames(p), ENST = acols(p)$ENST))
```



|         |Acc      |ENST            |
|:--------|:--------|:---------------|
|A4UGR9   |A4UGR9   |ENST00000409195 |
|A6H8Y1   |A6H8Y1   |ENST00000358731 |
|O43707   |O43707   |ENST00000252699 |
|O75369   |O75369   |ENST00000295956 |
|P00558   |P00558   |ENST00000373316 |
|P02545   |P02545   |ENST00000368300 |
|P04075   |P04075   |ENST00000338110 |
|P04075-2 |P04075-2 |ENST00000395248 |
|P60709   |P60709   |ENST00000331789 |

And each of these proteins have a set experimentally observed
peptides:


| A4UGR9| A6H8Y1| O43707| O75369| P00558| P02545| P04075| P04075-2| P60709|
|------:|------:|------:|------:|------:|------:|------:|--------:|------:|
|     36|     23|      6|     13|      5|     12|     21|       20|      1|

![plot of chunk pplot1](figure/pplot1-1.png) 

## Genomic coordinates of the transcripts/exons


```r
grl <- etrid2grl(acols(p)$ENST)
pcgrl <- proteinCoding(grl)
```

![plot of chunk gviz1](figure/gviz1-1.png) 
