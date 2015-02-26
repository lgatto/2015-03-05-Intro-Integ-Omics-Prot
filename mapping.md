# Mapping peptides to genomic coordinates



The **goal** is to map peptides from protein coordinates (1 to *L_p*)
to genomic coordinates.

<img src="figure/mapplot-1.png" title="plot of chunk mapplot" alt="plot of chunk mapplot" style="display: block; margin: auto;" />

Illustration with the
[`Pbase`](http://bioconductor.org/packages/devel/bioc/html/Pbase.html)
Bioconductor package (devel version).

## Data

We have an example data, named `p`, composed of 9
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

For example, P00558:

<img src="figure/pplot1-1.png" title="plot of chunk pplot1" alt="plot of chunk pplot1" style="display: block; margin: auto;" />

## Genomic coordinates of the transcripts/exons


```r
grl <- etrid2grl(acols(p)$ENST)
pcgrl <- proteinCoding(grl)
```

<img src="figure/gviz1-1.png" title="plot of chunk gviz1" alt="plot of chunk gviz1" style="display: block; margin: auto;" />

## Mapping peptides to the genome


```r
res <- pmapToGenome(p, pcgrl)
```

```
## Warning: Mapping failed. Returning an empty range.
```

```
## Warning: Mapping failed. Returning an empty range.
```

```
## Warning: Mapping failed. Returning an empty range.
```

```
## Warning: Mapping failed. Returning an empty range.
```


```
## Warning in plotTracks(list(ideoTrack, axisTrack, grTrack, pepTrack),
## groupAnnotation = "id", : The track chromosomes in 'trackList' differ.
## Setting all tracks to chromosome 'chrX'
```

<img src="figure/gviz2-1.png" title="plot of chunk gviz2" alt="plot of chunk gviz2" style="display: block; margin: auto;" />

## Detailed annotation tracks

Maintaining access to the raw MS data


```
## Error in .fillWithDefaults(data.frame(start = as.integer(start), end = as.integer(end)), : Number of elements in argument 'id' is invalid
```

```
## Error in plotTracks(list(ideoTrack, axisTrack, deTrack, grTrack), add53 = TRUE, : object 'deTrack' not found
```

See also
[Pang et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/24152167),
*Tools to covisualize and coanalyze proteomic data with genomes and
transcriptomes: validation of genes and alternative mRNA splicing.*