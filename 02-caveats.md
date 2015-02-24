# Integrating proteomics data: caveats

## Similarities and differences with high-throughput sequencing data

**Mapping** of
- peptides along protein sequences (although not explicitly considered
  a mapping task)
- short reads along genome coordinates

BUT

1. the protein database and the genome are _independent_, i.e. the
   proteins do not make explicitly reference to the genome they
   originate from.

2. **coverage**: % vs X


```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.000   6.813  38.890  39.870  67.620 100.000
```

This has an impact on **protein inference* and **missing values**.
