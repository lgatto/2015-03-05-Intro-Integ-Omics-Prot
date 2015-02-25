# Mapping peptides to genomic coordinates

Illustration with [`Pbase`](http://bioconductor.org/packages/devel/bioc/html/Pbase.html).


Different approaches to integration exist: (i) model- or network-based
approaches that identify common patterns in different data sources as
well as (ii) reference-based approaches, that map different sources of
data against a common reference. The former are very versatile and
rely on experiment-wide clustering/modelling
\cite{Le_Cao:2011,Haider:2013, Meng:2014} and crucially depend on
reliably linking features (explicitly via common identifiers or
through functional contextualisation). Transcript and protein
measurement have previously been combined and compared by linking the
respective features by a common gene identifier
\cite{Lundberg:2010,Nagaraj:2011}. Such approaches are often difficult
to track and are susceptible to inconsistencies in the relation
between different data sources when, for example, multiple transcripts
are compared to ambiguous protein groups. The latter approaches are a
natural choice for data stemming from genomics, transcriptomics,
epigenomics, etc that directly rely on mapping their data features
along a genome reference.

#### Mapping, database reduction, database expansion

- [Pang et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/24152167)
- [Evans et al. (2012)](http://www.ncbi.nlm.nih.gov/pubmed/23142869)
- [Sheynkman et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25149441)
- [Qeli and Ahrens (2010)](http://www.ncbi.nlm.nih.gov/pubmed/20622826)

- [Wang and Zhang (2013)](http://www.ncbi.nlm.nih.gov/pubmed/24058055) [`customProDB`](http://bioconductor.org/packages/release/bioc/html/customProDB.html)

- [Boekel et al (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25658277)

