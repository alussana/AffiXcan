# AffiXcan
AffiXcan is an R package that includes a set of functions to train and to apply
statistical models to estimate GReX (genetically regulated expression).

For usage, please see the package's [vignette](https://mega.nz/\
#!uX5y3IDR!cFClG3MazNDWs1IfXDLAtap3L6sV4hirwfLMWCNLZJA).

### What is GReX?
GReX is the component of gene expression (here defined as the transcript level,
e.g RPKM) explained by an individual's genetics.

The abundance of a transcript in a cell is determined by many factors, including
genetics, environmental factors, and disease. It can have an impact on the
cell's physiology and alter the expression of other transcripts or proteins,
their activity and regulation. Since transcription is initiated by the binding
of transcription factors to DNA, a portion of gene expression can be directly
explained by variants in _cis_ regulatory regions.

### Why GReX?
The estimation of GReX can be useful to perform TWAS when the real total
expression profile is unknown or can not be measured, for example in those
tissues - like the brain - that are inaccessible to _in vivo_ safe biopsies, or
in ancient genomes. 

GReX can be also exploited to estimate the constitutive susceptibility of a
genome to a certain status, the existence of which is at least partially
influenced by gene expression.

### Estimate GReX
Some efforts have been made to develop computational methods to predict GReX
from genotype data using mathematical models. 

[Gamazon et al.](http://www.nature.com/articles/ng.3367) developed a method
consisting of multiple-SNP prediction of expression levels, where the estimated
GReX for a gene is given by an additive model in which SNPs are the independent
variables.

__AffiXcan__ takes into account the contribution of all polymorphisms of given
genomic regions that are associated to the expression of a gene. This is done
using affinity scores - [TBA](https://journals.plos.org/plosone/\
article?id=10.1371/journal.pone.0143627) (Total Binding Affinity) - between
those regions and a set of transcription factors. A principal component analysis
(PCA) is performed on these scores and for each expressed gene a linear model is
fitted.

### AffiXcan Performance
We observed that the GReX of the majority of genes for which AffiXcan manages to
generate a significant model is not predictable by the method cited above.
Arguably, this is due to the nature of [TBA](https://journals.plos.org/plosone/\
article?id=10.1371/journal.pone.0143627) score, that allows to take into account
the additive small effect of all variants in a genomic region. Furthermore, the
goodness of prediction achieved by AffiXcan on both shared and non-shared genes
was significantly greater. For brief insights on AffiXcan's results in
preliminary tests, see AffiXcan Performance section in the package's
[vignette](https://mega.nz/\
#!uX5y3IDR!cFClG3MazNDWs1IfXDLAtap3L6sV4hirwfLMWCNLZJA).
