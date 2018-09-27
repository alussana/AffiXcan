# AffiXcan
AffiXcan is an R package that includes a set of functions to train and to apply
statistical models to estimate GReX (genetically regulated expression).

## Background
Understanding and predicting how genetic variation influences gene
expression in cells and tissues is of great interest in modern biological and
medical sciences.

The present methods to estimate the genetic contribution to gene expression do
not take into account functional information in identifying _expression
quantitative trait loci_ (eQTL), i.e. those genetic variants that contribute to
explaining the variation of gene expression. Relying on SNPs as predictors
allows to make significant models of gene expression only for those genes for
which SNPs with a fairly good effect size exist, but this condition is not
satisfied for the majority of genes, despite their expression having a non-zero
heritability (h<sup>2</sup>). To address this issue, new, different strategies
to analyze genetic variability of regulatory regions and their influence on
transcription are needed.

## General features
__AffiXcan__ (total binding AFFInity-eXpression sCANner) implements a functional
approach based on the
[TBA](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143627)
(Total Binding Affinity) score to make statistical models of gene expression,
being able to make significant predictions on genes for which SNPs with strong
effect size are absent. Furthermore, such a functional approach allows to make
mechanistic interpretations in terms of transcription factors binding events
that drive differential transcriptional expression. These features are of
considerable importance for eQTL discovery and to improve the capability to
estimate a [GReX](#grex) (genetically regulated expression) for a greater number
of genes, at the same time giving insights on the possible molecular mechanisms
that are involved in differential expression on genes.


