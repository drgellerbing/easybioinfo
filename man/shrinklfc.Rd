% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrinklfc.R
\name{shrinklfc}
\alias{shrinklfc}
\title{Shrink Log2 Fold Changes}
\usage{
shrinklfc(dds)
}
\arguments{
\item{dds}{a DESeqDataSet object, either after running DESeq from DESeq2 package or rundeseq function}
}
\value{
A list that contains the DESeqResults object with shrunken LFC, the full shrunken result dataframe, and the significant shrunken result dataframe.
Both dataframes can be saved into excel files
}
\description{
A simplified pipeline for the shrinkage of effect size (LFC estimates) which is useful for visualisation and ranking of genes.
Shrinks the DESeq2 results and filters the significant results based on your desired LFC and p-value that can be saved into excel files
Shrinkage of genes can be done using three different estimators: ashr/normal/apeglm(default).
Please visit http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.htmlhttp://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html for more information about the shrinkage estimators.
}
\examples{
expression = easybioinfo::deseqexpr
md = easybioinfo::deseqmd

exampledeseq = DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = md, design = ~condition, tidy = TRUE)
dds <- DESeq2::DESeq(exampledeseq)

library(easybioinfo)
dds <- easybioinfo::rundeseq(expression, md)
lfc <- shrinklfc(dds)
}
