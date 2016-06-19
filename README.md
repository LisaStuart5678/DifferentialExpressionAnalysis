---
title: "Differential Expression Analysis"
author: "Lisa Stuart"
date: "June 19, 2016"
output: html_document
---

Load the Bottomly data and perform a differential expression analysis using limma:

```{r}
library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
```

Build model matrix treating strain as the outcome but adjusting for lane as a factor:

```{r}
mod = model.matrix(~ pdata_bot$strain + pdata_bot$lane.number)
```

Fit model using limma and shrink using eBayes command:

```{r}
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
names(ebayes_limma)
```

Get topTable output to extract p-values:

```{r}
all = topTable(ebayes_limma, number = dim(edata)[1], sort.by="none")
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1], sort.by="none")$P.Value
```


Find genes significant at the 5% FDR rate using the Benjamini Hochberg correction to limit the p-values:

```{r}
fp_bh = p.adjust(limma_pvals,method="BH") 
sum(fp_bh < 0.05)
genes2 = as.numeric(fp_bh < 0.05)
names(genes2) = rownames(all)
```

Perform the gene set analysis with goseq following the protocol. Calculate the probability weight function by passing genes that are differentially expressed, the genome you are going to use, and that you are going to use ensGene for reference.

How many of the top 10 overrepresented categories are the same for the adjusted and unadjusted analysis?

```{r}
pwf=nullp(genes,"mm9","ensGene") # calc probability weight function by passing genes that are diff exp, genome going to use, and that going to use ensGene for references
pwf[1:100,]
GO.wall=goseq(pwf,"mm9","ensGene")
head(GO.wall)

pwf2=nullp(genes2,"mm9","ensGene") # calc probability weight function by passing genes that are diff exp, genome going to use, and that going to use ensGene for references
pwf2[1:100,]
GO.wall2=goseq(pwf2,"mm9","ensGene")
GO.wall2[1:10,]

```

Answer: 2

