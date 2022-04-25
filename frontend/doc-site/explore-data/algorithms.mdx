# Algorithms used by Cellxgene

## Differential expression

> We're [actively working on how to improve differential expression](https://github.com/chanzuckerberg/cellxgene/issues/2211). **The** [**current implementation**](https://github.com/chanzuckerberg/cellxgene/blob/6f6634a4d9766a93596674fe42efbcae6ffabea6/backend/czi_hosted/compute/diffexp_generic.py#L42) **assumes normally distributed values on a linear scale.**

Cellxgene's differential expression uses a [Welch's _t_-test](https://en.wikipedia.org/wiki/Welch%27s_t-test), which assumes that the two populations are normally distributed, but may have unequal variance. We use a two-sided t-test against the null hypothesis that the two populations have **equal** means. P-values are adjusted with the [Bonferroni corrrection](https://en.wikipedia.org/wiki/Bonferroni_correction).

To help avoid spurious results, we only return genes with log fold change greater than 0.01:

`|log2( mean(set1) / mean(set2) )| > 0.01`

The log fold change threshold can be configured on cellxgene Desktop with the [`--diffexp-lfc-cutoff`](../desktop/launch.md) command. We then sort genes by their associated `|t value|` and return the top 15 genes.

