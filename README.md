
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FixingTaxonTraceR

<!-- badges: start -->
<!-- badges: end -->

The goal of FixingTaxonTraceR is to check sets of taxon sets for fixing
taxon traceability (FTT) in R as introduced in Algorithm 1 in paper
“Perfect taxon sampling and fixing taxon traceability: Introducing a
class of phylogenetically decisive taxon sets” (add link to
publication/preprint server).

In the setting $c=4$, FTT indicates phylogenetic decisiveness as defined
by Sanderson and Steel.

## Installation

You can install the laterst version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("pottj/FixingTaxonTraceR")
```

## Overview

There are four basic functions implemented in this package:

1.  *FTT_createInput*: The function reads in a .txt file with taxon
    lists (c-tuples or more) as indicated in the *fn* parameter and
    transforms them to target format (only c-tuples, numbered, ordered).
    Each row in the text file defines one set. The parameter *sepSym*
    allows different ways to separate the taxa within each set.
2.  *FTT_initialChecks*: The function performs up to three initial
    checks: First, the bounds of the FTT algorithm are checked. If the
    number of c-tuples is above the upper or below the lower bound, then
    the FTT algorithm does not to have to be used as the result is fix.
    In case of for $c=4$, three properties of phylogenetic decisiveness
    are tested:
    1.  Are all triples covered?
    2.  Are triples with simple covered covered by different quadruples?
    3.  Are all tuples sufficiently covered?
3.  *FTT_algorithmGreen*: This function tests for FTT using the input
    c-tuples and checking if there are $c$ c-tuples so a cross c-tuple
    can be resolved.
4.  *FTT_algorithmRed*: This function tests for FTT using the cross
    c-tuples and checking if there are $c$ c-tuples so it can be
    resolved.

**Important Note 1**: In the manuscript, the algorithm refers to “white”
and “gray” quadruple, while the R-implementation uses “green” and “red”,
respectively. This is only due to figure coloring costs, and does not
change the algorithm itself.

**Important Note 2**: Both *FTT_algorithmGreen* and *FTT_algorithmRed*
lead to the same result, but depending on the input set, the computing
time can differ, e.g. if there are only a few cross c-tuples, it will be
more efficient to run *FTT_algorithmRed*.

To validate the results for $c=4$ regarding phylogenetic decisiveness,
we also implemented the algorithm proposed by [Parvini, Braught, and
Fernández-Baca](https://ieeexplore.ieee.org/document/9616390) in
the function *FTT_findNRC*, which is based on no-rainbow coloring.

## Example

All examples of the publication are stored within the package. Please
check out the vignette for more details on how to run the package.
