
# LineageHomology

- [Installation](#installation)
- [Introduction](#introduction)
- [Tutorial and gallery of plotting methods](#tutorial)

LineageHomology is an R-package that offers a set of functions for
analyzing the outputs of ancestral state reconstructions. The package
defines transmission lineages (TLs) and singletons based on state
probabilities at each node in the reconstruction. TLs are defined as
connected groups of tips where state transitions between ancestral and
descendant nodes have a probability lower than 50 percent, while
singletons are tips not connected to any other tips in this way. The
method used in LineageHomology is similar to that introduced by [du
Plessis et
al. (2021)](https://www.science.org/doi/10.1126/science.abf2946) (DOI:
10.1126/science.abf2946). The package also provides functions to
calculate the sizes of TLs, the number of singletons, estimates of
importation and local transmission, and other useful summaries, as well
as various plotting functions to visualize the results.

## Installation

You can install the latest version of LineageHomology by using devtools:

``` r
library(devtools)
devtools::install_github("magnusnosnes/LineageHomology")
```

## Introduction

This introduction provides a simple demonstration of how to apply
LineageHomology to analyze outputs of an ancestral state reconstruction
of simulated geographical data for two locations. We use LineageHomology
to describe transmission lineages and plot the sizes and estimated time
of the most recent common ancestor for each TL.

First, we simulate data and estimate the ancestral geographical states:

``` r
library(LineageHomology)
#Loading required packages.
library(ape)
#Loading other packages for simulating data. 
library(BactDating)

#Simulate data from Norway and rest of the world
tree_test = simdatedtree(nsam=10, dateroot=2000)
tree_test = ladderize(tree_test)
```

``` r
tiplabels=c("Sequence 7","Sequence 1","Sequence 5","Sequence 6","Sequence 4",
            "Sequence 3","Sequence 2","Sequence 9","Sequence 8","Sequence 10")
tree_test$tip.label=tiplabels
loc = c("Norway", "Norway","RoW","RoW", "RoW",
        "Norway", "Norway", "Norway", "RoW", "RoW")
names(loc) = tree_test$tip.label

#Reconstruct ancestral states using ace. 
fit1 = ace(x=loc, phy= tree_test, type="discrete", mod="ARD")
plot.phylo(tree_test,edge.width = 3,label.offset = 0.15, mar=c(0.2,0.2,0.2,0.2))
axisPhylo(root.time=2000, backward=F,lwd=2)
nodelabels(pie=fit1$lik.anc,cex=0.7,piecol=c("Red","Blue"))
tips = to.matrix(loc,seq=c("Norway", "RoW"))
tiplabels(pie=tips, cex=0.7,piecol=c("Red","Blue"))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" alt="Test" width="100%" />

**Figure 1:** The tree shows the reconstructed states using the ace
function in the ape package. Each node is coloured according as a
piechart with the probability of the location (red vs blue). For this
example red represents Norway, and blue represents the rest of the world
(RoW).

![LineageHomology_procedure](./Examples_and_plotting_methods/How_it_works.png)

**Figure 2:** This figure shows how LineageHomology estimates
transmission lineages, singletons, importation, exportation and local
transmission events. The green-shaded background represents transmission
lineages and encompasses the tips that are included in the group. The
transmission lineages are defined based on a 50 percent probability of
the same state on every node within the TL. The pink-shaded background
denotes singletons, which are unconnected to any other tips as per the
rule that defines TLs. The areas showing TLs and singletons extend back
to the importation date, which is the midpoint of the edge ancestral to
the most recent common ancestor (MRCA) for the TLs and the midpoint on
the ancestral edge leading to the first geographical transitions for
singletons. Exports are set to the estimated time at the midpoint of the
branch that leads to a different geographical location. Branching events
within TLs are used as estimates of local transmission. All types of
events are projected down as black points to individual axes to show the
events with respect to the focal location (red or blue / Norway or RoW).

Next we run LineageHomology on the output from ace according to the
approach outlined above.

``` r
Return = LineageHomology(tree_test, ace_nodes=fit1$lik.anc,
                        ace_tips = to.matrix(loc, seq=c("Norway", "RoW")), start_time=2000)
Return
#> $Import_LocalTrans
#> [1] 6 4
#> 
#> $Lineage_sizes
#> [1] 2 1 3 2 1 1
#> 
#> $Taxa_names
#> $Taxa_names$`Lineage no: 1`
#> [1] "Sequence 7" "Sequence 1"
#> 
#> $Taxa_names$`Lineage no: 2`
#> [1] "Sequence 5"
#> 
#> $Taxa_names$`Lineage no: 3`
#> [1] "Sequence 6"  "Sequence 4"  "Sequence 10"
#> 
#> $Taxa_names$`Lineage no: 4`
#> [1] "Sequence 3" "Sequence 2"
#> 
#> $Taxa_names$`Lineage no: 5`
#> [1] "Sequence 9"
#> 
#> $Taxa_names$`Lineage no: 6`
#> [1] "Sequence 8"
#> 
#> 
#> $`MRCA's`
#> [1] 2000.000 2004.735 2001.387 2005.682 2004.924 2007.846
#> 
#> $lineage_state
#> Norway    RoW    RoW Norway Norway    RoW 
#>      1      2      2      1      1      2 
#> 
#> $Halfedge_over_tmrca
#> [1] 2000.000 2003.343 2000.693 2004.679 2003.556 2006.385
```

In this example LineageHomology returned four 3 transmission lineages.
The taxa names of the tips included in each lineage is printed above
under “\$Taxa_names”. Lineages that only contain one tip are singletons.

The size distributions can be visualised by, e.g. using a treemap plot:

``` r
LineageHomology::treemap_lineagehomology(Return)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

**Figure 3:** The figure shows squares with areas representing the
transmission lineages’ sizes. The text on the squares gives the
estimated time of the MRCA and the number of tips in the transmission
lineages.

## Tutorial

See the full tutorial below for a full introduction to the package and
more methods for visualizing the results: [Tutorial and plotting
methods](https://github.com/magnusnosnes/LineageHomology/blob/master/Examples_and_plotting_methods/Simple_example/Basic_plotting.md)
