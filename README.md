
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GeoLineages

<!-- badges: start -->
<!-- badges: end -->

The goal of GeoLineages is to provide a set of functions that are used
to analyse the outputs from ancestral state analyses.

## Installation

You can install the latest version of GeoLineages this repository using

``` r
devtools::github_install("magnusnosnes/GeoLineages")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(GeoLineages)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
set.seed(400)
library(phytools)
#> Loading required package: ape
#> Loading required package: maps
library(ape)
library(phangorn)
library(BactDating)

#Simulate data from Norway and rest of the world
#set.seed(400)
tree_test = simdatedtree(nsam=10, dateroot=2000)
tree_test = ladderize(tree_test)
Q=matrix(c(0.5,0.5,0.5,0.5), nrow=2,ncol=2, byrow=F)
loc = c("Norway", "Norway","Norway","RoW", "RoW", "Norway", "Norway", "RoW", "RoW", "RoW")
names(loc) = tree_test$tip.label
mapped_Q = sim.history(tree_test, Q,nsim=10)
#> Some columns (or rows) of Q don't sum to 0.0. Fixing.
#> Done simulation(s).
colnames(Q)=c("Norway","RoW")

fit1 = ace(x=loc, phy= tree_test, type="discrete", mod="ARD")
plot.phylo(tree_test,lwd=2,label.offset = 0.15, mar=c(0.2,0.2,0.2,0.2))
axisPhylo(root.time=2000, backward=F)
nodelabels(pie=fit1$lik.anc,cex=0.7,piecol=c("Red","Blue"))
tips = to.matrix(loc,seq=c("Norway", "RoW"))
tiplabels(pie=tips, cex=0.7,piecol=c("Red","Blue"))
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
count_import_export_ace(tree_test, ace_nodes=fit1$lik.anc,
                        ace_tips = to.matrix(loc, seq=c("Norway", "RoW")), start_time=2000)
#> $Import_Export
#> [1] 4 6
#> 
#> $Lineage_sizes
#> [1] 3 2 2 3
#> 
#> $Taxa_names
#> $Taxa_names$`Lineage no: 1`
#> [1] "t6" "t7" "t1"
#> 
#> $Taxa_names$`Lineage no: 2`
#> [1] "t3" "t5"
#> 
#> $Taxa_names$`Lineage no: 3`
#> [1] "t4" "t2"
#> 
#> $Taxa_names$`Lineage no: 4`
#> [1] "t10" "t8"  "t9" 
#> 
#> 
#> $`MRCA's`
#> [1] 2000.000 2001.705 2007.369 2001.620
```
