
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Using LineageHomology for a given phylogeny and with tips observed in geographical locations.

#### Table of contents

-   [Estimating local transmission linages for phylogeographic data with
    two
    categories.](#estimating-local-transmission-linages-for-phylogeographic-data-with-two-categories)
-   [Simulate tip data with two
    locations](#simulate-phylogeographic-data-with-two-locations--norway-and-rest-of-the-world--row-)
-   [Running LineageHomology](#running-lineagehomology)
-   [Plotting transmission lineages](#plot-lineage-densities-over-time)
    -   [Color the lineages after
        location](#we-can-color-the-groups-by-the-state-they-are-in--and-restrict-the-plotted-groups-to-sizes-larger-than-1-4--and-10)
-   [Plot cumulative lineage size over
    time.](#plot-cumulative-lineage-size-over-time)
    -   [Color the lineages after
        location](#again-we-can-add-color-to-the-groups-by-specifying-it)
-   [Sample transmission lineages
    probabilistically](#probabilistic-counting-of-transmission-lineages)
-   [Produce confidence interval on the largest transmission
    lineage](#estimate-size-of-largest-transmission-lineage)

``` r
#Load needed packages
library(LineageHomology)
library(ggplot2)
library(scales)
library(lubridate)

#Loading other packages for simulating data. 
library(ape)
library(phytools)
library(phangorn)
library(BactDating)
```

##### Simulate tip data with two different state.

For the purpose of this example we let the geographical states represent
Norway and the rest of the World (RoW). The package was originally
developed geographical data with two locations.

``` r
set.seed(400)
tree_test = simdatedtree(nsam=300, dateroot=2000) #300 taxa and date of the root at year 2000-
tree_test = ladderize(tree_test) #Reorder the tree to make it look nice
Q=matrix(c(0.5,0.5,0.5,0.5), nrow=2,ncol=2, byrow=F) #Set up a transition matrix for trait simulation
colnames(Q)=c("Norway","RoW") #From state 1 to state 2. 
rownames(Q)=c("Norway","RoW") #From state 1 to state 2. 

trait = phytools::sim.Mk(tree=tree_test,Q=Q,nsim=1) #Simulate the traits on the phylogeny

#Reconstruct ancestral states using ace. 
fit1 = ace(x=trait, phy= tree_test, type="discrete", mod="ARD") #Estimate the ancestral history

plot.phylo(tree_test,lwd=2,label.offset = 0.15, mar=c(0.2,0.2,0.2,0.2),cex=0.3) #Ploy phylogeny

axisPhylo(root.time=2000, backward=F) #Add time axis

nodelabels(pie=fit1$lik.anc,cex=0.2,piecol=c("Red","Blue"))
tips = to.matrix(trait,seq=c("Norway", "RoW"))
tiplabels(pie=tips, cex=0.2,piecol=c("Red","Blue"))
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

The piecharts on the nodes of the phylogeny represents the probability
of the node being in the different geographical locations, where red
represents Norway and blue represents RoW. LineageHomology uses these
probabilities to count connected groups of tips, where all the nodes on
the paths between them have a probability that &gt; 50% for the same
geographical state. LineageHomology also returns the number of tips that
are not connected to any other tips in this way. Following du Plessis et
al. (2021) (DOI: 10.1126/science.abf2946) we refer to the groups as
transmission lineages (TLs) and isolated tips as singletons.
LineageHomology also returns other useful summaries

#### Running LineageHomology

Next we run LineageHomology. The inputs should take this format.

``` r
tree_test
#> 
#> Phylogenetic tree with 300 tips and 299 internal nodes.
#> 
#> Tip labels:
#>   t214, t151, t158, t56, t40, t13, ...
#> 
#> Rooted; includes branch lengths.
head(fit1$lik.anc)
#>          Norway       RoW
#> [1,] 0.19739205 0.8026080
#> [2,] 0.09659401 0.9034060
#> [3,] 0.09487369 0.9051263
#> [4,] 0.69688117 0.3031188
#> [5,] 0.43888708 0.5611129
#> [6,] 0.42001122 0.5799888
head( to.matrix(trait, seq=c("Norway", "RoW")))
#>      Norway RoW
#> t214      1   0
#> t151      0   1
#> t158      0   1
#> t56       0   1
#> t40       1   0
#> t13       0   1
```

Note that the node probabilities of the tips must have rownames that
match the tips of of the phylogeny.

``` r
Result = LineageHomology(tree_test, ace_nodes=fit1$lik.anc,
                        ace_tips = to.matrix(trait, seq=c("Norway", "RoW")), start_time=2000)
```

LineageHomology outputs a list of lists.

-   Import\_LocalTrans contains an estimate of the number of taxa that
    are have changed state wrt. the state of the closest common
    ancestral node shared with any other taxa.
-   Lineage\_sizes contains the number of taxa that belong to the mapped
    state. This means that for the entire group, all the interal nodes
    will be mapped to the same state with 50% or higher probability.
    \*Taxa names contains sub lists for each lineage, in the same order
    as the other lists.
    -   Each sub list contains the taxa names of all the taxa in the
        group.
-   MRCA’s contains the time of the most recent common ancestor of the
    groups.
-   Lineage\_state contains the state the groups are mapped to, the
    states are treated as numbers in the anlyses, so the relative
    ordering has to be looked up.

``` r
names(Result)
#> [1] "Import_LocalTrans"   "Lineage_sizes"       "Taxa_names"         
#> [4] "MRCA's"              "lineage_state"       "Halfedge_over_tmrca"
Result$Import_LocalTrans
#> [1] 120 180
Result$lineage_state
#> Norway    RoW Norway Norway Norway Norway    RoW    RoW Norway    RoW Norway 
#>      1      2      1      1      1      1      2      2      1      2      1 
#>    RoW Norway    RoW    RoW    RoW    RoW Norway Norway Norway Norway Norway 
#>      2      1      2      2      2      2      1      1      1      1      1 
#>    RoW Norway    RoW Norway Norway Norway    RoW Norway Norway Norway Norway 
#>      2      1      2      1      1      1      2      1      1      1      1 
#> Norway    RoW Norway    RoW Norway    RoW    RoW Norway Norway    RoW    RoW 
#>      1      2      1      2      1      2      2      1      1      2      2 
#> Norway    RoW Norway Norway Norway Norway    RoW Norway Norway Norway Norway 
#>      1      2      1      1      1      1      2      1      1      1      1 
#>    RoW    RoW Norway Norway Norway Norway Norway Norway Norway    RoW Norway 
#>      2      2      1      1      1      1      1      1      1      2      1 
#> Norway Norway    RoW Norway    RoW Norway Norway Norway    RoW Norway    RoW 
#>      1      1      2      1      2      1      1      1      2      1      2 
#> Norway    RoW Norway Norway    RoW Norway Norway    RoW    RoW    RoW Norway 
#>      1      2      1      1      2      1      1      2      2      2      1 
#>    RoW    RoW Norway    RoW Norway    RoW Norway    RoW    RoW Norway    RoW 
#>      2      2      1      2      1      2      1      2      2      1      2 
#> Norway Norway    RoW Norway    RoW    RoW    RoW Norway Norway    RoW    RoW 
#>      1      1      2      1      2      2      2      1      1      2      2 
#>    RoW    RoW Norway    RoW Norway    RoW    RoW Norway Norway    RoW 
#>      2      2      1      2      1      2      2      1      1      2
Result$Lineage_sizes
#>   [1]  1  9  1  2  1  8 18  1  3  1  1  1  2  1  5  7  1  3  1  8  2  1  2  1  1
#>  [26]  1  1  2  2  1  2  1  1  1  1  2  1  1  1  1 23  2  1 43  1  3  1  2  2  1
#>  [51]  2  1  2  1  8  1  1  3  1  1  1  1  1  4  1  1  1  2  1  1  1  1  2  1  1
#>  [76]  2  1  4  2  3  2  4  1  1  1  3  2  1  4  2  2  2  1  1  1  1  2  3  2  1
#> [101]  1  2  2  2  2  4  2  1  1  1  1  1  1  1  2  1  1  1  2  1
head(Result$Taxa_names)
#> $`Lineage no: 1`
#> [1] "t214"
#> 
#> $`Lineage no: 2`
#> [1] "t151" "t158" "t56"  "t13"  "t78"  "t67"  "t121" "t179" "t283"
#> 
#> $`Lineage no: 3`
#> [1] "t40"
#> 
#> $`Lineage no: 4`
#> [1] "t173" "t237"
#> 
#> $`Lineage no: 5`
#> [1] "t90"
#> 
#> $`Lineage no: 6`
#> [1] "t298" "t77"  "t142" "t141" "t18"  "t51"  "t12"  "t145"
Result$`MRCA's`
#>   [1] 2007.930 2003.575 2008.290 2008.563 2010.003 2002.453 2000.000 2012.839
#>   [9] 2008.899 2011.696 2015.350 2018.117 2017.801 2020.377 2016.547 2010.997
#>  [17] 2017.068 2014.607 2012.965 2013.090 2012.813 2012.869 2013.422 2014.381
#>  [25] 2014.017 2015.072 2012.170 2008.436 2009.600 2008.405 2010.026 2012.277
#>  [33] 2011.058 2008.982 2014.212 2007.166 2008.433 2010.349 2009.121 2010.671
#>  [41] 2000.831 2002.355 2002.036 2005.762 2015.719 2009.045 2013.681 2008.684
#>  [49] 2014.641 2017.424 2018.459 2023.799 2011.636 2013.835 2009.699 2013.235
#>  [57] 2014.098 2011.071 2013.331 2013.415 2013.012 2014.263 2021.411 2017.157
#>  [65] 2018.212 2009.985 2014.863 2016.249 2017.649 2011.083 2014.290 2016.099
#>  [73] 2011.329 2012.356 2013.407 2009.467 2010.733 2012.375 2002.329 2004.124
#>  [81] 2004.303 2006.608 2008.731 2009.043 2005.787 2012.925 2014.917 2014.136
#>  [89] 2015.214 2018.742 2017.105 2015.514 2013.682 2016.469 2015.963 2015.299
#>  [97] 2014.016 2010.969 2013.198 2017.861 2016.219 2009.868 2008.353 2011.624
#> [105] 2002.336 2007.166 2010.421 2008.355 2006.529 2007.034 2010.135 2012.856
#> [113] 2009.975 2009.423 2013.700 2013.159 2008.576 2011.045 2009.292 2012.001
Result$lineage_state
#> Norway    RoW Norway Norway Norway Norway    RoW    RoW Norway    RoW Norway 
#>      1      2      1      1      1      1      2      2      1      2      1 
#>    RoW Norway    RoW    RoW    RoW    RoW Norway Norway Norway Norway Norway 
#>      2      1      2      2      2      2      1      1      1      1      1 
#>    RoW Norway    RoW Norway Norway Norway    RoW Norway Norway Norway Norway 
#>      2      1      2      1      1      1      2      1      1      1      1 
#> Norway    RoW Norway    RoW Norway    RoW    RoW Norway Norway    RoW    RoW 
#>      1      2      1      2      1      2      2      1      1      2      2 
#> Norway    RoW Norway Norway Norway Norway    RoW Norway Norway Norway Norway 
#>      1      2      1      1      1      1      2      1      1      1      1 
#>    RoW    RoW Norway Norway Norway Norway Norway Norway Norway    RoW Norway 
#>      2      2      1      1      1      1      1      1      1      2      1 
#> Norway Norway    RoW Norway    RoW Norway Norway Norway    RoW Norway    RoW 
#>      1      1      2      1      2      1      1      1      2      1      2 
#> Norway    RoW Norway Norway    RoW Norway Norway    RoW    RoW    RoW Norway 
#>      1      2      1      1      2      1      1      2      2      2      1 
#>    RoW    RoW Norway    RoW Norway    RoW Norway    RoW    RoW Norway    RoW 
#>      2      2      1      2      1      2      1      2      2      1      2 
#> Norway Norway    RoW Norway    RoW    RoW    RoW Norway Norway    RoW    RoW 
#>      1      1      2      1      2      2      2      1      1      2      2 
#>    RoW    RoW Norway    RoW Norway    RoW    RoW Norway Norway    RoW 
#>      2      2      1      2      1      2      2      1      1      2
```

\#Make a treemap plot to get a quick overview of the lineages.

Treemap\_lineagehomology uses the treemap package to visualize
transmission lineages as rectangles with areas that reflect their
relative sizes. Additionally, the text on each square gives the name of
group number of the transmission lineage, the estimated time of the mrca
(TMRCA) of the transmission lineage and the number of tips in the group
(S)

``` r
LineageHomology::treemap_lineagehomology(Result)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

##### Plot lineage densities over time

Before plotting, we process the results using the function
lineage\_info. lineage\_info takes the results from LineageHomology and
names and dates observed dates of the tips in the tree. The formats are
shown below.

``` r
#Set up matrix with taxa info
name_date = data.frame(name = names(trait), dates= BactDating::leafDates(tree_test))
head(name_date)
#>   name    dates
#> 1 t214 2007.930
#> 2 t151 2005.666
#> 3 t158 2008.698
#> 4  t56 2008.041
#> 5  t40 2008.290
#> 6  t13 2006.730
Result_lineage_info = LineageHomology::lineage_info(Result,name_date)
```

Another plotting method is density plots of the transmission lineages.
Ridgeplot\_lineagedensities takes the input from the function
lineage\_info and produces density plots using for each TL over time.
The density plots are made using the ggridges R-package. The TLs are
sorted from largest to smallest, and can be colored according to the
state that defines the TL.

``` r
LineageHomology::ridgeplot_lineagedensities(Result_lineage_info=Result_lineage_info,groups_larger_than = 1,datelims=c("2000-01-01","2025-01-01","3 year"),color_by_state = F)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

##### We can color the groups by the state they are in, and restrict the plotted groups to sizes larger than 1,4, and 10

``` r
LineageHomology::ridgeplot_lineagedensities(Result_lineage_info=Result_lineage_info,groups_larger_than = 1,datelims=c("2000-01-01","2025-01-01","3 year"),color_by_state = T)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

``` r
LineageHomology::ridgeplot_lineagedensities(Result_lineage_info=Result_lineage_info,groups_larger_than = 4,datelims=c("2000-01-01","2025-01-01","3 year"),color_by_state = T)
```

<img src="man/figures/README-unnamed-chunk-10-2.png" width="100%" />

``` r
LineageHomology::ridgeplot_lineagedensities(Result_lineage_info=Result_lineage_info,groups_larger_than = 10,datelims=c("2000-01-01","2025-01-01","3 year"),color_by_state = T)
```

<img src="man/figures/README-unnamed-chunk-10-3.png" width="100%" />

##### Plot cumulative lineage size over time.

``` r
LineageHomology::lineage_growth_cumulative(Result_lineage_info = Result_lineage_info, datelims=c("2000-01-01","2025-01-01","3 year"))
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

##### Again we can add color to the groups by specifying it.

``` r
LineageHomology::lineage_growth_cumulative(Result_lineage_info = Result_lineage_info, datelims=c("2000-01-01","2025-01-01","3 year"),color_by_state = T)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

### Probabilistic counting of transmission lineages

For all the figures shown above, the transmission lineage division was
done by division of lineages at points where a transition is estimated
with more than 50 percent probability. If one wishes to rather count
transmission lineages by sampling from the observed probabilities, one
can use the function LineageHomology\_w\_uncertainty.

``` r
Result1 = LineageHomology_w_uncertainty_v2(tree_test, ace_nodes=fit1$lik.anc,
                        ace_tips = to.matrix(trait, seq=c("Norway", "RoW")), start_time=2000)
LineageHomology::treemap_lineagehomology(Result1)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

``` r
Result2 =  LineageHomology_w_uncertainty_v2(tree_test, ace_nodes=fit1$lik.anc,
                        ace_tips = to.matrix(trait, seq=c("Norway", "RoW")), start_time=2000)
LineageHomology::treemap_lineagehomology(Result2)
```

<img src="man/figures/README-unnamed-chunk-13-2.png" width="100%" />

As this figures show, the division of lineages will vary according to
the sampled states of some of the nodes. This can be useful for
producing confidence limits on e.g. the size of the largest lineage.

#### Estimate size of largest transmission lineage

We can use pbreplicate to counts transmission lineages probabilistically
many times to obtain confidence intervals (or, in this case, a
distribution) for the size of the largest transmission lineage.

``` r
library(pbapply)
multi_counts = pbreplicate(
  1000,
  LineageHomology::LineageHomology_w_uncertainty_v2(
    tree=tree_test,
    ace_nodes=fit1$lik.anc,
    ace_tips = to.matrix(trait, seq=c("Norway", "RoW")),
    start_time = 2000
  )
)

largest_lineage = c()
for(i in 1:1000) {
  largest_lineage = c(largest_lineage, max(multi_counts[,i]$Lineage_sizes))
}

hist(largest_lineage, breaks=100, xlim=c(0,40), probability = F, main="Size of largest transmission lineage",xlab="Size", ylab="Probability")
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />
