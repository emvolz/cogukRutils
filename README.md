# Scripts for analysis of COGUK sequence data 

## Requirements

- Development versions of 
	- `treedater` >=0.5.1 [https://github.com/emvolz/treedater](https://github.com/emvolz/treedater)
		- Install using `devtools::install_github( 'emvolz/treedater' )`
	- `skygrowth` >=0.3.1 [https://github.com/mrc-ide/skygrowth](https://github.com/mrc-ide/skygrowth)
		- `devtools::install_github( 'mrc-ide/skygrowth' )`
	- `treestructure` >= 0.1.1 [https://github.com/emvolz-phylodynamics/treestructure](https://github.com/emvolz-phylodynamics/treestructure)
		- `devtools::install_github( 'emvolz-phylodynamics/treestructure' )`
- MAFFT 
- IQTree
- `sarscov2` package: [https://github.com/emvolz-phylodynamics/sarscov2Rutils](https://github.com/emvolz-phylodynamics/sarscov2Rutils)
	- This contains miscelaneous tools for SARS CoV 2 developed by the phylo team at Imperial College 
	- Dependencies which should be installed automatically: ape,lubridate,ggtree,ggplot2,treeio,knitr,coda,phangorn,Hmisc,yaml,glue,seqinr,treestructure. 

These functions will make system calls from R which will work on linux and OSX. Windows users will need to use virtualization or alter the internals of several functions. 



## Spike 614 phylodynamics

The script [r/spike0.R](r/spike0.R) contains an example for for inferring growth rates of spike 614 D & G clades using subsets of data. 
This requires a ML tree, metadata in correct format, and and aligment with names matching the tree. 

```{r}
source('r/spike0.R') # load the script 

# Compute a maximum likelihood tree and read into R 
mltree = read.tree( 'mltree.newick' )

# load the coguk metadata 
# NOTE this table must include a column 'tip.label' which matches labels in the mltree 
maj = read.table( 'majora.20200515.metadata.tsv', header=TRUE, stringsAs=FALSE , sep = '\t')

# run the analysis for Edinburgh
result = spike0( maj = maj
 , tr1 = mltree
 , cogfasta = 'cog_2020-05-15_all_alignment.fasta'
 , r = 'EDINBURGH'
) 

# plot effective size 
pne = with( result$spike, add_ggNe_sarscov2skygrowth( sgD, sgG ) )
# plot growth rates
pgr = with( result$spike, add_ggGR_sarscov2skygrowth( sgD, sgG ) )
# plot sample times 
spl = plot_sample_distribution_spike614( result$s614 )

print( pne )
```

The `spike0` function will do the following

- Extract a section of metadata where `adm2` matches the requested region
- Prune the tree to only include matching samples 
- Estimate multiple time trees: 
	- Polytomies in the ML tree are resolved randomly n times
	- Each tree is passed to treedater which estimates a time tree with a strict clock and a rate bounded 0.0008-0.0015 subst / site / year 
	- Each of these is passed through a function that samples node times under the constraint of the ML topology and clock rates estimated in treedater. This is to 'smooth out' the distribution of node times since treedater will often return trees with many short branch lengths corresponding to polytomies in the ML tree.
- The genotype for each sample in the trees are computed from the given alignment. Make sure this aligment has labels which match the tree. 
- Each tree is broken into two clades by pruning out G or D tips. These trees might share nodes, but those generally occur before March 1. 
- For each time tree, Ne(t) and growthrates are estimated using skygrowth
	- The final results are averaged across skygrowth fits to each time tree 
	- This only uses nodes in the tree __after March 1__ by default. This date can be customized with `tstart` option. The G & D trees do not generally share any nodes after that date 
	- `tau0` sets the scale of an exponential prior on the 'smoothing parameter'. Larger values give smoother trajectories and for very large values will force exponential growth. The default value implies that there is 50% probability that growth rates will change less than 4.5 units (1/year) over the course of 1 day. For perspective, a growth rate of 50 / year corresponds to a doubling time of 5.06 days, and reducing the rate to 50-5 = 45 changes the doubling time to 5.6 days. 


To see numerical results see the data frames `result$sgD$growthrate` and `result$sgG$growthrate`
